#include "../include/hmc.hpp"
#include "../include/all_hamils.hpp"
#include "../include/leapfrog.hpp"
#include "../include/mklrand.hpp"
#include "../include/constants.hpp"
#include <cmath>
#include <exception>
#include <iostream>
#define _USE_MATH_DEFINES

bool hmc::accept_trial( const double e, const double e_trial,
                         mklrand::mkl_drand &rng )
{
    double e_diff = e_trial - e;
    double acceptance_threshold = std::exp( -e_diff );
    return ( rng.gen() < acceptance_threshold );
}

double hmc::kinetic_energy(
    const std::valarray<double> velocity )
{
    double energy=0;
    for( auto &v : velocity )
        energy += v * v;
    return energy/2.0;
}

std::vector<std::valarray<double> > hmc::hmc(
    std::valarray<double> &energy,
    const std::valarray<double> &initial_state,
    const double leapfrog_eps,
    const size_t leapfrog_steps,
    const size_t samples,
    const std::function<double(const std::valarray<double>&)> &f_energy,
    const std::function<void(std::valarray<double>&, const std::valarray<double>&)> &f_energy_grad,
    const std::function<std::valarray<double>(const std::valarray<double>&)> &reduce )
{
    // system size
    size_t system_size = initial_state.size();

    // initialise vector of results
    std::vector<std::valarray<double> > trace;

    // allocate memory for work
    std::valarray<double> current_state( initial_state );
    std::valarray<double> current_velocity( system_size );
    std::valarray<double> temp_state( system_size );
    std::valarray<double> temp_velocity( system_size );
    std::valarray<double> trial_state( system_size );
    std::valarray<double> trial_velocity( system_size );
    std::valarray<double> work( system_size );

    // Allocate RNGs
    mklrand::mkl_drand uniform_rng( 100000, 1001 );
    mklrand::mkl_nrand normal_rng( 0, 1, 100000, 555555 );

    // Total energy
    std::function<double(const std::valarray<double>&, const std::valarray<double>&)>
        total_energy = [&f_energy](const std::valarray<double>& state, const std::valarray<double>& velocities)
        { return f_energy( state ) + kinetic_energy( velocities ); };

    // Run a monte carlo step until we get specific number of samples
    for( unsigned int sample=0; sample<samples; sample++ )
    {
        // Get the initial state and a random choice of velocity
        temp_state = current_state;
        for( unsigned int i=0; i<system_size; i++ )
            temp_velocity[i] = normal_rng.gen();
        current_velocity = temp_velocity;


        // Run a number of leapfrog steps
        for( unsigned int n=0; n<leapfrog_steps; n++ )
        {
            leapfrog::lfs( trial_state, trial_velocity, work,
                           temp_state, temp_velocity,
                           f_energy_grad,
                           leapfrog_eps );

            // Update arrays
            temp_state = trial_state;
            temp_velocity = trial_velocity;
        }

        // Compute energies
        double current_energy = total_energy( current_state, current_velocity );
        double trial_energy = total_energy( trial_state, trial_velocity );

        // check acceptance
        bool accept = accept_trial( current_energy, trial_energy, uniform_rng );
        if( accept )
        {
            current_state = trial_state;
            current_energy = trial_energy;
        }


        // Store the energy
        energy[sample] = f_energy(current_state);

        // Store the reduced parameters
        trace.push_back( reduce( current_state ) );
    }

    return trace;
}

std::vector<std::valarray<double> > hmc::nuts(
    std::valarray<double> &sample_energy,
    const std::valarray<double> &initial_state,
    const double leapfrog_eps,
    const size_t samples,
    const std::function<double(const std::valarray<double>&)> &f_energy,
    const std::function<void(std::valarray<double>&, const std::valarray<double>&)> &f_energy_grad,
    const std::function<std::valarray<double>(const std::valarray<double>&)> &reduce
)
{
    // system size
    size_t system_size = initial_state.size();

    // max tree height
    size_t max_tree_height = 15;

    // initialise vector of results
    std::vector<std::valarray<double> > trace;

    // allocate memory for work
    std::valarray<double> current_state( initial_state );
    std::valarray<double> current_velocity( system_size );
    std::valarray<double> work( system_size );
    std::valarray<double> temp_state(system_size);
    std::valarray<double> next_state(system_size);
    std::vector<std::valarray<double> > fb_state(2, std::valarray<double>(system_size));
    std::vector<std::valarray<double> > fb_velocity(2, std::valarray<double>(system_size));
    std::vector<std::valarray<double> > dud_state_tree(max_tree_height, std::valarray<double>(system_size));
    std::vector<std::valarray<double> > dud_vel_tree(max_tree_height, std::valarray<double>(system_size));
    std::vector<std::valarray<double> > pos_state_tree(max_tree_height, std::valarray<double>(system_size));

    // Allocate RNGs
    mklrand::mkl_irand int_rng(100000, 666);
    mklrand::mkl_drand uniform_rng( 100000, 1001 );
    mklrand::mkl_nrand normal_rng( 0, 1, 100000, 555555 );

    // Total energy
    std::function<double(const std::valarray<double>&, const std::valarray<double>&)>
        total_energy = [&f_energy](const std::valarray<double>& state, const std::valarray<double>& velocities)
        { return f_energy( state ) + kinetic_energy( velocities ); };

    // Run a monte carlo step until we get specific number of samples
    for( unsigned int sample=0; sample<samples; sample++ )
    {
        // Get the initial state and a random choice of velocity
        for( unsigned int i=0; i<system_size; i++ )
            {current_velocity[i] = normal_rng.gen();}
        fb_velocity[0] = current_velocity;
        fb_velocity[1] = current_velocity;
        fb_state[0] = current_state;
        fb_state[1] = current_state;
        next_state = current_state;

        // Init tree size
        int n = 1;
        int tree_height = 0;

        // Find acceptance slice
        double current_energy = total_energy( current_state, current_velocity );
        double r = uniform_rng.gen();
        double lu = std::log(r) - current_energy;

        bool check1 = true;
        // Begin building tree
        while(check1)
        {
            // Init n for this tree
            int temp_n = 0;

            // Pick random direction
            int dir_choice = int_rng.gen();
            int dir = 2*dir_choice - 1;
            double temp_eps = leapfrog_eps*dir;

            // Run leapfrog in that direction
            if(temp_eps > 0)
            {
                build_tree(fb_state[1], fb_velocity[1], dud_state_tree[tree_height+1], dud_vel_tree[tree_height+1], fb_state[1], fb_velocity[1],
                           temp_state, dud_state_tree, dud_vel_tree, pos_state_tree, lu, tree_height, temp_eps, check1, temp_n, f_energy_grad, total_energy, work, uniform_rng);
            }
            else
            {
                build_tree(fb_state[0], fb_velocity[0], fb_state[0], fb_velocity[0], dud_state_tree[tree_height+1], dud_vel_tree[tree_height+1],
                           temp_state, dud_state_tree, dud_vel_tree, pos_state_tree, lu, tree_height, temp_eps, check1, temp_n, f_energy_grad, total_energy, work, uniform_rng);
            }

            if (check1 && uniform_rng.gen() < temp_n/float(n)) {next_state = temp_state;}
            n += temp_n;
            check1 *= (((fb_state[1] - fb_state[0]) * fb_velocity[0]).sum() >= 0);
            check1 *= (((fb_state[1] - fb_state[0]) * fb_velocity[1]).sum() >= 0);
            tree_height++;
        }

        // Set next state
        current_state = next_state;

        // Store the energy
        sample_energy[sample] = f_energy(current_state);

        // Store the reduced parameters
        trace.push_back( reduce( current_state ) );
    }
    return trace;
}

void hmc::build_tree(
    std::valarray<double> &in_state,
    std::valarray<double> &in_vel,
    std::valarray<double> &back_state,
    std::valarray<double> &back_vel,
    std::valarray<double> &for_state,
    std::valarray<double> &for_vel,
    std::valarray<double> &out_state,
    std::vector<std::valarray<double> > &dud_state_tree,
    std::vector<std::valarray<double> > &dud_vel_tree,
    std::vector<std::valarray<double> > &pos_state_tree,
    const double slice,
    const int tree_height,
    const double eps,
    bool &break_check,
    int &n_check,
    const std::function<void(std::valarray<double>&,const std::valarray<double>&)> &energy_grads,
    const std::function<double(const std::valarray<double>&, const std::valarray<double>&)> &energy,
    std::valarray<double> &work,
    mklrand::mkl_drand &rng
)
{
    if(tree_height == 0)
    {
        leapfrog::lfs(dud_state_tree[tree_height], dud_vel_tree[tree_height], work, in_state, in_vel, energy_grads, eps);

        back_state = dud_state_tree[tree_height];
        back_vel = dud_vel_tree[tree_height];
        for_state = dud_state_tree[tree_height];
        for_vel = dud_vel_tree[tree_height];
        out_state = dud_state_tree[tree_height];

        double temp_E = -energy(dud_state_tree[tree_height], dud_vel_tree[tree_height]);
        n_check = (temp_E >= slice);
        break_check *= (temp_E >= slice - 1000);
        return;
    }
    else
    {
        n_check = 0;
        build_tree(in_state, in_vel, back_state, back_vel, for_state, for_vel,
                   out_state, dud_state_tree, dud_vel_tree, pos_state_tree, slice, tree_height-1, eps, break_check, n_check,
                   energy_grads, energy, work, rng);
        if(break_check)
        {
            int temp_n = 0;
            if(eps > 0)
            {
                build_tree(for_state, for_vel, dud_state_tree[tree_height], dud_vel_tree[tree_height], for_state,
                           for_vel, pos_state_tree[tree_height-1], dud_state_tree, dud_vel_tree, pos_state_tree, slice, tree_height-1, eps,
                           break_check, temp_n, energy_grads, energy, work, rng);
            }
            else
            {
                build_tree(back_state, back_vel, back_state, back_vel, dud_state_tree[tree_height],
                           dud_vel_tree[tree_height], pos_state_tree[tree_height-1], dud_state_tree, dud_vel_tree, pos_state_tree, slice, tree_height-1, eps,
                           break_check, temp_n, energy_grads, energy, work, rng);
            }
            if(rng.gen() < temp_n/float(temp_n+n_check)) {out_state = pos_state_tree[tree_height-1];}
            break_check *= (((for_state - back_state) * back_vel).sum() >= 0);
            break_check *= (((for_state - back_state) * for_vel).sum() >= 0);
            n_check += temp_n;
        }
        return;
    }
}

/// Interface to Heisenberg hmc
void hmc::heisenberg_model(
    std::valarray<double> &sample_energy,
    std::valarray<double> &sample_magnetisation,
    const std::vector<size_t> system_dimensions,
    const HamiltonianOptions options,
    const double beta,
    const double leapfrog_eps,
    const size_t nsamples,
    const long initial_state_seed )
{
    // Compute the size of the state vector
    // theta and phi for every element in system
    const double state_size = 2 * std::accumulate(
        system_dimensions.begin(),
        system_dimensions.end(),
        1, std::multiplies<double>() );

    // Random initial state is controlled with the initial_state_seed
    std::valarray<double> initial_state( state_size );
    mklrand::mkl_drand rng( initial_state_seed );
    for( unsigned int i=0; i<state_size; i++ )
        initial_state[i] = rng.gen() * 2 * M_PI;

    // Get the Hamiltonian and gradient functions
    size_t ndim = system_dimensions.size();
    auto energy_function = gen_total_energy( options, beta, ndim );
    auto grad_function = gen_total_grad( options, ndim );

    // Reduction to compute the magnetisation
    std::function<std::valarray<double>(const std::valarray<double>&)>
        reduce = []( const std::valarray<double>& state )
        {
            std::valarray<double> res = { magnetisation( state ) };
            return res;
        };

    // EXECUTE HMC
    auto trace = hmc::nuts( sample_energy, initial_state, leapfrog_eps, nsamples,
                          energy_function, grad_function, reduce );

    // Energy is returned normalised so we turn it into real energy
    sample_energy = sample_energy / beta;

    // Magnetisation is stored in the trace
    for( size_t n=0; n<nsamples; n++ )
        sample_magnetisation[n] = trace[n][0];
}

double hmc::magnetisation( const std::valarray<double>& state )
{
    int halfsize = state.size() / 2;
    std::slice tslice(0, halfsize, 1);
    std::slice pslice(halfsize, halfsize, 1);
    double x = ( cos( state[tslice] ) * sin( state[pslice] ) ).sum();
    double y = ( sin( state[tslice] ) * sin( state[pslice] ) ).sum();
    double z = cos( state[pslice] ).sum();
    return std::sqrt( x*x + y*y + z*z );
}
