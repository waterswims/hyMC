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
        energy[sample] = current_energy;

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

    // initialise vector of results
    std::vector<std::valarray<double> > trace;

    // allocate memory for work
    std::valarray<double> current_state( initial_state );
    std::valarray<double> current_velocity( system_size );
    std::valarray<double> temp_state( system_size );
    std::valarray<double> temp_velocity( system_size );
    std::valarray<double> work( system_size );
    std::vector<std::valarray<double> > state_tree(100, std::valarray<double>(system_size));
    std::vector<std::valarray<double> > velocity_tree(100, std::valarray<double>(system_size));
    std::vector<std::valarray<double> > fb_state(2, std::valarray<double>(system_size));
    std::vector<std::valarray<double> > fb_velocity(2, std::valarray<double>(system_size));
    fb_state[0] = current_state;
    fb_state[1] = current_state;
    fb_velocity[0] = current_velocity;
    fb_velocity[1] = current_velocity;

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
            current_velocity[i] = normal_rng.gen();
        fb_velocity[0] = current_velocity;
        fb_velocity[1] = current_velocity;

        // Init tree size
        int tree_size = 1;
        int tree_count = 1;
        state_tree[0] = current_state;
        velocity_tree[0] = current_velocity;
        int tree_height = 1;

        // Find acceptance slice
        double current_energy = total_energy( current_state, current_velocity );
        int u = uniform_rng.gen() * std::exp( -current_energy );

        bool check1 = true, check2 = true;
        // Begin building tree
        while(check1 && check2)
        {
            // Pick random direction
            int dir_choice = int_rng.gen();
            int other_choice = (dir_choice+1)%2;
            temp_state = fb_state[dir_choice];
            temp_velocity = fb_state[dir_choice];
            int dir = dir_choice*2 - 1;

            // Run leapfrog in that direction
            for (int i=0; i < tree_height; i++)
            {
                leapfrog::lfs( state_tree[tree_count],
                               velocity_tree[tree_count], work,
                               temp_state, temp_velocity,
                               f_energy_grad,
                               dir*leapfrog_eps );

                // Update arrays
                temp_state = state_tree[tree_count];
                temp_velocity = velocity_tree[tree_count];
                tree_count++;

                // Check break
                if(i == tree_height - 1) {continue;}
                check1 *= -total_energy(temp_state, temp_velocity) > u - 1000;
                check2 *= dir * ((temp_state - fb_state[other_choice]) *
                         temp_velocity).sum() >= 0;
                check2 *= dir * ((temp_state - fb_state[other_choice]) *
                         fb_velocity[other_choice]).sum() >= 0;
            }

            if (check1 && check2)
            {
                // Add to tree size
                tree_size += tree_height;
                tree_height *= 2;

                fb_state[dir_choice] = temp_state;
                fb_velocity[dir_choice] = temp_velocity;

                check1 *= -total_energy(temp_state, temp_velocity) > u - 1000;
                check2 *= ((fb_state[0] - fb_state[1]) * fb_velocity[0]).sum() >= 0;
                check2 *= ((fb_state[0] - fb_state[1]) * fb_velocity[1]).sum() >= 0;
            }
        }

        // Randomly Sample from the tree
        int choice = int(uniform_rng.gen()*tree_size);
        current_state = state_tree[choice];
        current_velocity = velocity_tree[choice];
        fb_state[0] = current_state;
        fb_state[1] = current_state;

        // Store the energy
        sample_energy[sample] = f_energy(current_state);

        // Store the reduced parameters
        trace.push_back( reduce( current_state ) );
    }

    return trace;
}

/// Interface to Heisenberg hmc
void hmc::heisenberg_model(
    std::valarray<double> &sample_energy,
    std::valarray<double> &sample_magnetisation,
    const std::vector<size_t> system_dimensions,
    const HamiltonianOptions options,
    const double beta,
    const double leapfrog_eps,
    const size_t leapfrog_steps,
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
    auto trace = hmc::hmc( sample_energy, initial_state, leapfrog_eps, leapfrog_steps, nsamples,
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
