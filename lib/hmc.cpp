#include "../include/hmc.hpp"
#include "../include/leapfrog.hpp"
#include "../include/mklrand.hpp"
#include "../include/array_alloc.hpp"
#include <cmath>

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
