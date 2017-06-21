#ifndef HMC_TEST
#define HMC_TEST
#include "../include/hmc.hpp"
#include "../include/mklrand.hpp"
#include <cmath>

TEST( hmc, acceptance_always )
{
    // Create rng
    int seed = 1001;
    mkl_drand rng( 1, seed );

    // Always accept when the new energy is lower
    double before = 0.2;
    double after = -0.15;

    bool accept = hmc::accept_energy( before, after, rng );
    EXPECT_TRUE( accept );

    after = 0.0;
    accept = hmc::accept_energy( before, after, rng );
    EXPECT_TRUE( accept );

    after = 0.1999;
    accept = hmc::accept_energy( before, after, rng );
    EXPECT_TRUE( accept );
}

TEST( hmc, accept_probability )
{
    // Create rng
    int seed = 1001;
    mkl_drand rng( 1000, seed );

    // Accept with prob exp(-(after-before))
    double before = 2;
    double after = 2.2;

    double accept_rate = 0;
    size_t trials = 100000;

    for( unsigned int i=0; i<trials; i++ )
        accept_rate += (
            (double) hmc::accept_energy( before, after, rng )
            ) / trials;
    EXPECT_NEAR( std::exp( before - after), accept_rate, accept_rate * 1e-3 );

    accept_rate = 0;
    after = 2.2;
    trials = 5000000;
    for( unsigned int i=0; i<trials; i++ )
        accept_rate += (
            (double) hmc::accept_energy( before, after, rng )
            ) / trials;
    EXPECT_NEAR( std::exp( before - after), accept_rate, accept_rate * 1e-3 );
}
#endif
