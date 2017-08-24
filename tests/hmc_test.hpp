#ifndef HMC_TEST
#define HMC_TEST

#include "../include/hmc.hpp"
#include "../include/mklrand.hpp"
#include "test_funcs.hpp"
#include <functional>
#include <cmath>
#include <valarray>
#include <iostream>

TEST( hmc, kinetic_energy )
{
    std::valarray<double> vels = { 1.5, 2.5, -4.0 };
    EXPECT_DOUBLE_EQ( 12.25, hmc::kinetic_energy( vels ) );
}

TEST( hmc, acceptance_always )
{
    // Create rng
    int seed = 1001;
    mklrand::mkl_drand rng( 1, seed );

    // Always accept when the new energy is lower
    double before = 0.2;
    double after = -0.15;

    bool accept = hmc::accept_trial( before, after, rng );
    EXPECT_TRUE( accept );

    after = 0.0;
    accept = hmc::accept_trial( before, after, rng );
    EXPECT_TRUE( accept );

    after = 0.1999;
    accept = hmc::accept_trial( before, after, rng );
    EXPECT_TRUE( accept );
}

TEST( hmc, accept_probability )
{
    // Create rng
    int seed = 1001;
    mklrand::mkl_drand rng( 1000, seed );

    // Accept with prob exp(-(after-before))
    double before = 2;
    double after = 2.2;

    double accept_rate = 0;
    size_t trials = 100000;

    for( unsigned int i=0; i<trials; i++ )
        accept_rate += (
            (double) hmc::accept_trial( before, after, rng )
            ) / trials;
    EXPECT_NEAR( std::exp( before - after), accept_rate, accept_rate * 1e-3 );

    accept_rate = 0;
    after = 2.2;
    trials = 5000000;
    for( unsigned int i=0; i<trials; i++ )
        accept_rate += (
            (double) hmc::accept_trial( before, after, rng )
            ) / trials;
    EXPECT_NEAR( std::exp( before - after), accept_rate, accept_rate * 1e-3 );
}

TEST( hmc, swap_pts )
{
    double *x = (double*) malloc(2*sizeof(double));
    x[0]=1; x[1]=2;
    double *y = (double*) malloc(2*sizeof(double));
    y[0]=3; y[1]=4;

    hmc::_swap_ptrs( x, y );
    EXPECT_DOUBLE_EQ( 1, y[0] );
    EXPECT_DOUBLE_EQ( 2, y[1] );
    EXPECT_DOUBLE_EQ( 3, x[0] );
    EXPECT_DOUBLE_EQ( 4, x[1] );

    free(x);
    free(y);
}

TEST( hmc, hmc_sampling_normal_distribution )
{
    // Sampling from a univariate normal distribution

    // initial state
    std::valarray<double> x_init = {1.0};

    // normal parameters
    double mu = 1.0;
    double std = 0.8;

    // energy function
    std::function<double(const std::valarray<double>&)>
        energy = [mu,std]( const std::valarray<double>&x)
        { return (x[0]-mu)*(x[0]-mu)/(2*std*std); };

    // energy gradient
    std::function<void(std::valarray<double>&,const std::valarray<double>&)>
        energy_grad = [mu,std](std::valarray<double>&grad, const std::valarray<double>&x )
        { grad[0] = (x[0]-mu)/(std*std); };

    // reduce function
    std::function<std::valarray<double>(const std::valarray<double>&)>
        reduce = [](const std::valarray<double>&x)
        { return std::valarray<double>( x ); };

    // Run hmc
    size_t N = 1000000;
    std::valarray<double> energies( N );
    double eps = 0.18;
    size_t steps = 20;

    // Run HMC!
    auto trace = hmc::hmc( energies, x_init, eps, steps, N, energy, energy_grad, reduce );

    // Are the dimensions correct?
    EXPECT_EQ( trace.size(), N );
    EXPECT_EQ( trace[0].size(), 1 );

    // Compute sample mean and variance
    double sample_mean = 0.0;
    for( unsigned int i=N/2.; i<N; i++ )
    {
        sample_mean += trace[i][0];
    }
    sample_mean /= N/2.;

    double expected_squares = 0;
    for( unsigned int i=N/2.; i<N; i++ )
        expected_squares += trace[i][0] * trace[i][0];
    expected_squares /= N/2;
    double sample_variance = expected_squares - sample_mean*sample_mean;
    double sample_std = std::sqrt( sample_variance );

    EXPECT_NEAR( mu, sample_mean, 0.001 );
    EXPECT_NEAR( std, sample_std, 0.0015 );

    int N_bins = 100;
    std::vector<int> bins(N_bins, 0);
    for(int i = 0; i < N; i++)
    {
        double g = trace[i][0];
        if (g >= 3 || g < -1) {continue;}
        int b_num = int((g+1)*(N_bins/4.0));
        bins[b_num]++;
    }
    std::vector<double> expect(N_bins);
    for(int i=0; i < N_bins; i++)
    {
        double lower = (i)*4.0/float(N_bins) - 1;
        double higher = (i+1)*4.0/float(N_bins) - 1;
        double temp = erf((higher - sample_mean)/(sample_std*pow(2, 0.5))) - erf((lower - sample_mean)/(sample_std*pow(2, 0.5)));
        expect[i] = N*0.5*temp;
    }
    double chi2_test = chi2(bins, expect);
    EXPECT_GT(chi2_test, 0.9);
    EXPECT_LT(chi2_test, 1.3);
}

TEST( hmc, nuts_sampling_normal_distribution )
{
    // Sampling from a univariate normal distribution

    // initial state
    std::valarray<double> x_init = {1.0};

    // normal parameters
    double mu = 1.0;
    double std = 0.8;

    // energy function
    std::function<double(const std::valarray<double>&)>
        energy = [mu,std]( const std::valarray<double>&x)
        { return (x[0]-mu)*(x[0]-mu)/(2*std*std); };

    // energy gradient
    std::function<void(std::valarray<double>&,const std::valarray<double>&)>
        energy_grad = [mu,std](std::valarray<double>&grad, const std::valarray<double>&x )
        { grad[0] = (x[0]-mu)/(std*std); };

    // reduce function
    std::function<std::valarray<double>(const std::valarray<double>&)>
        reduce = [](const std::valarray<double>&x)
        { return std::valarray<double>( x ); };

    // Run hmc
    size_t N = 1000000;
    std::valarray<double> energies( N );
    double eps = 0.25;

    // Run HMC!
    auto trace = hmc::nuts( energies, x_init, eps, N, energy, energy_grad, reduce );

    // Are the dimensions correct?
    EXPECT_EQ( trace.size(), N );
    EXPECT_EQ( trace[0].size(), 1 );

    // Compute sample mean and variance
    double sample_mean = 0.0;
    for( unsigned int i=0; i<N; i++ )
    {
        sample_mean += trace[i][0];
    }
    sample_mean /= float(N);

    double expected_squares = 0;
    for( unsigned int i=0; i<N; i++ )
        expected_squares += trace[i][0] * trace[i][0];
    expected_squares /= float(N);
    double sample_variance = expected_squares - sample_mean*sample_mean;
    double sample_std = std::sqrt( sample_variance );

    EXPECT_NEAR( mu, sample_mean, 0.001 );
    EXPECT_NEAR( std, sample_std, 0.0015 );

    int N_bins = 100;
    std::vector<int> bins(N_bins, 0);
    for(int i = 0; i < N; i++)
    {
        double g = trace[i][0];
        if (g > 3 || g < -1) {continue;}
        int b_num = int((g+1)*(N_bins/4.0));
        bins[b_num]++;
    }
    std::vector<double> expect(N_bins);
    for(int i=0; i < N_bins; i++)
    {
        double lower = (i)*4.0/float(N_bins) - 1;
        double higher = (i+1)*4.0/float(N_bins) - 1;
        double temp = erf((higher - sample_mean)/(sample_std*pow(2, 0.5))) - erf((lower - sample_mean)/(sample_std*pow(2, 0.5)));
        expect[i] = N*0.5*temp;
    }
    double chi2_test = chi2(bins, expect);
    EXPECT_GT(chi2_test, 0.9);
    EXPECT_LT(chi2_test, 1.3);
}

TEST( hmc, magnetisaton )
{
    std::valarray<double> state = { 0.2, 0.2, 1.1, 1.1 };
    hmc::set_slices(4, 1);
    EXPECT_DOUBLE_EQ( 2.0, hmc::magnetisation( state ) );
}
#endif
