#include "../include/hmc.hpp"
#include <cmath>

bool hmc::accept_energy( const double e, const double e_trial,
                         mkl_drand &rng )
{
    double e_diff = e_trial - e;
    double acceptance_threshold = std::exp( -e_diff );
    return ( rng.gen() < acceptance_threshold );
}
