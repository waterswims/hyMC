#ifndef HMC_H
#define HMC_H
#include "./mklrand.hpp"
namespace hmc {

    bool accept_energy(
        const double current_energy,
        const double trial_energy,
        mkl_drand &rng );

}

#endif
