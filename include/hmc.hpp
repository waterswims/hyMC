#ifndef HMC_H
#define HMC_H
#include "./mklrand.hpp"
#include <functional>
#include <valarray>
#include <vector>

namespace hmc {

    template <typename T>
    void _swap_ptrs( T* &a, T* &b )
    {
        T *tmp = a;
        a = b;
        b = tmp;
    }

    bool accept_trial(
        const double current_energy,
        const double trial_energy,
        mkl_drand &rng );

    double kinetic_energy(
        const std::valarray<double> velocity );

    std::vector<std::valarray<double> > hmc(
        std::valarray<double> &sample_energy,
        const std::valarray<double> &initial_state,
        const double leapfrog_eps,
        const size_t leapfrog_steps,
        const size_t samples,
        const std::function<double(const std::valarray<double>&)> &f_energy,
        const std::function<void(std::valarray<double>&, const std::valarray<double>&)> &f_energy_grad,
        const std::function<std::valarray<double>(const std::valarray<double>&)> &reduce );

}

#endif
