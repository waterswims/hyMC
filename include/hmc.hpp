#ifndef HMC_H
#define HMC_H
#include "./mklrand.hpp"
#include <functional>
#include <valarray>
#include <vector>
#include <unordered_map>

namespace hmc {

    struct HamiltonianOptions {
        double J;
        double H;
    };

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
        mklrand::mkl_drand &rng );

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

    void heisenberg_model(
        std::valarray<double> &sample_energy,
        std::valarray<double> &sample_magnetisation,
        const std::vector<size_t> system_dimensions,
        const HamiltonianOptions options,
        const double beta,
        const double leapfrog_eps,
        const size_t leapfrog_steps,
        const size_t nsamples,
        const long initial_state_seed );

    std::function<void(std::valarray<double>&, const std::valarray<double>&)>
    heisenberg_gradients( const size_t ndims, const HamiltonianOptions options );

    std::function<double(const std::valarray<double>&)>
    heisenberg_hamiltonian( const size_t ndims, const HamiltonianOptions options, const double beta );

    double magnetisation( const std::valarray<double>& state );
}

#endif
