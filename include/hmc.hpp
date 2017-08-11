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

    std::vector<std::valarray<double> > nuts(
        std::valarray<double> &sample_energy,
        const std::valarray<double> &initial_state,
        const double leapfrog_eps,
        const size_t samples,
        const std::function<double(const std::valarray<double>&)> &f_energy,
        const std::function<void(std::valarray<double>&, const std::valarray<double>&)> &f_energy_grad,
        const std::function<std::valarray<double>(const std::valarray<double>&)> &reduce
    );

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Builds the sample tree for the NUTS.
    ///
    /// \param back_state Input state
    /// \param back_vel Input Velocity
    /// \param back_state State at the back of the tree
    /// \param back_vel Velocity at the back of the tree
    /// \param for_state State at the front of the tree
    /// \param for_vel Velocity at the front of the tree
    /// \param slice Random number which decides the accepted energies
    /// \param dir_choice The choice of direction
    /// \param tree_height The current height of the energy tree
    /// \param eps Time step for the leapfrog integrator
    /// \param state_tree The constructed tree of states
    /// \param vel_tree The constructed tree of velocities
    /// \param tree_size The size of the constructed tree
    /// \param break_check Check for the completion of the tree
    /// \param energy_grads A function which calculates the energy gradient
    /// \param work Valarray to work in
    ///////////////////////////////////////////////////////////////////////////
    void build_tree(
        std::valarray<double> &in_state,
        std::valarray<double> &in_vel,
        std::valarray<double> &back_state,
        std::valarray<double> &back_vel,
        std::valarray<double> &for_state,
        std::valarray<double> &for_vel,
        const double slice,
        const int dir_choice,
        const int tree_height,
        const double eps,
        std::vector<std::valarray<double> > &state_tree,
        std::vector<std::valarray<double> > &vel_tree,
        int &tree_count,
        int &tree_added,
        bool &break_check,
        const std::function<void(std::valarray<double>&,const std::valarray<double>&)> &energy_grads,
        const std::function<double(const std::valarray<double>&, const std::valarray<double>&)> &energy,
        std::valarray<double> &work
    );

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

    double magnetisation( const std::valarray<double>& state );
}

#endif
