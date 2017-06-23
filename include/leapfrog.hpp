#ifndef LEAPFROG_H
#define LEAPFROG_H
#include <functional>
#include <valarray>

namespace leapfrog {

    /// New velocity of the system states after half a step
    void half_step_velocity(
        std::valarray<double> &new_vels,
        const std::valarray<double> &vels,
        const std::valarray<double> &energy_gads,
        const double eps );

    /// New states of the system after a full step
    void step_state( std::valarray<double> &new_state,
                     const std::valarray<double> &state,
                     const std::valarray<double> &vels,
                     const double eps );

    /// New state and velocity of the system after one leap frog step
    void lfs(
        std::valarray<double> &new_state,
        std::valarray<double> &new_velocity,
        std::valarray<double> &energy_grads_work,
        const std::valarray<double> &state,
        const std::valarray<double> &velocity,
        const std::function<void(std::valarray<double>&,const std::valarray<double> &)> &energy_grads,
        const double eps );

} // end namespace

#endif
