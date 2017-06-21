#ifndef LEAPFROG_H
#define LEAPFROG_H
#include <functional>

namespace leapfrog {

    /// New velocity of the system states after half a step
    void half_step_velocity(
        double *new_vels,
        const double *vels,
        const double *energy_gads,
        const size_t system_size,
        const double eps );

    /// New states of the system after a full step
    void step_state( double *new_state,
                     const double *state,
                     const double *vels,
                     const size_t system_size,
                     const double eps );

    /// New state and velocity of the system after one leap frog step
    void lfs(
        double *new_state,
        double *new_velocity,
        double *energy_grads_work,
        const double *state,
        const double *velocity,
        const std::function<void(double*,const double*, const size_t)> &energy_grads,
        const size_t system_size,
        const double eps );

} // end namespace

#endif
