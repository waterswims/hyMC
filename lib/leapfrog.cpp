#include "../include/leapfrog.hpp"

void leapfrog::half_step_velocity(
    double *new_vels,
    const double *vels,
    const double *energy_grads,
    const size_t N,
    const double eps )
{
    #pragma simd
    for( unsigned int n=0; n<N; n++ )
        new_vels[n] = vels[n] - eps/2.0 * energy_grads[n];
}

void leapfrog::step_state(
    double *new_state,
    const double *state,
    const double *vels,
    const size_t N,
    const double eps )
{
    #pragma simd
    for( unsigned int n=0; n<N; n++ )
        new_state[n] = state[n] + eps * vels[n];
}

void leapfrog::lfs(
    double *new_state,
    double *new_velocity,
    double *energy_grads_work,
    const double *state,
    const double *velocity,
    const std::function<void(double*,const double*,const size_t)> &energy_grads,
    const size_t N,
    const double eps )
{
    // Compute the gradient of the energy landscape
    energy_grads( energy_grads_work, state, N );

    // Compute the half-step velocity and next state
    half_step_velocity( new_velocity, velocity, energy_grads_work, N, eps );
    step_state( new_state, state, new_velocity, N, eps );

    // Finish by computing the full-step velocity
    energy_grads( energy_grads_work, new_state, N );
    half_step_velocity( new_velocity, new_velocity, energy_grads_work, N, eps );
}
