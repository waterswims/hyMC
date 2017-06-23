#include "../include/leapfrog.hpp"

void leapfrog::half_step_velocity(
    std::valarray<double> &new_vels,
    const std::valarray<double> &vels,
    const std::valarray<double> &energy_grads,
    const double eps )
{
    new_vels = vels - eps / 2.0 * energy_grads;
}

void leapfrog::step_state(
    std::valarray<double> &new_state,
    const std::valarray<double> &state,
    const std::valarray<double> &vels,
    const double eps )
{
    new_state = state + eps * vels;
}

void leapfrog::lfs(
    std::valarray<double> &new_state,
    std::valarray<double> &new_velocity,
    std::valarray<double> &energy_grads_work,
    const std::valarray<double> &state,
    const std::valarray<double> &velocity,
    const std::function<void(std::valarray<double>&,const std::valarray<double>&)> &energy_grads,
    const double eps )
{
    // Compute the gradient of the energy landscape
    energy_grads( energy_grads_work, state);

    // Compute the half-step velocity and next state
    half_step_velocity( new_velocity, velocity, energy_grads_work, eps );
    step_state( new_state, state, new_velocity, eps );

    // Finish by computing the full-step velocity
    energy_grads( energy_grads_work, new_state );
    half_step_velocity( new_velocity, new_velocity, energy_grads_work, eps );
}
