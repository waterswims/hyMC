#include "../include/all_hamils.hpp"

double hmc::zeeman_energy(
    const std::vector<std::valarray<double> >& trig_angles, double H)
{
    int ind = trig_angles.size() / 2;
    return (H * trig_angles[ind]).sum();
}
