#include "../include/all_hamils.hpp"
#include <iostream>

double hmc::zeeman_energy(
    const std::vector<std::valarray<double> >& trig_angles, double H)
{
    int ind = trig_angles.size() / 2;
    return -(H * trig_angles[ind]).sum();
}

void hmc::zeeman_grad(std::valarray<double>& grad_out,
    const std::vector<std::valarray<double> >& trig_angles, double H)
{
    int halfsize = trig_angles[0].size();
    std::slice pslice(halfsize, halfsize, 1);
    int d = trig_angles.size() / 4;

    grad_out[pslice] += H * trig_angles[d*3];
}
