#include "../include/data_types.hpp"

hmc::hberg_lattice::hberg_lattice(const size_t length)
{
    N = length;
    thetas = std::valarray<double>(length);
    phis = std::valarray<double>(length);
}

double hmc::hberg_lattice::kinetic_energy()
{
    return (pow(thetas, 2.) + pow(phis, 2.)).sum() / 2.0;
}

double hmc::hberg_lattice_1d::total_energy(d_type* vels)
{
    double E = vels->kinetic_energy();

    if(E_flags[0])
    {
        E += exchange_energy();
    }

    return E;
}

void hmc::hberg_lattice_1d::total_energy_grad(d_type* grad_out)
{
    grad_out->zero();

    if(E_flags[0])
    {
        exchange_energy_grad(grad_out);
    }
}

double hmc::hberg_lattice_1d::exchange_energy()
{
    std::valarray<double> cos_theta = cos(thetas);
    std::valarray<double> sin_theta = sin(thetas);
    std::valarray<double> cos_phi = cos(phis);
    std::valarray<double> sin_phi = sin(phis);

    std::valarray<double> cos_theta_left = cos_theta.cshift(1);
    std::valarray<double> sin_theta_left = sin_theta.cshift(1);
    std::valarray<double> cos_phi_left = cos_phi.cshift(1);
    std::valarray<double> sin_phi_left = sin_phi.cshift(1);

    double t1 = (cos_phi * cos_phi_left).sum();
    double t2 = (cos_theta * cos_theta_left * sin_phi * sin_phi_left).sum();
    double t3 = (sin_theta * sin_phi * sin_theta_left * sin_phi_left).sum();

    return -(t1 + t2 + t3);
}

void hmc::hberg_lattice_1d::exchange_energy_grad(d_type* grad_out)
{
    std::valarray<double> cos_theta = cos(thetas);
    std::valarray<double> sin_theta = sin(thetas);
    std::valarray<double> cos_phi = cos(phis);
    std::valarray<double> sin_phi = sin(phis);

    std::valarray<double> cos_theta_left = cos_theta.cshift(1);
    std::valarray<double> sin_theta_left = sin_theta.cshift(1);
    std::valarray<double> cos_phi_left = cos_phi.cshift(1);
    std::valarray<double> sin_phi_left = sin_phi.cshift(1);

    std::valarray<double> cos_theta_right = cos_theta.cshift(-1);
    std::valarray<double> sin_theta_right = sin_theta.cshift(-1);
    std::valarray<double> cos_phi_right = cos_phi.cshift(-1);
    std::valarray<double> sin_phi_right = sin_phi.cshift(-1);

    std::valarray<double> t1, t2, t3, t4, t5, t6;
    t1 = -cos_theta_left * sin_theta * sin_phi * sin_phi_left;
    t2 = cos_theta * sin_phi * sin_theta_left * sin_phi_left;
    t3 = -cos_theta_right * sin_theta * sin_phi * sin_phi_right;
    t4 = cos_theta * sin_phi * sin_theta_right * sin_phi_right;
    grad_out->add_to_d1(t1 + t2 + t3 + t4);

    t1 = -cos_phi_left * sin_phi;
    t2 = cos_theta * cos_phi * cos_theta_left * sin_phi_left;
    t3 = cos_phi * sin_theta * sin_theta_left * sin_phi_left;
    t4 = -cos_phi_right * sin_phi;
    t5 = cos_theta * cos_phi * cos_theta_right * sin_phi_right;
    t6 = cos_phi * sin_theta * sin_theta_right * sin_phi_right;
    grad_out->add_to_d2(t1 + t2 + t3 + t4 + t5 + t6);
}
