#include "../include/all_hamils.hpp"
#include <iostream>

namespace hmc
{
    int halfsize;
    std::vector<std::gslice> large_slice, small_slice, opp_large_slice;
    std::vector<std::gslice> opp_small_slice;
    std::slice tslice, pslice;
    std::valarray<double> temp, small_temp, cos_the, sin_the, cos_phi, sin_phi;
}

double hmc::zeeman_energy(
    const double H)
{
    return -H * cos_the.sum();
}

void hmc::zeeman_grad(
    std::valarray<double>& grad_out,
    const double H)
{
    grad_out[tslice] += H * sin_the;
}

double hmc::exchange_energy(
    const double J,
    const int d)
{
    double temp_sum = 0;

    // mult xs
    small_temp = 0;
    for(int i = 0; i < d; i++)
    {
        small_temp[opp_large_slice[i]] += cos_phi[large_slice[i]]*sin_the[large_slice[i]];
        small_temp[opp_small_slice[i]] += cos_phi[small_slice[i]]*sin_the[small_slice[i]];
    }
    temp_sum += (cos_phi*sin_the*small_temp).sum();

    // mult ys
    small_temp = 0;
    for(int i = 0; i < d; i++)
    {
        small_temp[opp_large_slice[i]] += sin_phi[large_slice[i]]*sin_the[large_slice[i]];
        small_temp[opp_small_slice[i]] += sin_phi[small_slice[i]]*sin_the[small_slice[i]];
    }
    temp_sum += (sin_phi*sin_the*small_temp).sum();

    // mult zs
    small_temp = 0;
    for(int i = 0; i < d; i++)
    {
        small_temp[opp_large_slice[i]] += cos_the[large_slice[i]];
        small_temp[opp_small_slice[i]] += cos_the[small_slice[i]];
    }
    temp_sum += (cos_the*small_temp).sum();

    return -J * temp_sum;
}

void hmc::exchange_grad(
    std::valarray<double>& grad_out,
    const double J,
    const int d)
{
    // mult xs
    small_temp = 0;
    for(int i = 0; i < d; i++)
    {
        small_temp[opp_large_slice[i]] += cos_phi[large_slice[i]]*sin_the[large_slice[i]];
        small_temp[opp_small_slice[i]] += cos_phi[small_slice[i]]*sin_the[small_slice[i]];
        small_temp[large_slice[i]] += cos_phi[opp_large_slice[i]]*sin_the[opp_large_slice[i]];
        small_temp[small_slice[i]] += cos_phi[opp_small_slice[i]]*sin_the[opp_small_slice[i]];
    }
    grad_out[pslice] -= J*sin_phi*sin_the*small_temp;
    grad_out[tslice] += J*cos_phi*cos_the*small_temp;

    // mult ys
    small_temp = 0;
    for(int i = 0; i < d; i++)
    {
        small_temp[opp_large_slice[i]] += sin_phi[large_slice[i]]*sin_the[large_slice[i]];
        small_temp[opp_small_slice[i]] += sin_phi[small_slice[i]]*sin_the[small_slice[i]];
        small_temp[large_slice[i]] += sin_phi[opp_large_slice[i]]*sin_the[opp_large_slice[i]];
        small_temp[small_slice[i]] += sin_phi[opp_small_slice[i]]*sin_the[opp_small_slice[i]];
    }
    grad_out[pslice] += J*cos_phi*sin_the*small_temp;
    grad_out[tslice] += J*sin_phi*cos_the*small_temp;

    // mult zs
    small_temp = 0;
    for(int i = 0; i < d; i++)
    {
        small_temp[opp_large_slice[i]] += cos_the[large_slice[i]];
        small_temp[opp_small_slice[i]] += cos_the[small_slice[i]];
        small_temp[large_slice[i]] += cos_the[opp_large_slice[i]];
        small_temp[small_slice[i]] += cos_the[opp_small_slice[i]];
    }
    grad_out[tslice] -= J*sin_the*small_temp;
}

std::function<double(const std::valarray<double>&)> hmc::gen_total_energy(
    const HamiltonianOptions options,
    const double beta,
    const int d,
    const int size)
{
    // Init blank energy function
    std::function<double(const std::valarray<double>&)> init_f =
        [options](const std::valarray<double> &data)
        {
            calc_trig(data);
            return 0;
        };
    std::function<double(const std::valarray<double>&)> new_f;

    // Add exchange energy
    if ( options.J != 0)
    {
        new_f =[init_f, options, d](const std::valarray<double> &data)
            {return init_f(data) + exchange_energy(options.J, d);};
        init_f = new_f;
    }

    // Add zeeman energy
    if ( options.H != 0 )
    {
        new_f = [init_f, options](const std::valarray<double> &data)
            {return init_f(data) + zeeman_energy(options.H);};
        init_f = new_f;
    }

    // Set work arrs and slices
    set_slices(size, d);

    return init_f;
}

std::function<void(std::valarray<double>&, const std::valarray<double>&)>
hmc::gen_total_grad(
    const HamiltonianOptions options,
    const int d,
    const int size)
{
    // Init blank grad function
    std::function<void(std::valarray<double>&, const std::valarray<double>&)>
        init_f, new_f;
    init_f = [options](
        std::valarray<double>& grad_out,
        const std::valarray<double>& data)
        {
            grad_out = 0;
            calc_trig(data);
        };

    // Add Exchange gradient
    if( options.J != 0)
    {
        new_f = [init_f, options, d](
            std::valarray<double>& grad_out,
            const std::valarray<double>& data)
        {
            init_f(grad_out, data);
            exchange_grad(grad_out, options.J, d);
        };
        init_f = new_f;
    }

    // Add Zeeman gradient
    if( options.H != 0)
    {
        new_f = [init_f, options](
            std::valarray<double>& grad_out,
            const std::valarray<double>& data)
        {
            init_f(grad_out, data);
            zeeman_grad(grad_out, options.H);
        };
        init_f = new_f;
    }

    // Set work arrs and slices
    set_slices(size, d);

    return init_f;
}

void hmc::calc_trig(const std::valarray<double> &data)
{
    cos_the = cos(data[tslice]);
    sin_the = sin(data[tslice]);
    cos_phi = cos(data[pslice]);
    sin_phi = sin(data[pslice]);
}

void hmc::set_slices(int size, int dim)
{
    // Set arrays
    large_slice.resize(dim);
    small_slice.resize(dim);
    opp_small_slice.resize(dim);
    opp_large_slice.resize(dim);
    temp.resize(size);

    // Set slices
    halfsize = size / 2;
    int sidesize = pow(halfsize, 1./float(dim));
    int sidesize2 = sidesize*sidesize;
    tslice = std::slice(0, halfsize, 1);
    pslice = std::slice(halfsize, halfsize, 1);
    small_temp.resize(halfsize);
    cos_the.resize(halfsize);
    sin_the.resize(halfsize);
    cos_phi.resize(halfsize);
    sin_phi.resize(halfsize);
    switch(dim)
    {
        case 1:
        large_slice[0] = std::gslice(0, {halfsize-1}, {1});
        small_slice[0] = std::gslice(halfsize-1, {1}, {halfsize});
        opp_large_slice[0] = std::gslice(1, {halfsize-1}, {1});
        opp_small_slice[0] = std::gslice(0, {1}, {halfsize});
        break;

        case 2:
        large_slice[0] = std::gslice(0, {sidesize, sidesize-1}, {sidesize, 1});
        opp_large_slice[0] = std::gslice(1, {sidesize, sidesize-1}, {sidesize, 1});
        opp_small_slice[0] = std::gslice(0, {sidesize}, {sidesize});
        small_slice[0] = std::gslice(sidesize-1, {sidesize}, {sidesize});
        large_slice[1] = std::gslice(0, {sidesize2-sidesize}, {1});
        opp_large_slice[1] = std::gslice(sidesize, {sidesize2-sidesize}, {1});
        small_slice[1] = std::gslice(sidesize2-sidesize, {sidesize}, {1});
        opp_small_slice[1] = std::gslice(0, {sidesize}, {1});
        break;

        case 3:
        large_slice[0] = std::gslice(0, {sidesize2, sidesize-1}, {sidesize, 1});
        opp_large_slice[0] = std::gslice(1, {sidesize2, sidesize-1}, {sidesize, 1});
        opp_small_slice[0] = std::gslice(0, {sidesize2}, {sidesize});
        small_slice[0] = std::gslice(sidesize-1, {sidesize2}, {sidesize});
        large_slice[1] = std::gslice(0, {sidesize, sidesize2-sidesize}, {sidesize2, 1});
        opp_large_slice[1] = std::gslice(sidesize, {sidesize, sidesize2-sidesize}, {sidesize2, 1});
        opp_small_slice[1] = std::gslice(0, {sidesize, sidesize}, {sidesize2, 1});
        small_slice[1] = std::gslice(sidesize2-sidesize, {sidesize, sidesize}, {sidesize2, 1});
        large_slice[2] = std::gslice(0, {halfsize-sidesize2}, {1});
        opp_large_slice[2] = std::gslice(sidesize2, {halfsize-sidesize2}, {1});
        small_slice[2] = std::gslice(halfsize-sidesize2, {sidesize2}, {1});
        opp_small_slice[2] = std::gslice(0, {sidesize2}, {1});
        break;
    }
}
