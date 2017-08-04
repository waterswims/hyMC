#include "../include/all_hamils.hpp"
#include "../include/1d_hamils.hpp"
#include "../include/2d_hamils.hpp"
#include "../include/3d_hamils.hpp"
#include <iostream>

double hmc::zeeman_energy(
    const std::vector<std::valarray<double> >& trig_angles,
    const double H)
{
    int ind = trig_angles.size() / 2;
    return -(H * trig_angles[ind]).sum();
}

void hmc::zeeman_grad(
    std::valarray<double>& grad_out,
    const std::vector<std::valarray<double> >& trig_angles,
    const double H)
{
    int halfsize = trig_angles[0].size();
    std::slice pslice(halfsize, halfsize, 1);
    int d = trig_angles.size() / 4;

    grad_out[pslice] += H * trig_angles[d*3];
}

double hmc::exchange_energy(
    const std::vector<std::valarray<double> >& trig_angles,
    const double J,
    const int d)
{
    std::valarray<double> temp(trig_angles[0].size());
    temp = 0;
    for(int i = 0; i < d; i++)
    {
        temp += trig_angles[2*(d+1)+1+i];
    }
    double t1 = (trig_angles[2*(d+1)] * temp).sum();
    temp = 0;
    for(int i = 0; i < d; i++)
    {
        temp += trig_angles[1+i] * trig_angles[3*(d+1)+1+i];
    }
    double t2 = (trig_angles[0] * trig_angles[3*(d+1)] * temp).sum();
    temp = 0;
    for(int i = 0; i < d; i++)
    {
        temp += trig_angles[d+2+i] * trig_angles[3*(d+1)+1+i];
    }
    double t3 = (trig_angles[d+1] * trig_angles[3*(d+1)] * temp).sum();

    return -J * (t1 + t2 + t3);
}

void hmc::exchange_grad(
    std::valarray<double>& grad_out,
    const std::vector<std::valarray<double> >& trig_angles,
    const double J,
    const int d)
{
    int halfsize = trig_angles[0].size();
    std::slice tslice(0, halfsize, 1);
    std::slice pslice(halfsize, halfsize, 1);
    std::valarray<double> temp(trig_angles[0].size());

    temp = 0;
    for(int i = 0; i < 2*d; i++)
    {
        temp += trig_angles[1+i] * trig_angles[3*(2*d+1)+1+i];
    }
    grad_out[tslice] += -trig_angles[2*d+1] * trig_angles[3*(2*d+1)] * temp;
    grad_out[pslice] += trig_angles[0] * trig_angles[2*(2*d+1)] * temp;
    temp = 0;
    for(int i = 0; i < 2*d; i++)
    {
        temp += trig_angles[2*d+2+i] * trig_angles[3*(2*d+1)+1+i];
    }
    grad_out[tslice] += trig_angles[0] * trig_angles[3*(2*d+1)] * temp;
    grad_out[pslice] += trig_angles[2*d+1] * trig_angles[2*(2*d+1)] * temp;

    temp = 0;
    for(int i = 0; i < 2*d; i++)
    {
        temp += trig_angles[2*(2*d+1)+1+i];
    }
    grad_out[pslice] += -trig_angles[3*(2*d+1)] * temp;
    temp = 0;
}

std::function<double(const std::valarray<double>&)> hmc::gen_total_energy(
    const HamiltonianOptions options,
    const double beta,
    const int d)
{
    // Init blank energy function
    std::function<double(const std::vector<std::valarray<double> >&)> init_f =
        [](const std::vector<std::valarray<double> >& data){return 0;};
    std::function<double(const std::vector<std::valarray<double> >&)> new_f;

    // Add exchange energy
    if ( options.J != 0)
    {
        new_f =
            [init_f, options, d](const std::vector<std::valarray<double> >& data)
            {return init_f(data) + exchange_energy(data, options.J, d);};
        init_f = new_f;
    }

    // Add zeeman energy
    if ( options.H != 0 )
    {
        new_f =
            [init_f, options](const std::vector<std::valarray<double> >& data)
            {return init_f(data) + zeeman_energy(data, options.H);};
        init_f = new_f;
    }

    // Find trigs and feed into energy funcs for different dimensions
    std::function<double(const std::valarray<double>&)> final_f;
    switch(d)
    {
    case 1:
        final_f =
            [init_f, beta](const std::valarray<double>& data)
            {
                std::vector<std::valarray<double> > trigs = trig_left(data);
                return init_f(trigs) * beta;
            };
        break;
    case 2:
        final_f =
            [init_f, beta](const std::valarray<double>& data)
            {
                std::vector<std::valarray<double> > trigs = trig_left_up(data);
                return init_f(trigs) * beta;
            };
        break;
    case 3:
        final_f =
            [init_f, beta](const std::valarray<double>& data)
            {
                std::vector<std::valarray<double> > trigs = trig_luf(data);
                return init_f(trigs) * beta;
            };
        break;
    default:
        std::cerr << "Dimension of the problem should be between 1 and 3"
                  << std::endl;
    }

    return final_f;
}

std::function<void(std::valarray<double>&, const std::valarray<double>&)>
hmc::gen_total_grad( const HamiltonianOptions options, const int d)
{
    // Init blank grad function
    std::function<void(std::valarray<double>&, const std::vector<std::valarray<double> >&)> init_f, new_f;
    init_f = [](std::valarray<double>& grad_out, const std::vector<std::valarray<double> >& data)
        {grad_out = 0;};

    // Add Exchange gradient
    if( options.J != 0)
    {
        new_f = [init_f, options, d](std::valarray<double>& grad_out, const std::vector<std::valarray<double> >& data)
        {
            init_f(grad_out, data);
            exchange_grad(grad_out, data, options.J, d);
        };
        init_f = new_f;
    }

    // Add Zeeman gradient
    if( options.H != 0)
    {
        new_f = [init_f, options](std::valarray<double>& grad_out, const std::vector<std::valarray<double> >& data)
        {
            init_f(grad_out, data);
            zeeman_grad(grad_out, data, options.H);
        };
        init_f = new_f;
    }

    // Find trigs and finalise for different dims
    std::function<void(std::valarray<double>&, const std::valarray<double>&)> final_f;
    switch(d)
    {
    case 1:
        final_f = [init_f](std::valarray<double>& grad_out, const std::valarray<double>& data)
            {
                std::vector<std::valarray<double> > trigs = trig_lr(data);
                init_f(grad_out, trigs);
            };
        break;
    case 2:
        final_f = [init_f](std::valarray<double>& grad_out, const std::valarray<double>& data)
            {
                std::vector<std::valarray<double> > trigs = trig_lrud(data);
                init_f(grad_out, trigs);
            };
        break;
    case 3:
        final_f = [init_f](std::valarray<double>& grad_out, const std::valarray<double>& data)
            {
                std::vector<std::valarray<double> > trigs = trig_lrudfb(data);
                init_f(grad_out, trigs);
            };
        break;
    default:
        std::cerr << "Dimension of the problem should be between 1 and 3"
                  << std::endl;
    }

    return final_f;
}
