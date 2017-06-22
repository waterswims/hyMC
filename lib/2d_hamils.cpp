#include "../include/2d_hamils.hpp"
#include "../include/hmc.hpp"
#include <cmath>

double hmc::exchange_energy_2d(const std::vector<std::valarray<double> >& trig_angles)
{
    double t1 = (trig_angles[6] * (trig_angles[7] + trig_angles[8])).sum();
    double t2 = (trig_angles[0] * trig_angles[9] *
                 (trig_angles[1] * trig_angles[10] +
                  trig_angles[2] * trig_angles[11])).sum();
    double t3 = (trig_angles[3] * trig_angles[9] *
                 (trig_angles[4] * trig_angles[10] +
                  trig_angles[5] * trig_angles[11])).sum();

    return -(t1 + t2 + t3);
}

void hmc::exchange_grad_2d(std::valarray<double>& grad_out,
    const std::vector<std::valarray<double> >& trig_angles)
{
    int halfsize = trig_angles[0].size();
    std::slice tslice(0, halfsize, 1);
    std::slice pslice(halfsize, halfsize, 1);
    std::valarray<double> t1, t2, t3;

    t1 = -trig_angles[5] * trig_angles[15] * (trig_angles[1] * trig_angles[16] +
        trig_angles[2] * trig_angles[17] + trig_angles[3] * trig_angles[18] +
        trig_angles[4] * trig_angles[19]);
    t2 = trig_angles[0] * trig_angles[15] * (trig_angles[6] * trig_angles[16] +
        trig_angles[7] * trig_angles[17] + trig_angles[8] * trig_angles[18] +
        trig_angles[9] * trig_angles[19]);

    grad_out[tslice] += t1 + t2;

    t1 = -trig_angles[15] * (trig_angles[11] + trig_angles[12] + trig_angles[13] +
                             trig_angles[14]);
    t2 = trig_angles[0] * trig_angles[10] * (trig_angles[1] * trig_angles[16] +
        trig_angles[2] * trig_angles[17] + trig_angles[3] * trig_angles[18] +
        trig_angles[4] * trig_angles[19]);
    t3 = trig_angles[10] * trig_angles[5] * (trig_angles[6] * trig_angles[16] +
        trig_angles[7] * trig_angles[17] + trig_angles[8] * trig_angles[18] +
        trig_angles[9] * trig_angles[19]);

    grad_out[pslice] += t1 + t2 + t3;
}

std::function<void(std::valarray<double>&, const std::valarray<double>&)>
hmc::gen_total_grad_2d( const HamiltonianOptions options )
{
    // Init blank grad function
    std::function<void(std::valarray<double>&, const std::vector<std::valarray<double> >&)> init_f, new_f;
    init_f = [](std::valarray<double>& grad_out, const std::vector<std::valarray<double> >& data)
        {grad_out = 0;};

    // Add Exchange gradient
    if( options.J != 0)
    {
        new_f = [init_f](std::valarray<double>& grad_out, const std::vector<std::valarray<double> >& data)
            {
                init_f(grad_out, data);
                exchange_grad_2d(grad_out, data);
            };
        init_f = new_f;
    }

    // Find trigs and finalise
    std::function<void(std::valarray<double>&, const std::valarray<double>&)> final_f =
        [init_f](std::valarray<double>& grad_out, const std::valarray<double>& data)
        {
            std::vector<std::valarray<double> > trigs = trig_lrud(data);
            init_f(grad_out, trigs);
        };

    return final_f;
}

std::vector<std::valarray<double> > hmc::trig_left_up(const std::valarray<double> data)
{
    int halfsize = data.size() / 2;
    int size = sqrt(halfsize);
    std::slice tslice(0, halfsize, 1);
    std::slice pslice(halfsize, halfsize, 1);

    std::gslice full_left_slice(0, {size, size-1}, {size, 1});
    std::gslice full_right_slice(1, {size, size-1}, {size, 1});
    std::gslice left_slice(0, {size}, {size});
    std::gslice right_slice(size-1, {size}, {size});

    std::vector<std::valarray<double> > E_input(12, std::valarray<double>(halfsize));
    E_input[0] = cos(data[tslice]);
    E_input[3] = sin(data[tslice]);
    E_input[6] = cos(data[pslice]);
    E_input[9] = sin(data[pslice]);
    for(int i = 0; i < 4; i++)
    {
        E_input[3*i+1][full_right_slice] = E_input[3*i][full_left_slice];
        E_input[3*i+1][left_slice] = E_input[3*i][right_slice];
        E_input[3*i+2] = E_input[3*i].cshift(-size);
    }

    return E_input;
}

std::vector<std::valarray<double> > hmc::trig_lrud(const std::valarray<double> data)
{
    int halfsize = data.size() / 2;
    int size = sqrt(halfsize);
    std::slice tslice(0, halfsize, 1);
    std::slice pslice(halfsize, halfsize, 1);

    std::gslice full_left_slice(0, {size, size-1}, {size, 1});
    std::gslice full_right_slice(1, {size, size-1}, {size, 1});
    std::gslice left_slice(0, {size, 1}, {size, 1});
    std::gslice right_slice(size-1, {size, 1}, {size, 1});

    std::vector<std::valarray<double> > E_input(20, std::valarray<double>(halfsize));
    E_input[0] = cos(data[tslice]);
    E_input[5] = sin(data[tslice]);
    E_input[10] = cos(data[pslice]);
    E_input[15] = sin(data[pslice]);
    for(int i = 0; i < 4; i++)
    {
        E_input[5*i+1][full_right_slice] = E_input[5*i][full_left_slice];
        E_input[5*i+1][left_slice] = E_input[5*i][right_slice];
        E_input[5*i+2][full_left_slice] = E_input[5*i][full_right_slice];
        E_input[5*i+2][right_slice] = E_input[5*i][left_slice];
        E_input[5*i+3] = E_input[5*i].cshift(-size);
        E_input[5*i+4] = E_input[5*i].cshift(size);
    }

    return E_input;
}

std::function<double(const std::valarray<double>&)>
hmc::gen_total_energy_2d( const HamiltonianOptions options, const double beta )
{
    // Init blank energy function
    std::function<double(const std::vector<std::valarray<double> >&)> init_f =
        [](const std::vector<std::valarray<double> >& data){return 0;};
    std::function<double(const std::vector<std::valarray<double> >&)> new_f;

    // Add exchange energy
    if ( options.J != 0)
    {
        new_f =
            [init_f, beta](const std::vector<std::valarray<double> >& data)
            {return init_f(data) + exchange_energy_2d(data) * beta;};
        init_f = new_f;
    }

    // Find trigs and feed into energy funcs
    std::function<double(const std::valarray<double>&)> final_f =
        [init_f](const std::valarray<double>& data)
        {
            std::vector<std::valarray<double> > trigs = trig_left_up(data);
            return init_f(trigs);
        };

    return final_f;
}
