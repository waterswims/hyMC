#include "../include/2d_hamils.hpp"
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

std::vector<std::valarray<double> > hmc::trig_left_up(const std::valarray<double> data)
{
    int halfsize = data.size() / 2;
    int size = sqrt(halfsize);
    std::slice tslice(0, halfsize, 1);
    std::slice pslice(halfsize, halfsize, 1);

    std::gslice full_left_slice(0, {size, size-1}, {size, 1});
    std::gslice full_right_slice(1, {size, size-1}, {size, 1});
    std::gslice left_slice(0, {size, 1}, {size, 1});
    std::gslice right_slice(size-1, {size, 1}, {size, 1});

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

std::function<double(const std::valarray<double>&)>
    hmc::gen_total_energy_2d(std::vector<bool> E_flags)
{
    // Init blank energy function
    std::function<double(const std::vector<std::valarray<double> >&)> init_f =
        [](const std::vector<std::valarray<double> >& data){return 0;};
    std::function<double(const std::vector<std::valarray<double> >&)> new_f;

    // Add exchange energy
    if (E_flags[0])
    {
        new_f =
            [init_f](const std::vector<std::valarray<double> >& data)
            {return init_f(data) + exchange_energy_2d(data);};
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
