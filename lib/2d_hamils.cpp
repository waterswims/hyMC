#include "../include/2d_hamils.hpp"
#include "../include/hmc.hpp"
#include <cmath>
#include <iostream>

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
