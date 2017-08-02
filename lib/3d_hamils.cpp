#include "../include/3d_hamils.hpp"
#include <cmath>

std::vector<std::valarray<double> > hmc::trig_luf(const std::valarray<double> data)
{
    int halfsize = data.size() / 2;
    int size = pow(halfsize, 1./3.);
    int size2 = size*size;
    std::slice tslice(0, halfsize, 1);
    std::slice pslice(halfsize, halfsize, 1);

    std::gslice full_left_slice(0, {size2, size-1}, {size, 1});
    std::gslice full_right_slice(1, {size2, size-1}, {size, 1});
    std::gslice left_slice(0, {size2}, {size});
    std::gslice right_slice(size-1, {size2}, {size});
    std::gslice full_up_slice(0, {size, size2-size}, {size2, 1});
    std::gslice full_down_slice(size, {size, size2-size}, {size2, 1});
    std::gslice up_slice(0, {size, size}, {size2, 1});
    std::gslice down_slice(size2-size, {size, size}, {size2, 1});

    std::vector<std::valarray<double> > E_input(16, std::valarray<double>(halfsize));
    E_input[0] = cos(data[tslice]);
    E_input[4] = sin(data[tslice]);
    E_input[8] = cos(data[pslice]);
    E_input[12] = sin(data[pslice]);
    for(int i = 0; i < 4; i++)
    {
        E_input[4*i+1][full_right_slice] = E_input[4*i][full_left_slice];
        E_input[4*i+1][left_slice] = E_input[4*i][right_slice];
        E_input[4*i+2][full_down_slice] = E_input[4*i][full_up_slice];
        E_input[4*i+2][up_slice] = E_input[4*i][down_slice];
        E_input[4*i+3] = E_input[4*i].cshift(-size2);
    }

    return E_input;
}

std::vector<std::valarray<double> > hmc::trig_lrudfb(const std::valarray<double> data)
{
    int halfsize = data.size() / 2;
    int size = pow(halfsize, 1./3.);
    int size2 = size*size;
    std::slice tslice(0, halfsize, 1);
    std::slice pslice(halfsize, halfsize, 1);

    std::gslice full_left_slice(0, {size2, size-1}, {size, 1});
    std::gslice full_right_slice(1, {size2, size-1}, {size, 1});
    std::gslice left_slice(0, {size2}, {size});
    std::gslice right_slice(size-1, {size2}, {size});
    std::gslice full_up_slice(0, {size, size2-size}, {size2, 1});
    std::gslice full_down_slice(size, {size, size2-size}, {size2, 1});
    std::gslice up_slice(0, {size, size}, {size2, 1});
    std::gslice down_slice(size2-size, {size, size}, {size2, 1});

    std::vector<std::valarray<double> > E_input(28, std::valarray<double>(halfsize));
    E_input[0] = cos(data[tslice]);
    E_input[7] = sin(data[tslice]);
    E_input[14] = cos(data[pslice]);
    E_input[21] = sin(data[pslice]);
    for(int i = 0; i < 4; i++)
    {
        E_input[7*i+1][full_right_slice] = E_input[7*i][full_left_slice];
        E_input[7*i+1][left_slice] = E_input[7*i][right_slice];
        E_input[7*i+2][full_left_slice] = E_input[7*i][full_right_slice];
        E_input[7*i+2][right_slice] = E_input[7*i][left_slice];
        E_input[7*i+3][full_down_slice] = E_input[7*i][full_up_slice];
        E_input[7*i+3][up_slice] = E_input[7*i][down_slice];
        E_input[7*i+4][full_up_slice] = E_input[7*i][full_down_slice];
        E_input[7*i+4][down_slice] = E_input[7*i][up_slice];
        E_input[7*i+5] = E_input[7*i].cshift(-size2);
        E_input[7*i+6] = E_input[7*i].cshift(size2);
    }

    return E_input;
}
