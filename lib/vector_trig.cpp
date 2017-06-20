#include "../include/vector_trig.hpp"
#include <cmath>

void cos_1dvec(double* inp_vec, double* out_vec, int size)
{
    #pragma simd
    for(int i = 0; i < size; i++)
    {
        out_vec[i] = cos(inp_vec[i]);
    }
}

void sin_1dvec(double* inp_vec, double* out_vec, int size)
{
    #pragma simd
    for(int i = 0; i < size; i++)
    {
        out_vec[i] = sin(inp_vec[i]);
    }
}
