#include "../include/spins_class.hpp"
#include "../include/array_alloc.hpp"

spin_lattice_1d::spin_lattice_1d(const int length)
{
    thetas = alloc_1darr<double>(length);
    phis = alloc_1darr<double>(length);
    N = length;
    mem_manage = true;
}

spin_lattice_1d::spin_lattice_1d(const int length, double *t, double *p)
{
    thetas = t;
    phis = p;
    N = length;
    mem_manage = false;
}

spin_lattice_1d::spin_lattice_1d(const spin_lattice_1d &cop)
{
    N = cop.N;
    thetas = deep_copy_1darr<double>(N, cop.thetas);
    phis = deep_copy_1darr<double>(N, cop.phis);
}

spin_lattice_1d::~spin_lattice_1d()
{
    if( mem_manage )
    {
        dealloc_1darr<double>(thetas);
        dealloc_1darr<double>(phis);
    }
}
