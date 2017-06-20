#include "../include/spins_class.hpp"
#include "../include/array_alloc.hpp"

spin_lattice_1d::spin_lattice_1d(int length)
{
    thetas = alloc_1darr(length);
    phis = alloc_1darr(length);
    N = length;
}

spin_lattice_1d::~spin_lattice_1d()
{
    dealloc_1darr(thetas);
    dealloc_1darr(phis);
}
