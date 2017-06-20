#include "spin_lattice_1d.hpp"
#include "array_alloc.hpp"

spin_lattice_1d::spin_lattice_1d(int length)
{
    thetas = alloc_1darr(length);
    phis = alloc_1darr(length);
}
