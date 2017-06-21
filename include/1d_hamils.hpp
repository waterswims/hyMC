#ifndef HAMIL_1D
#define HAMIL_1D

#include "spins_class.hpp"

double total_energy_1d(const spin_lattice_1d &spins,
    const spin_lattice_1d &vels,
    const bool* energy_flags);

double exchange_1d(const spin_lattice_1d &spins);

double kinetic_1d(const spin_lattice_1d &vels);

#endif
