#ifndef HAMIL_1D
#define HAMIL_1D

#include "spins_class.hpp"

double total_energy_1d(const spin_lattice_1d &spins,
    const spin_lattice_1d &vels,
    const bool* energy_flags);

void total_energy_grad_1d(const spin_lattice_1d &spins,
    spin_lattice_1d &grad_out,
    const bool* energy_flags);

double exchange_1d(const spin_lattice_1d &spins);

void exchange_grad_1d(const spin_lattice_1d &spins, spin_lattice_1d &grad_out);

double kinetic_1d(const spin_lattice_1d &vels);

#endif
