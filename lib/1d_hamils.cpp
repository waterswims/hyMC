#include "../include/1d_hamils.hpp"
#include "../include/array_alloc.hpp"
#include "../include/vector_trig.hpp"
#include "../include/spins_class.hpp"
#include <cmath>

double total_energy_1d(const spin_lattice_1d &spins,
    const spin_lattice_1d &vels,
    const bool* energy_flags)
{
    // Start with kinetic energy
    double E = kinetic_1d(vels);

    // Add Exchange
    if (energy_flags[0])
    {
        E += exchange_1d(spins);
    }

    return E;
}

double exchange_1d(const spin_lattice_1d &spins)
{
    double* cos_theta = alloc_1darr<double>(spins.N);
    double* cos_phi = alloc_1darr<double>(spins.N);
    double* sin_theta = alloc_1darr<double>(spins.N);
    double* sin_phi = alloc_1darr<double>(spins.N);

    cos_1dvec(spins.thetas, cos_theta, spins.N);
    cos_1dvec(spins.phis, cos_phi, spins.N);
    sin_1dvec(spins.thetas, sin_theta, spins.N);
    sin_1dvec(spins.phis, sin_phi, spins.N);

    double E = 0;

    #pragma simd
    for(int i = 0; i < spins.N; i++)
    {
        int right = (i + 1)%spins.N;
        double t1 = cos_phi[i] * cos_phi[right];
        double t2 = cos_theta[i] * cos_theta[right] * sin_phi[i] * sin_phi[right];
        double t3 = sin_theta[i] * sin_phi[i] * sin_theta[right] * sin_phi[right];
        E -= t1 + t2 + t3;
    }

    dealloc_1darr<double>(cos_theta);
    dealloc_1darr<double>(cos_phi);
    dealloc_1darr<double>(sin_theta);
    dealloc_1darr<double>(sin_phi);

    return E;
}

void exchange_grad_1d(const spin_lattice_1d &spins, spin_lattice_1d &grad_out)
{
    double* cos_theta = alloc_1darr<double>(spins.N);
    double* cos_phi = alloc_1darr<double>(spins.N);
    double* sin_theta = alloc_1darr<double>(spins.N);
    double* sin_phi = alloc_1darr<double>(spins.N);

    cos_1dvec(spins.thetas, cos_theta, spins.N);
    cos_1dvec(spins.phis, cos_phi, spins.N);
    sin_1dvec(spins.thetas, sin_theta, spins.N);
    sin_1dvec(spins.phis, sin_phi, spins.N);

    for(int i = 0; i < spins.N; i++)
    {
        int left = (spins.N + ((i - 1) % spins.N)) % spins.N;
        int right = (i + 1) % spins.N;

        double dt_t1 = -cos_theta[left] * sin_theta[i] * sin_phi[i] * sin_phi[left];
        double dt_t2 = cos_theta[i] * sin_phi[i] * sin_theta[left] * sin_phi[left];
        double dt_t3 = -cos_theta[right] * sin_theta[i] * sin_phi[i] * sin_phi[right];
        double dt_t4 = cos_theta[i] * sin_phi[i] * sin_theta[right] * sin_phi[right];
        grad_out.thetas[i] = dt_t1 + dt_t2 + dt_t3 + dt_t4;

        double dp_t1 = -cos_phi[left] * sin_phi[i];
        double dp_t2 = cos_theta[i] * cos_phi[i] * cos_theta[left] * sin_phi[left];
        double dp_t3 = cos_phi[i] * sin_theta[i] * sin_theta[left] * sin_phi[left];
        double dp_t4 = -cos_phi[right] * sin_phi[i];
        double dp_t5 = cos_theta[i] * cos_phi[i] * cos_theta[right] * sin_phi[right];
        double dp_t6 = cos_phi[i] * sin_theta[i] * sin_theta[right] * sin_phi[right];
        grad_out.phis[i] = dp_t1 + dp_t2 + dp_t3 + dp_t4 + dp_t5 + dp_t6;
    }

    dealloc_1darr<double>(cos_theta);
    dealloc_1darr<double>(cos_phi);
    dealloc_1darr<double>(sin_theta);
    dealloc_1darr<double>(sin_phi);
}

double kinetic_1d(const spin_lattice_1d &vels)
{
    double s = 0;
    #pragma simd
    for(int i = 0; i < vels.N; i++)
    {
        s += pow(vels.thetas[i], 2);
        s += pow(vels.phis[i], 2);
    }
    double E = s / 2.;
    return E;
}
