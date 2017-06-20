#include "../include/1d_hamils.hpp"
#include "../include/array_alloc.hpp"
#include "../include/vector_trig.hpp"
#include "../include/spins_class.hpp"

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
