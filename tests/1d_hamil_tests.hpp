#ifndef HAMIL_1D_TEST
#define HAMIL_1D_TEST

#include "../include/1d_hamils.hpp"
#include "../include/spins_class.hpp"
#include "../include/array_alloc.hpp"
#include <gtest/gtest.h>

TEST(Hamiltonian_1d, exchange_alligned)
{
    int size = 20;
    spin_lattice_1d test_spins(size);
    for(int i = 0; i < size; i++)
    {
        test_spins.thetas[i] = 4.5;
        test_spins.phis[i] = 1.2;
    }
    EXPECT_NEAR(exchange_1d(test_spins), -size, size*1e-15);
}

TEST(Hamiltonian_1d, exchange_antialligned)
{
    int size = 20;
    spin_lattice_1d test_spins(size);
    for(int i = 0; i < size; i++)
    {
        test_spins.thetas[i] = (i%2)*pi + 1.2;
        test_spins.phis[i] = ((i%2)*2-1)*0.4 + pi/2;
    }
    EXPECT_NEAR(exchange_1d(test_spins), size, size*1e-15);
}

TEST(Hamiltonian_1d, exchange_specific)
{
    int size = 3;
    spin_lattice_1d test_spins(size);
    test_spins.thetas[0] = 0.3;
    test_spins.phis[0] = 2.3;
    test_spins.thetas[1] = 1.6;
    test_spins.phis[1] = 0.1;
    test_spins.thetas[2] = 5.2;
    test_spins.phis[2] = 1.1;
    EXPECT_NEAR(exchange_1d(test_spins), 0.44975793200288505, 0.44975793200288505*1e-15);
}

TEST(Hamiltonian_1d, exchange_grad)
{
    int size = 3;
    spin_lattice_1d test_spins(size);
    test_spins.thetas[0] = 0.3;
    test_spins.phis[0] = 2.3;
    test_spins.thetas[1] = 1.6;
    test_spins.phis[1] = 0.1;
    test_spins.thetas[2] = 5.2;
    test_spins.phis[2] = 1.1;

    spin_lattice_1d grad(size);
    exchange_grad_1d(test_spins, grad);
    EXPECT_NEAR(grad.thetas[0], -0.58118303, 0.58118303*1e-7);
    EXPECT_NEAR(grad.thetas[1], -0.11110539, 0.11110539*1e-7);
    EXPECT_NEAR(grad.thetas[2], 0.69228842, 0.69228842*1e-7);
    EXPECT_NEAR(grad.phis[0], -1.2087711, 1.2087711*1e-7);
    EXPECT_NEAR(grad.phis[1], -0.57549375, 0.57549375*1e-7);
    EXPECT_NEAR(grad.phis[2], -0.27048617, 0.27048617*1e-7);
}

TEST(Hamiltonian_1d, kinetic)
{
    int size = 10;
    spin_lattice_1d test_vels(size);
    for(int i = 0; i < size; i++)
    {
        test_vels.thetas[i] = i + 2.1;
        test_vels.phis[i] = -i + 4.6;
    }
    EXPECT_NEAR(kinetic_1d(test_vels), 300.34999999999997, 300.34999999999997*1e-15);
}

TEST(Hamiltonian_1d, total)
{
    int size = 3;
    spin_lattice_1d test_spins(size);
    test_spins.thetas[0] = 0.3;
    test_spins.phis[0] = 2.3;
    test_spins.thetas[1] = 1.6;
    test_spins.phis[1] = 0.1;
    test_spins.thetas[2] = 5.2;
    test_spins.phis[2] = 1.1;

    spin_lattice_1d test_vels(size);
    test_vels.thetas[0] = 1.1;
    test_vels.phis[0] = 0.2;
    test_vels.thetas[1] = 2.1;
    test_vels.phis[1] = 0.6;
    test_vels.thetas[2] = -0.2;
    test_vels.phis[2] = -1;

    bool* E_flags = alloc_1darr<bool>(1);
    E_flags[0] = true;
    EXPECT_NEAR(total_energy_1d(test_spins, test_vels, E_flags),
        3.9797579320028853, 3.9797579320028853*1e-15);

    E_flags[0] = false;
    EXPECT_NEAR(total_energy_1d(test_spins, test_vels, E_flags),
        3.5300000000000002, 3.5300000000000002*1e-15);

    dealloc_1darr<bool>(E_flags);
}

#endif
