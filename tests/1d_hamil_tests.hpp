#ifndef HAMIL_1D_TEST
#define HAMIL_1D_TEST

#include "../include/1d_hamils.hpp"
#include "../include/spins_class.hpp"
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

#endif
