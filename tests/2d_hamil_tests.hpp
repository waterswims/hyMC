#ifndef HAMIL_2D_TEST
#define HAMIL_2D_TEST

#include "../include/2d_hamils.hpp"
#include <gtest/gtest.h>

TEST(Hamiltonian_2d, trig_left_up)
{
    int tsize = 9;
    double pi = 3.14159265359;
    std::valarray<double> test_spins(tsize*2);
    test_spins = 1.5;
    test_spins[4] = pi / 2;
    test_spins[16] = pi / 2;

    std::vector<std::valarray<double> > test = hmc::trig_left_up(test_spins);
    EXPECT_NEAR(test[0][4], 0, 1e-6);
    EXPECT_NEAR(test[1][5], 0, 1e-6);
    EXPECT_NEAR(test[2][7], 0, 1e-6);
    EXPECT_NEAR(test[3][4], 1, 1e-6);
    EXPECT_NEAR(test[4][5], 1, 1e-6);
    EXPECT_NEAR(test[5][7], 1, 1e-6);
    EXPECT_NEAR(test[6][7], 0, 1e-6);
    EXPECT_NEAR(test[7][8], 0, 1e-6);
    EXPECT_NEAR(test[8][1], 0, 1e-6);
    EXPECT_NEAR(test[9][7], 1, 1e-6);
    EXPECT_NEAR(test[10][8], 1, 1e-6);
    EXPECT_NEAR(test[11][1], 1, 1e-6);
}

TEST(Hamiltonian_2d, exchange_alligned)
{
    int size = 4;
    int tsize = size*size;
    std::valarray<double> test_spins(tsize*2);
    for(int i = 0; i < tsize; i++)
    {
        test_spins[i] = 4.5;
        test_spins[tsize+i] = 1.2;
    }

    std::vector<std::valarray<double> > E_input = hmc::trig_left_up(test_spins);

    EXPECT_NEAR(hmc::exchange_energy_2d(E_input), -2*tsize, 2*tsize*1e-15);
}

TEST(Hamiltonian_2d, exchange_antialligned)
{
    int size = 4;
    int tsize = size*size;
    std::valarray<double> test_spins(tsize*2);
    for(int i = 0; i < tsize; i++)
    {
        test_spins[i] = ((i/4)%2+(i%2))*pi + 1.2;
        test_spins[tsize+i] = (((i/4)%2)*2-1)*((i%2)*2-1)*0.4 + pi/2;
    }

    std::vector<std::valarray<double> > E_input = hmc::trig_left_up(test_spins);

    EXPECT_NEAR(hmc::exchange_energy_2d(E_input), 2*tsize, 2*tsize*1e-15);
}

TEST(Hamiltonian_2d, exchange_specific)
{
    int size = 4;
    std::valarray<double> test_spins(size*2);

    test_spins[0] = 0.3;
    test_spins[4] = 2.3;

    test_spins[1] = 1.6;
    test_spins[5] = 0.1;

    test_spins[2] = 5.2;
    test_spins[6] = 1.1;

    test_spins[3] = 3.4;
    test_spins[7] = 2.1;

    std::vector<std::valarray<double> > E_input = hmc::trig_left_up(test_spins);
    EXPECT_NEAR(hmc::exchange_energy_2d(E_input), 3.4939748074272869,
                                        3.4939748074272869*1e-15);
}

TEST(Hamiltonian_2d, exchange_specific)
{
    int size = 4;
    std::valarray<double> test_spins(size*2);

    test_spins[0] = 0.3;
    test_spins[4] = 2.3;

    test_spins[1] = 1.6;
    test_spins[5] = 0.1;

    test_spins[2] = 5.2;
    test_spins[6] = 1.1;

    test_spins[3] = 3.4;
    test_spins[7] = 2.1;

    std::vector<bool> E_flags(1);
    std::function<double(const std::valarray<double>&)> E_func;

    E_flags[0] = true;
    E_func = hmc::gen_total_energy_2d(E_flags);
    EXPECT_NEAR(E_func(test_spins), 3.4939748074272869,
                                        3.4939748074272869*1e-15);

    E_flags[0] = false;
    E_func = hmc::gen_total_energy_2d(E_flags);
    EXPECT_EQ(E_func(test_spins), 0);
}

#endif
