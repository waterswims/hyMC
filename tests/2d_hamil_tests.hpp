#ifndef HAMIL_2D_TEST
#define HAMIL_2D_TEST

#include "../include/all_hamils.hpp"
#include "../include/hmc.hpp"
#include <gtest/gtest.h>
#include <iostream>

TEST(Hamiltonian_2d, exchange_alligned)
{
    int size = 4;
    int tsize = size*size;
    std::valarray<double> test_spins(tsize*2);
    for(int i = 0; i < tsize; i++)
    {
        test_spins[i] = 0.25;
        test_spins[tsize+i] = 0.5;
    }
    hmc::set_slices(tsize*2, 2);
    hmc::calc_trig(test_spins);
    EXPECT_NEAR(hmc::exchange_energy(1, 2), -2*tsize, 2*tsize*1e-15);
}

TEST(Hamiltonian_2d, exchange_antialligned)
{
    int size = 4;
    int tsize = size*size;
    std::valarray<double> test_spins(tsize*2);

    for(int i = 0; i < tsize; i++)
    {
        test_spins[tsize+i] = ((i/4)%2+(i%2))*pi + 1.2;
        test_spins[i] = (((i/4)%2)*2-1)*((i%2)*2-1)*0.4 + pi/2;
    }
    hmc::set_slices(tsize*2, 2);
    hmc::calc_trig(test_spins);
    EXPECT_NEAR(hmc::exchange_energy(1, 2), 2*tsize, 2*tsize*1e-15);
}

TEST(Hamiltonian_2d, exchange_specific)
{
    int size = 4;
    std::valarray<double> test_spins(size*2);

    test_spins[0] = 2.3;
    test_spins[1] = 0.1;
    test_spins[2] = 1.1;
    test_spins[3] = 2.1;
    test_spins[4] = 0.3;
    test_spins[5] = 1.6;
    test_spins[6] = 5.2;
    test_spins[7] = 3.4;

    hmc::set_slices(size*2, 2);
    hmc::calc_trig(test_spins);
    EXPECT_NEAR(hmc::exchange_energy(1, 2), 3.4939748074272869,
                                        3.4939748074272869*1e-8);
}

TEST(Hamiltonian_2d, total_energy)
{
    int size = 4;
    std::valarray<double> test_spins(size*2);

    test_spins[0] = 2.3;
    test_spins[1] = 0.1;
    test_spins[2] = 1.1;
    test_spins[3] = 2.1;
    test_spins[4] = 0.3;
    test_spins[5] = 1.6;
    test_spins[6] = 5.2;
    test_spins[7] = 3.4;

    struct hmc::HamiltonianOptions options;
    double beta = 1.0;
    std::function<double(const std::valarray<double>&)> E_func;

    options.J = 1.0;
    options.H = 0.0;
    E_func = hmc::gen_total_energy( options, beta, 2, size*2 );
    EXPECT_NEAR(E_func(test_spins), 3.4939748074272869,
                                        3.4939748074272869*1e-9);

    options.J = 0.0;
    E_func = hmc::gen_total_energy( options, beta, 2, size*2 );
    EXPECT_EQ(E_func(test_spins), 0);

    options.J = 1.0;
    options.H = 2.1;
    E_func = hmc::gen_total_energy( options, beta, 2, size*2 );
    EXPECT_NEAR(E_func(test_spins), 2.91127067, 2.91127067*1e-9);
}

TEST(Hamiltonian_2d, exchange_grad)
{
    int size = 4;
    std::valarray<double> test_spins(size*2);

    test_spins[0] = 2.3;
    test_spins[1] = 0.1;
    test_spins[2] = 1.1;
    test_spins[3] = 2.1;
    test_spins[4] = 0.3;
    test_spins[5] = 1.6;
    test_spins[6] = 5.2;
    test_spins[7] = 3.4;

    std::valarray<double> grad(size*2);
    grad = 0;
    hmc::set_slices(size*2, 2);
    hmc::calc_trig(test_spins);
    hmc::exchange_grad(grad, 1, 2);
    EXPECT_FLOAT_EQ(grad[0], -2.4175421944456765);
    EXPECT_FLOAT_EQ(grad[1], 0.24050534024646605);
    EXPECT_FLOAT_EQ(grad[2], 2.0356793154147046);
    EXPECT_FLOAT_EQ(grad[3], -2.2735417704169287);
    EXPECT_FLOAT_EQ(grad[4], -1.1623660509440559);
    EXPECT_FLOAT_EQ(grad[5], 0.024380126112986563);
    EXPECT_FLOAT_EQ(grad[6], -0.19252639008304462);
    EXPECT_FLOAT_EQ(grad[7], 1.3305123149141138);
}

TEST(Hamiltonian_2d, total_grad)
{
    int size = 4;
    std::valarray<double> test_spins(size*2);

    test_spins[0] = 2.3;
    test_spins[1] = 0.1;
    test_spins[2] = 1.1;
    test_spins[3] = 2.1;
    test_spins[4] = 0.3;
    test_spins[5] = 1.6;
    test_spins[6] = 5.2;
    test_spins[7] = 3.4;

    std::valarray<double> grad(size*2);
    struct hmc::HamiltonianOptions options;
    std::function<void(std::valarray<double>&, const std::valarray<double>&)> g_func;

    options.J = 1.0;
    options.H = 0;
    g_func = hmc::gen_total_grad( options, 2, size*2 );
    g_func(grad, test_spins);
    EXPECT_FLOAT_EQ(grad[0], -2.4175421944456765);
    EXPECT_FLOAT_EQ(grad[1], 0.24050534024646605);
    EXPECT_FLOAT_EQ(grad[2], 2.0356793154147046);
    EXPECT_FLOAT_EQ(grad[3], -2.2735417704169287);
    EXPECT_FLOAT_EQ(grad[4], -1.1623660509440559);
    EXPECT_FLOAT_EQ(grad[5], 0.024380126112986563);
    EXPECT_FLOAT_EQ(grad[6], -0.19252639008304462);
    EXPECT_FLOAT_EQ(grad[7], 1.3305123149141138);

    options.J = 0.0;
    g_func = hmc::gen_total_grad( options, 2, size*2 );
    g_func(grad, test_spins);
    EXPECT_EQ(grad[0], 0);
    EXPECT_EQ(grad[1], 0);
    EXPECT_EQ(grad[2], 0);
    EXPECT_EQ(grad[3], 0);
    EXPECT_EQ(grad[4], 0);
    EXPECT_EQ(grad[5], 0);
    EXPECT_EQ(grad[6], 0);
    EXPECT_EQ(grad[7], 0);

    options.J = 1.0;
    options.H = 2.1;
    g_func = hmc::gen_total_grad( options, 2, size*2 );
    g_func(grad, test_spins);
    EXPECT_NEAR(grad[0], -0.8515612488, 0.8515612488*1e-9);
    EXPECT_NEAR(grad[1], 0.4501555152, 0.4501555152*1e-9);
    EXPECT_NEAR(grad[2], 3.907214772, 3.907214772*1e-9);
    EXPECT_NEAR(grad[3], -0.4608021006, 0.4608021006*1e-9);
    EXPECT_FLOAT_EQ(grad[4], -1.1623660509440559);
    EXPECT_FLOAT_EQ(grad[5], 0.024380126112986563);
    EXPECT_FLOAT_EQ(grad[6], -0.19252639008304462);
    EXPECT_FLOAT_EQ(grad[7], 1.3305123149141138);
}

#endif
