#ifndef HAMIL_1D_TEST
#define HAMIL_1D_TEST

#include "../include/all_hamils.hpp"
#include "../include/1d_hamils.hpp"
#include "../include/hmc.hpp"
#include <gtest/gtest.h>

TEST(Hamiltonian_1d, exchange_alligned)
{
    int size = 20;
    std::valarray<double> test_spins(size*2);
    for(int i = 0; i < size; i++)
    {
        test_spins[i] = 4.5;
        test_spins[size+i] = 1.2;
    }
    std::vector<std::valarray<double> > E_input = hmc::trig_left(test_spins);

    EXPECT_NEAR(hmc::exchange_energy(E_input, 1, 1), -size, size*1e-15);
}

TEST(Hamiltonian_1d, exchange_antialligned)
{
    int size = 20;
    std::valarray<double> test_spins(size*2);
    for(int i = 0; i < size; i++)
    {
        test_spins[i] = (i%2)*pi + 1.2;
        test_spins[size+i] = ((i%2)*2-1)*0.4 + pi/2;
    }
    std::vector<std::valarray<double> > E_input = hmc::trig_left(test_spins);

    EXPECT_NEAR(hmc::exchange_energy(E_input, 1, 1), size, size*1e-15);
}

TEST(Hamiltonian_1d, exchange_specific)
{
    int size = 3;
    std::valarray<double> test_spins(size*2);
    test_spins[0] = 0.3;
    test_spins[3] = 2.3;
    test_spins[1] = 1.6;
    test_spins[4] = 0.1;
    test_spins[2] = 5.2;
    test_spins[5] = 1.1;
    std::vector<std::valarray<double> > E_input = hmc::trig_left(test_spins);

    EXPECT_NEAR(hmc::exchange_energy(E_input, 1, 1), 0.44975793200288505, 0.44975793200288505*1e-15);
}

TEST(Hamiltonian_1d, total_energy)
{
    int size = 3;
    std::valarray<double> test_spins(size*2);
    test_spins[0] = 0.3;
    test_spins[3] = 2.3;
    test_spins[1] = 1.6;
    test_spins[4] = 0.1;
    test_spins[2] = 5.2;
    test_spins[5] = 1.1;

    struct hmc::HamiltonianOptions options;
    double beta = 1.0;
    std::function<double(const std::valarray<double>&)> E_func;

    options.J = 1.0;
    options.H = 0.0;
    E_func = hmc::gen_total_energy( options, beta, 1);
    EXPECT_NEAR(E_func(test_spins), 0.44975793200288505, 0.44975793200288505*1e-15);

    options.J = 0.0;
    E_func = hmc::gen_total_energy( options, beta, 1 );
    EXPECT_EQ(E_func(test_spins), 0);
}

TEST(Hamiltonian_1d, exchange_grad)
{
    int size = 3;
    std::valarray<double> test_spins(size*2);
    test_spins[0] = 0.3;
    test_spins[3] = 2.3;
    test_spins[1] = 1.6;
    test_spins[4] = 0.1;
    test_spins[2] = 5.2;
    test_spins[5] = 1.1;

    std::valarray<double> grad(size*2);
    grad = 0;

    std::vector<std::valarray<double> > g_input = hmc::trig_lr(test_spins);
    hmc::exchange_grad(grad, g_input, 1, 1);
    EXPECT_NEAR(grad[0], -0.58118303, 0.58118303*1e-7);
    EXPECT_NEAR(grad[1], -0.11110539, 0.11110539*1e-7);
    EXPECT_NEAR(grad[2], 0.69228842, 0.69228842*1e-7);
    EXPECT_NEAR(grad[3], -1.2087711, 1.2087711*1e-7);
    EXPECT_NEAR(grad[4], -0.57549375, 0.57549375*1e-7);
    EXPECT_NEAR(grad[5], -0.27048617, 0.27048617*1e-7);
}

TEST(Hamiltonian_1d, total_grad)
{
    int size = 3;
    std::valarray<double> test_spins(size*2);
    test_spins[0] = 0.3;
    test_spins[3] = 2.3;
    test_spins[1] = 1.6;
    test_spins[4] = 0.1;
    test_spins[2] = 5.2;
    test_spins[5] = 1.1;

    std::valarray<double> grads(size*2);
    struct hmc::HamiltonianOptions options;
    std::function<void(std::valarray<double>&, const std::valarray<double>&)> g_func;

    options.J = 1.0;
    options.H = 0;
    g_func = hmc::gen_total_grad( options, 1 );
    g_func(grads, test_spins);
    EXPECT_NEAR(grads[0], -0.58118303, 0.58118303*1e-7);
    EXPECT_NEAR(grads[1], -0.11110539, 0.11110539*1e-7);
    EXPECT_NEAR(grads[2], 0.69228842, 0.69228842*1e-7);
    EXPECT_NEAR(grads[3], -1.2087711, 1.2087711*1e-7);
    EXPECT_NEAR(grads[4], -0.57549375, 0.57549375*1e-7);
    EXPECT_NEAR(grads[5], -0.27048617, 0.27048617*1e-7);

    options.J = 0.0;
    g_func = hmc::gen_total_grad(options, 1);
    g_func(grads, test_spins);
    EXPECT_EQ(grads[0], 0);
    EXPECT_EQ(grads[1], 0);
    EXPECT_EQ(grads[2], 0);
    EXPECT_EQ(grads[3], 0);
    EXPECT_EQ(grads[4], 0);
    EXPECT_EQ(grads[5], 0);
}

#endif
