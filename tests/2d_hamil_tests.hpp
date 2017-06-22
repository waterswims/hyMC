#ifndef HAMIL_2D_TEST
#define HAMIL_2D_TEST

#include "../include/2d_hamils.hpp"
#include "../include/hmc.hpp"
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

TEST(Hamiltonian_2d, trig_lrud)
{
    int tsize = 9;
    double pi = 3.14159265359;
    std::valarray<double> test_spins(tsize*2);
    test_spins = 1.5;
    test_spins[4] = pi / 2;
    test_spins[16] = pi / 2;

    std::vector<std::valarray<double> > test = hmc::trig_lrud(test_spins);
    EXPECT_NEAR(test[0][4], 0, 1e-6);
    EXPECT_NEAR(test[1][5], 0, 1e-6);
    EXPECT_NEAR(test[2][3], 0, 1e-6);
    EXPECT_NEAR(test[3][7], 0, 1e-6);
    EXPECT_NEAR(test[4][1], 0, 1e-6);

    EXPECT_NEAR(test[5][4], 1, 1e-6);
    EXPECT_NEAR(test[6][5], 1, 1e-6);
    EXPECT_NEAR(test[7][3], 1, 1e-6);
    EXPECT_NEAR(test[8][7], 1, 1e-6);
    EXPECT_NEAR(test[9][1], 1, 1e-6);

    EXPECT_NEAR(test[10][7], 0, 1e-6);
    EXPECT_NEAR(test[11][8], 0, 1e-6);
    EXPECT_NEAR(test[12][6], 0, 1e-6);
    EXPECT_NEAR(test[13][1], 0, 1e-6);
    EXPECT_NEAR(test[14][4], 0, 1e-6);

    EXPECT_NEAR(test[15][7], 1, 1e-6);
    EXPECT_NEAR(test[16][8], 1, 1e-6);
    EXPECT_NEAR(test[17][6], 1, 1e-6);
    EXPECT_NEAR(test[18][1], 1, 1e-6);
    EXPECT_NEAR(test[19][4], 1, 1e-6);
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

TEST(Hamiltonian_2d, total_energy)
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

    struct hmc::HamiltonianOptions options;
    double beta = 1.0;
    std::function<double(const std::valarray<double>&)> E_func;

    options.J = 1.0;
    E_func = hmc::gen_total_energy_2d( options, beta );
    EXPECT_NEAR(E_func(test_spins), 3.4939748074272869,
                                        3.4939748074272869*1e-15);

    options.J = 0.0;
    E_func = hmc::gen_total_energy_2d( options, beta );
    EXPECT_EQ(E_func(test_spins), 0);
}

TEST(Hamiltonian_2d, exchange_grad)
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

    std::valarray<double> grad(size*2);
    grad = 0;
    std::vector<std::valarray<double> > g_input = hmc::trig_lrud(test_spins);
    hmc::exchange_grad_2d(grad, g_input);
    EXPECT_NEAR(grad[0], -1.1623660509440559, 1.1623660509440559*1e-15);
    EXPECT_NEAR(grad[1], 0.024380126112986563, 0.024380126112986563*1e-15);
    EXPECT_NEAR(grad[2], -0.19252639008304462, 0.19252639008304462*1e-15);
    EXPECT_NEAR(grad[3], 1.3305123149141138, 1.3305123149141138*1e-15);
    EXPECT_NEAR(grad[4], -2.4175421944456765, 2.4175421944456765*1e-15);
    EXPECT_NEAR(grad[5], 0.24050534024646605, 0.24050534024646605*1e-15);
    EXPECT_NEAR(grad[6], 2.0356793154147046, 2.0356793154147046*1e-15);
    EXPECT_NEAR(grad[7], -2.2735417704169287, 2.2735417704169287*1e-15);
}

TEST(Hamiltonian_2d, total_grad)
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

    std::valarray<double> grad(size*2);
    struct hmc::HamiltonianOptions options;
    std::function<void(std::valarray<double>&, const std::valarray<double>&)> g_func;

    options.J = 1.0;
    g_func = hmc::gen_total_grad_2d( options );
    g_func(grad, test_spins);
    EXPECT_NEAR(grad[0], -1.1623660509440559, 1.1623660509440559*1e-15);
    EXPECT_NEAR(grad[1], 0.024380126112986563, 0.024380126112986563*1e-15);
    EXPECT_NEAR(grad[2], -0.19252639008304462, 0.19252639008304462*1e-15);
    EXPECT_NEAR(grad[3], 1.3305123149141138, 1.3305123149141138*1e-15);
    EXPECT_NEAR(grad[4], -2.4175421944456765, 2.4175421944456765*1e-15);
    EXPECT_NEAR(grad[5], 0.24050534024646605, 0.24050534024646605*1e-15);
    EXPECT_NEAR(grad[6], 2.0356793154147046, 2.0356793154147046*1e-15);
    EXPECT_NEAR(grad[7], -2.2735417704169287, 2.2735417704169287*1e-15);

    options.J = 0.0;
    g_func = hmc::gen_total_grad_2d( options );
    g_func(grad, test_spins);
    EXPECT_EQ(grad[0], 0);
    EXPECT_EQ(grad[1], 0);
    EXPECT_EQ(grad[2], 0);
    EXPECT_EQ(grad[3], 0);
    EXPECT_EQ(grad[4], 0);
    EXPECT_EQ(grad[5], 0);
    EXPECT_EQ(grad[6], 0);
    EXPECT_EQ(grad[7], 0);
}

#endif
