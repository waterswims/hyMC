#ifndef HAMIL_GEN_TEST
#define HAMIL_GEN_TEST

#include "../include/all_hamils.hpp"

TEST(Hamiltonian_Gen, Zeeman)
{
    double H = 2.1;
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

    EXPECT_NEAR(hmc::zeeman_energy(E_input, H), -0.5827041377, 0.5827041377*1e-9);
}

TEST(Hamiltonian_Gen, Zeeman_Grad)
{
    double H = 2.1;
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
    hmc::zeeman_grad(grad, g_input, H);
    for (int i = 0; i < 4; i++)
    {
        EXPECT_EQ(grad[i], 0);
    }
    EXPECT_NEAR(grad[4], H*0.7457052122, H*0.7457052122*1e-9);
    EXPECT_NEAR(grad[5], H*0.09983341665, H*0.09983341665*1e-9);
    EXPECT_NEAR(grad[6], H*0.8912073601, H*0.8912073601*1e-9);
    EXPECT_NEAR(grad[7], H*0.8632093666, H*0.8632093666*1e-9);
}

#endif
