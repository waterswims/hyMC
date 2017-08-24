#ifndef HAMIL_GEN_TEST
#define HAMIL_GEN_TEST

#include "../include/all_hamils.hpp"

TEST(Hamiltonian_Gen, Zeeman)
{
    double H = 2.1;
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
    EXPECT_FLOAT_EQ(hmc::zeeman_energy(H), -0.5827041377);
}

TEST(Hamiltonian_Gen, Zeeman_Grad)
{
    double H = 2.1;
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

    hmc::zeeman_grad(grad, H);
    for (int i = 4; i < 8; i++)
    {
        EXPECT_EQ(grad[i], 0);
    }
    EXPECT_FLOAT_EQ(grad[0], H*0.7457052122);
    EXPECT_FLOAT_EQ(grad[1], H*0.09983341665);
    EXPECT_FLOAT_EQ(grad[2], H*0.8912073601);
    EXPECT_FLOAT_EQ(grad[3], H*0.8632093666);
}

#endif
