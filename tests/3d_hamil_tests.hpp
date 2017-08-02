#ifndef HAMIL_3D_TEST
#define HAMIL_3D_TEST

#include "../include/3d_hamils.hpp"
#include <gtest/gtest.h>

TEST(Hamiltonian_3d, trig_luf)
{
    int tsize = 27;
    double pi = 3.14159265359;
    std::valarray<double> test_spins(tsize*2);
    test_spins = 1.5;
    test_spins[13] = pi / 2;
    test_spins[27] = pi / 2;

    std::vector<std::valarray<double> > test = hmc::trig_luf(test_spins);
    EXPECT_NEAR(test[0][13], 0, 1e-6);
    EXPECT_NEAR(test[1][14], 0, 1e-6);
    EXPECT_NEAR(test[2][16], 0, 1e-6);
    EXPECT_NEAR(test[3][22], 0, 1e-6);

    EXPECT_NEAR(test[4][13], 1, 1e-6);
    EXPECT_NEAR(test[5][14], 1, 1e-6);
    EXPECT_NEAR(test[6][16], 1, 1e-6);
    EXPECT_NEAR(test[7][22], 1, 1e-6);

    EXPECT_NEAR(test[8][0], 0, 1e-6);
    EXPECT_NEAR(test[9][1], 0, 1e-6);
    EXPECT_NEAR(test[10][3], 0, 1e-6);
    EXPECT_NEAR(test[11][9], 0, 1e-6);

    EXPECT_NEAR(test[12][0], 1, 1e-6);
    EXPECT_NEAR(test[13][1], 1, 1e-6);
    EXPECT_NEAR(test[14][3], 1, 1e-6);
    EXPECT_NEAR(test[15][9], 1, 1e-6);
}

TEST(Hamiltonian_3d, trig_lrudfb)
{
    int tsize = 27;
    double pi = 3.14159265359;
    std::valarray<double> test_spins(tsize*2);
    test_spins = 1.5;
    test_spins[13] = pi / 2;
    test_spins[27] = pi / 2;

    std::vector<std::valarray<double> > test = hmc::trig_lrudfb(test_spins);
    EXPECT_NEAR(test[0][13], 0, 1e-6);
    EXPECT_NEAR(test[1][14], 0, 1e-6);
    EXPECT_NEAR(test[2][12], 0, 1e-6);
    EXPECT_NEAR(test[3][16], 0, 1e-6);
    EXPECT_NEAR(test[4][10], 0, 1e-6);
    EXPECT_NEAR(test[5][22], 0, 1e-6);
    EXPECT_NEAR(test[6][4], 0, 1e-6);

    EXPECT_NEAR(test[7][13], 1, 1e-6);
    EXPECT_NEAR(test[8][14], 1, 1e-6);
    EXPECT_NEAR(test[9][12], 1, 1e-6);
    EXPECT_NEAR(test[10][16], 1, 1e-6);
    EXPECT_NEAR(test[11][10], 1, 1e-6);
    EXPECT_NEAR(test[12][22], 1, 1e-6);
    EXPECT_NEAR(test[13][4], 1, 1e-6);

    EXPECT_NEAR(test[14][0], 0, 1e-6);
    EXPECT_NEAR(test[15][1], 0, 1e-6);
    EXPECT_NEAR(test[16][2], 0, 1e-6);
    EXPECT_NEAR(test[17][3], 0, 1e-6);
    EXPECT_NEAR(test[18][6], 0, 1e-6);
    EXPECT_NEAR(test[19][9], 0, 1e-6);
    EXPECT_NEAR(test[20][18], 0, 1e-6);

    EXPECT_NEAR(test[21][0], 1, 1e-6);
    EXPECT_NEAR(test[22][1], 1, 1e-6);
    EXPECT_NEAR(test[23][2], 1, 1e-6);
    EXPECT_NEAR(test[24][3], 1, 1e-6);
    EXPECT_NEAR(test[25][6], 1, 1e-6);
    EXPECT_NEAR(test[26][9], 1, 1e-6);
    EXPECT_NEAR(test[27][18], 1, 1e-6);
}

#endif
