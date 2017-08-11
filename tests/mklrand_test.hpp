#ifndef _RANDTEST
#define _RANDTEST

#include "../include/mklrand.hpp"
#include "test_funcs.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>

mklrand::mkl_irand st_rand_int(1e5, 1);
mklrand::mkl_drand st_rand_double(1e5, 2);
mklrand::mkl_lnrand rand_ln(0, 0.25, 1e5, 3);
mklrand::mkl_nrand rand_n(2, 1, 1e5, 4);

///////////////////////////////////////////////////////
// Random Number Tests
///////////////////////////////////////////////////////

TEST(Random_Numbers, Integer_Test)
{
    int N_atts = 1e6, N_toss = 100;
    std::vector<int> bins(N_toss/2, 0);
    for(int i=0; i < N_atts; i++)
    {
        int N_zero = 0;
        for (int j=0; j < N_toss; j++)
        {
            if(st_rand_int.gen() == 0){N_zero++;}
        }
        bins[abs(N_zero - N_toss/2)]++;
    }
    std::vector<double> expect(N_toss/2);
    for(int i = 0; i < N_toss/2; i++)
    {
        expect[i] = binomial(N_toss, (N_toss/2)-i) * N_atts / float(pow(2, (N_toss - 1)));
    }
    expect[0] = expect[0] / 2;
    double chi2_test = chi2(bins, expect);
    EXPECT_GT(chi2_test, 0.9);
    EXPECT_LT(chi2_test, 1.3);
}

TEST(Random_Numbers, Double_Test)
{
    int N_bins = 100, N_atts = 1e6;
    std::vector<int> bins(N_bins, 0);
    for(int i = 0; i < N_atts; i++)
    {
        bins[int(st_rand_double.gen()*N_bins)]++;
    }
    std::vector<double> expect(N_bins, (N_atts/N_bins));
    double chi2_test = chi2(bins, expect);
    EXPECT_GT(chi2_test, 0.9);
    EXPECT_LT(chi2_test, 1.3);
}

TEST(Random_Numbers, Ln_Test)
{
    int N_bins = 100, N_atts = 1e6;
    std::vector<int> bins(N_bins, 0);
    for(int i = 0; i < N_atts; i++)
    {
        double b_num = int(rand_ln.gen()*(N_bins/2.0));
        if (b_num >= N_bins) {i--; continue;}
        bins[b_num]++;
    }
    std::vector<double> expect(N_bins);
    for(int i=0; i < N_bins; i++)
    {
        double x = (i+0.5)*2.0/float(N_bins);
        double temp = pow(log(x), 2) / 0.125;
        expect[i] = N_atts * exp(-temp) / (50.08490379760093 * x * 0.25 * pow(2*pi, 0.5));
    }
    double chi2_test = chi2(bins, expect);
    EXPECT_GT(chi2_test, 0.9);
    EXPECT_LT(chi2_test, 1.3);
}

TEST(Random_Numbers, Normal_Test)
{
    int N_bins = 100, N_atts = 1e6;
    std::vector<int> bins(N_bins, 0);
    for(int i = 0; i < N_atts; i++)
    {
        double g = rand_n.gen();
        if (g >= 4 || g < 0) {continue;}
        int b_num = int(g*(N_bins/4.0));
        bins[b_num]++;
    }
    std::vector<double> expect(N_bins);
    for(int i=0; i < N_bins; i++)
    {
        double lower = (i)*4.0/float(N_bins);
        double higher = (i+1)*4.0/float(N_bins);
        double x = (i+0.5)*4.0/float(N_bins);
        expect[i] = -N_atts*0.5*(erf((lower - 2.)/pow(2, 0.5)) - erf((higher - 2.)/pow(2, 0.5)));
    }
    double chi2_test = chi2(bins, expect);
    EXPECT_GT(chi2_test, 0.9);
    EXPECT_LT(chi2_test, 1.3);
}

TEST(Random_Numbers, Read_Checkpoint)
{
    st_rand_double.load("tests/test_files/double_load_test.cp");
    st_rand_double.fill();
    double loaded;
    std::ifstream in;
    in.open("tests/test_files/double_load_test.result");
    for (int i = 0; i < 100; i++)
    {
        in >> loaded;
        EXPECT_NEAR(st_rand_double.gen(), loaded, 1e-6);
    }
    in.close();
}

TEST(Random_Numbers, Save_Checkpoint)
{
    std::vector<double> first(100);
    st_rand_double.save("tests/test_files/double_save_test.cp");
    st_rand_double.fill();
    for(int i = 0; i < 100; i++)
    {
        first[i] = st_rand_double.gen();
    }
    st_rand_double.load("tests/test_files/double_save_test.cp");
    st_rand_double.fill();
    for(int i = 0; i < 100; i++)
    {
        EXPECT_EQ(first[i], st_rand_double.gen());
    }
}

#endif
