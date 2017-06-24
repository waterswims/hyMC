#ifndef _RANDTEST
#define _RANDTEST

#include "../include/stdrand.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>

///////////////////////////////////////////////////////
// Non-test functions
///////////////////////////////////////////////////////

stdrand::std_d_unirand st_rand_double(2);
stdrand::std_normrand rand_n(2, 1, 4);

double chi2(std::vector<int>& count, std::vector<double>& expect)
{
    double ans = 0;
    int skip = 0;
    for(int i=0; i < count.size(); i++)
    {
        if(count[i] == 0)
        {
            skip++;
            continue;
        }
        ans += pow((count[i] - expect[i]), 2) / expect[i];
    }
    return ans / (count.size() - 1 - skip);
}

double beta(int a, int b)
{
    return tgamma(a)*tgamma(b)/(tgamma(a+b));
}

double binomial(int n, int k)
{
    return 1 / (k * beta(k, n-k+1));
}

///////////////////////////////////////////////////////
// Random Number Tests
///////////////////////////////////////////////////////

TEST(Std_Random_Numbers, Double_Test)
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

TEST(Std_Random_Numbers, Normal_Test)
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
        double temp = pow((x - 2), 2)/2.;
        expect[i] = -N_atts*0.5*(erf((lower - 2.)/pow(2, 0.5)) - erf((higher - 2.)/pow(2, 0.5)));
        // expect[i] = 0.01595343685 * N_atts * exp(-temp) / (pow(2*pi, 0.5));
    }
    double chi2_test = chi2(bins, expect);
    EXPECT_GT(chi2_test, 0.9);
    EXPECT_LT(chi2_test, 1.3);
}

#endif
