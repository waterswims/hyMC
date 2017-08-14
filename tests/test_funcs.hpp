#ifndef _TESTFUNCS
#define _TESTFUNCS

#include <vector>
#include <cmath>

///////////////////////////////////////////////////////
// Non-test functions
///////////////////////////////////////////////////////

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

#endif
