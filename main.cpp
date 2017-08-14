#include "include/hmc.hpp"
#include <valarray>

int main()
{
    int N = 10000;
    std::valarray<double> mag(N), E(N);
    std::vector<size_t> dims = {10, 10, 10};
    hmc::HamiltonianOptions opt;
    opt.H = 0;
    opt.J = 1;
    double beta = 1;
    double eps = 0.01;

    hmc::heisenberg_model(E, mag, dims, opt, beta, eps, N, 1000);
}
