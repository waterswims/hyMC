#include "include/hmc.hpp"
#include <valarray>
#include <iostream>
#include <fstream>

int main()
{
    int N;
    std::vector<int> dims = {1, 1, 1};
    hmc::HamiltonianOptions opt;
    double beta;
    double eps;

    std::ifstream ifs;
    ifs.open("testin.txt");
    ifs >> N;
    ifs >> dims[0];
    ifs >> dims[1];
    ifs >> dims[2];
    ifs >> opt.H;
    ifs >> opt.J;
    ifs >> beta;
    ifs >> eps;
    ifs.close();

    std::valarray<double> mag(N), E(N);

    hmc::heisenberg_model(E, mag, dims, opt, beta, eps, N, 1001);
    std::ofstream ofs;
    ofs.open("testout.txt");
    for(int i=0; i < N; i++)
    {
        ofs << E[i] << " " << mag[i] << std::endl;
    }
    ofs.close();
}
