#ifndef _SPINSCLASS
#define _SPINSCLASS

class spin_lattice_1d
{
public:
    double* thetas;
    double* phis;
    int N;
    bool mem_manage;
    spin_lattice_1d(const int length);
    spin_lattice_1d(const int length, double* thetas, double* phis);
    spin_lattice_1d(const spin_lattice_1d &cop);
    ~spin_lattice_1d();
};

#endif
