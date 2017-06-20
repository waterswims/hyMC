#ifndef _SPINSCLASS
#define _SPINSCLASS

class spin_lattice_1d
{
public:
    double* thetas;
    double* phis;
    int N;
    spin_lattice_1d(const int length);
    ~spin_lattice_1d();
};

#endif
