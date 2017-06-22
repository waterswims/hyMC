#ifndef _DTYPES
#define _DTYPES

#include <valarray>
#include <array>

namespace hmc
{
    // Basic data type
    class d_type
    {
    public:
        d_type(){}
        ~d_type(){}
        virtual double total_energy(d_type* vels){return 0;}
        virtual void total_energy_grad(d_type* &grad_out){}
        virtual double kinetic_energy(){return 0;}

        // needed to make maths work
        virtual void zero(){}
        virtual void add_to_d1(std::valarray<double> arr){}
        virtual void add_to_d2(std::valarray<double> arr){}
    };

    // General Hberg lattice, for storing velocities, grads etc.
    class hberg_lattice: public d_type
    {
    public:
        std::valarray<double> thetas;
        std::valarray<double> phis;
        size_t N;
        hberg_lattice(){}
        hberg_lattice(const size_t length);
        ~hberg_lattice(){}
        void zero(){thetas = 0; phis = 0;}
        void add_to_d1(std::valarray<double> arr){thetas += arr;}
        void add_to_d2(std::valarray<double> arr){phis += arr;}

        // kinetic energy here since will be the same no matter the dimension of
        // the lattice
        double kinetic_energy();
    };

    // Contains 1D energy methods
    class hberg_lattice_1d : public hberg_lattice
    {
    public:
        std::array<bool, 1> E_flags;
        hberg_lattice_1d() : hberg_lattice() {}
        hberg_lattice_1d(const size_t length, const std::array<bool, 1> E_flags_in)
            : hberg_lattice(length) {E_flags = E_flags_in;}
        ~hberg_lattice_1d(){}
        double total_energy(d_type* vels);
        void total_energy_grad(d_type* grad_out);

        double exchange_energy();
        void exchange_energy_grad(d_type* grad_out);
    };
}

#endif
