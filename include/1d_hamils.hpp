#ifndef HAMIL_1D
#define HAMIL_1D

#include <functional>
#include <valarray>
#include <vector>

namespace hmc
{
    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the gradient of the exchange energy of a 1D system
    ///////////////////////////////////////////////////////////////////////////
    void exchange_grad_1d(std::valarray<double>& grad_out,
                          const std::valarray<double>& data);

    std::function<void(std::valarray<double>&, const std::valarray<double>&)>
        gen_total_grad_1d(std::vector<bool> E_flags);

    double exchange_energy_1d(const std::valarray<double>& data);

    std::function<double(const std::valarray<double>&)>
        gen_total_energy_1d(std::vector<bool> E_flags);
}

#endif
