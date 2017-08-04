#ifndef HAMIL_1D
#define HAMIL_1D

#include "../include/hmc.hpp"
#include <functional>
#include <valarray>
#include <vector>

namespace hmc
{
    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the cosine and sine of the spin angles and then shifts
    ///        them left.
    ///
    /// \param data Reference to the valarray where the input data is stored.
    ///             The first half of the array is the ordered list of
    ///             azimuthal angles \f$\theta\f$ and the second half is the
    ///             ordered list of polar angles \f$\phi\f$ of each spin.
    ///
    /// \return Vector of the trigonometric transformations of the
    ///         angles of the spins ordered as:
    ///         cos(theta): (asis, left),
    ///         sin(theta): (asis, left),
    ///         cos(phi): (asis, left),
    ///         sin(phi): (asis, left)
    ///////////////////////////////////////////////////////////////////////////
    std::vector<std::valarray<double> > trig_left(const std::valarray<double> data);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the cosine and sine of the spin angles and then shifts
    ///        them left and right.
    ///
    /// \param data Reference to the valarray where the input data is stored.
    ///             The first half of the array is the ordered list of
    ///             azimuthal angles \f$\theta\f$ and the second half is the
    ///             ordered list of polar angles \f$\phi\f$ of each spin.
    ///
    /// \return Vector of the trigonometric transformations of the
    ///         angles of the spins ordered as:
    ///         cos(theta): (asis, left, right),
    ///         sin(theta): (asis, left, right),
    ///         cos(phi): (asis, left, right),
    ///         sin(phi): (asis, left, right)
    ///////////////////////////////////////////////////////////////////////////
    std::vector<std::valarray<double> > trig_lr(const std::valarray<double> data);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the gradient of the exchange energy of a 1D system.
    ///
    /// Calculates the gradient \f$\partial/\partial \theta\f$ and
    /// \f$\partial/\partial \phi\f$ of the standard Heisenberg exchange energy
    /// \f$\sum_{ij} J_{ij} \mathbf{s}_i\cdot\mathbf{s}_j \f$ where
    /// \f$J_{ij} =\f$ -1 for neighbouring spins and 0 for other pairs.
    ///
    /// \param data Reference to the valarray where the input data is stored.
    ///             The first half of the array is the ordered list of
    ///             azimuthal angles \f$\theta\f$ and the second half is the
    ///             ordered list of polar angles \f$\phi\f$ of each spin.
    /// \param grad_out Reference to the valarray where the output gradient
    ///                 will be stored. The gradient is stored as the element
    ///                 wise derivitive of data.
    ///////////////////////////////////////////////////////////////////////////
    void exchange_grad_1d(std::valarray<double>& grad_out,
                          const std::valarray<double>& data);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Generates the total gradient function for a 1D system.
    ///
    /// Returns a std::function object of the form void
    /// f(std::valarray<double>&, const std::valarray<double>&) which calculates
    /// the total gradient of a 1D system. A reference to the output gradient is
    /// passed as the first parameter and a reference to the input data is
    /// passed as the second parameter.
    ///
    /// \param E_flags A vector of switches which determines the which energy
    ///                gradients to turn on. The first boolean eneables the
    ///                exchange energy gradient.
    ///////////////////////////////////////////////////////////////////////////
    std::function<void(std::valarray<double>&, const std::valarray<double>&)>
    gen_total_grad_1d( const HamiltonianOptions );

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the exchange energy of a 1D system.
    ///
    /// Calculates the standard Heisenberg exchange energy \f$\sum_{ij} J_{ij}
    /// \mathbf{s}_i\cdot\mathbf{s}_j \f$ where \f$J_{ij} =\f$ -1 for
    /// neighbouring spins and 0 for other pairs.
    ///
    /// \param data Reference to the valarray where the input data is stored.
    ///             The first half of the array is the ordered list of
    ///             azimuthal angles \f$\theta\f$ and the second half is the
    ///             ordered list of polar angles \f$\phi\f$ of each spin.
    ///////////////////////////////////////////////////////////////////////////
    double exchange_energy_1d(const std::valarray<double>& data);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Generates the total energy function for a 1D system.
    ///
    /// Returns a std::function object of the form double
    /// f(const std::valarray<double>&) which calculates and returns
    /// the total energy of a 1D system. A reference to the input data is
    /// passed as the parameter.
    ///
    /// \param E_flags A vector of switches which determines the which energy
    ///                terms to turn on. The first boolean eneables the
    ///                exchange energy.
    ///////////////////////////////////////////////////////////////////////////
    std::function<double(const std::valarray<double>&)>
    gen_total_energy_1d( const HamiltonianOptions, const double beta );
}

#endif
