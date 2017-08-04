#ifndef HAMIL_SHARE
#define HAMIL_SHARE

#include "hmc.hpp"
#include <valarray>
#include <vector>
#include <functional>

namespace hmc
{
    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the zeeman energy of a system.
    ///
    /// Calculates the standard zeeman energy \f$ H \sum_{i}
    /// \mathbf{s}_i\f$ where \f$H =\f$ is actually the field strength
    /// multiplied by the moment of an atom.
    ///
    /// \param trig_angles Vector of the trigonometric transformations of the
    ///         angles of the spins ordered as:
    ///         cos(theta): (asis, left(1d), up(2d), forward(3d)),
    ///         sin(theta): (asis, left(1d), up(2d), forward(3d)),
    ///         cos(phi): (asis, left(1d), up(2d), forward(3d)),
    ///         sin(phi): (asis, left(1d), up(2d), forward(3d))
    /// \param H The field strength multiplied by the moment of an atom
    ///////////////////////////////////////////////////////////////////////////
    double zeeman_energy(
        const std::vector<std::valarray<double> >& trig_angles,
        const double H);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the gradient of the zeeman energy of a system.
    ///
    /// Calculates the gradient \f$\partial/\partial \theta\f$ and
    /// \f$\partial/\partial \phi\f$ of the standard zeeman energy \f$ H \sum_{i}
    /// \mathbf{s}_i\f$ where \f$H =\f$ is actually the field strength
    /// multiplied by the moment of an atom.
    ///
    /// \param trig_angles Vector of the trigonometric transformations of the
    ///         angles of the spins ordered as:
    ///         cos(theta): (asis, left(1d), right, up(2d), down, forward(3d)..),
    ///         sin(theta): (asis, left(1d), right, up(2d), down, forward(3d)..),
    ///         cos(phi): (asis, left(1d), right, up(2d), down, forward(3d)..),
    ///         sin(phi): (asis, left(1d), right, up(2d), down, forward(3d)..)
    /// \param grad_out Reference to the valarray where the output gradient
    ///                 will be stored. The gradient is stored as the element
    ///                 wise derivitive of data.
    /// \param H The field strength multiplied by the moment of an atom
    ///////////////////////////////////////////////////////////////////////////
    void zeeman_grad(
        std::valarray<double>& grad_out,
        const std::vector<std::valarray<double> >& trig_angles,
        const double H);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the exchange energy of a system.
    ///
    /// Calculates the standard Heisenberg exchange energy \f$\sum_{ij} J_{ij}
    /// \mathbf{s}_i\cdot\mathbf{s}_j \f$ where \f$J_{ij} =\f$ -1 for
    /// neighbouring spins and 0 for other pairs.
    ///
    /// \param trig_angles Vector of the trigonometric transformations of the
    ///         angles of the spins ordered as:
    ///         cos(theta): (asis, left(1d), up(2d), forward(3d)),
    ///         sin(theta): (asis, left(1d), up(2d), forward(3d)),
    ///         cos(phi): (asis, left(1d), up(2d), forward(3d)),
    ///         sin(phi): (asis, left(1d), up(2d), forward(3d))
    /// \param J Exchange constant
    /// \param d Dimension of the lattice
    ///////////////////////////////////////////////////////////////////////////
    double exchange_energy(
        const std::vector<std::valarray<double> >& trig_angles,
        const double J,
        const int d);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the gradient of the exchange energy of a system.
    ///
    /// Calculates the gradient \f$\partial/\partial \theta\f$ and
    /// \f$\partial/\partial \phi\f$ of the standard Heisenberg exchange energy
    /// \f$\sum_{ij} J_{ij} \mathbf{s}_i\cdot\mathbf{s}_j \f$ where
    /// \f$J_{ij} =\f$ -1 for neighbouring spins and 0 for other pairs.
    ///
    /// \param trig_angles Vector of the trigonometric transformations of the
    ///         angles of the spins ordered as:
    ///         cos(theta): (asis, left(1d), up(2d), forward(3d)),
    ///         sin(theta): (asis, left(1d), up(2d), forward(3d)),
    ///         cos(phi): (asis, left(1d), up(2d), forward(3d)),
    ///         sin(phi): (asis, left(1d), up(2d), forward(3d))
    /// \param grad_out Reference to the valarray where the output gradient
    ///                 will be stored. The gradient is stored as the element
    ///                 wise derivitive of data.
    /// \param J Exchange constant
    /// \param d Dimension of the lattice
    ///////////////////////////////////////////////////////////////////////////
    void exchange_grad(
        std::valarray<double>& grad_out,
        const std::vector<std::valarray<double> >& trig_angles,
        const double J,
        const int d);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Generates the total energy function for a system.
    ///
    /// Returns a std::function object of the form double
    /// f(const std::valarray<double>&) which calculates and returns
    /// the total energy of a system. A reference to the input data is
    /// passed as the parameter.
    ///
    /// \param options A struct containing the Exchange constant and
    ///                           the external field strength
    /// \param beta The relative temperature
    /// \param d The dimension of the lattice
    ///////////////////////////////////////////////////////////////////////////
    std::function<double(const std::valarray<double>&)>gen_total_energy(
        const HamiltonianOptions options,
        const double beta,
        const int d);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Generates the total gradient function for a system.
    ///
    /// Returns a std::function object of the form void
    /// f(std::valarray<double>&, const std::valarray<double>&) which calculates
    /// the total gradient of a system. A reference to the output gradient is
    /// passed as the first parameter and a reference to the input data is
    /// passed as the second parameter.
    ///
    /// \param options A struct containing the Exchange constant and
    ///                           the external field strength
    /// \param d The dimension of the lattice
    ///////////////////////////////////////////////////////////////////////////
    std::function<void(std::valarray<double>&, const std::valarray<double>&)>
    gen_total_grad( const HamiltonianOptions options, const int d);
}

#endif
