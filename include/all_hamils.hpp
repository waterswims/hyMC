#ifndef HAMIL_SHARE
#define HAMIL_SHARE

#include <valarray>
#include <vector>

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
    ///                    angles of the spins ordered as: cos(theta),
    ///                    cos(theta_left), cos(theta_up), sin(theta),
    ///                    sin(theta_left), sin(theta_up), cos(phi),
    ///                    cos(phi_left), cos(phi_up), sin(phi), sin(phi_left),
    ///                    sin(phi_up)
    /// \param H The field strength multiplied by the moment of an atom
    ///////////////////////////////////////////////////////////////////////////
    double zeeman_energy(const std::vector<std::valarray<double> >& trig_angles, double H);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the gradient of the zeeman energy of a system.
    ///
    /// Calculates the gradient \f$\partial/\partial \theta\f$ and
    /// \f$\partial/\partial \phi\f$ of the standard zeeman energy \f$ H \sum_{i}
    /// \mathbf{s}_i\f$ where \f$H =\f$ is actually the field strength
    /// multiplied by the moment of an atom.
    ///
    /// \param trig_angles Vector of the trigonometric transformations of the
    ///                    angles of the spins ordered as: cos(theta),
    ///                    sin(theta), cos(phi), sin(phi), each angle given as
    ///                    current, left, right, up, down.
    /// \param grad_out Reference to the valarray where the output gradient
    ///                 will be stored. The gradient is stored as the element
    ///                 wise derivitive of data.
    /// \param H The field strength multiplied by the moment of an atom
    ///////////////////////////////////////////////////////////////////////////
    void zeeman_grad(std::valarray<double>& grad_out,
        const std::vector<std::valarray<double> >& trig_angles, double H);
}

#endif
