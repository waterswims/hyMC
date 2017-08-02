#ifndef HAMIL_SHARE
#define HAMIL_SHARE

#include <valarray>
#include <vector>

namespace hmc
{
    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the zeeman energy of a system.
    ///
    /// Calculates the standard Heisenberg exchange energy \f$ H \sum_{i}
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
}

#endif
