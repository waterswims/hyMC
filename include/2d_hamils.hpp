#ifndef HAMIL_2D
#define HAMIL_2D

#include "../include/hmc.hpp"

#include <functional>
#include <valarray>
#include <vector>

namespace hmc
{
    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the cosine and sine of the spin angles and then shifts
    ///        them left and up.
    ///
    /// \param data Reference to the valarray where the input data is stored.
    ///             The first half of the array is the ordered list of
    ///             azimuthal angles \f$\theta\f$ and the second half is the
    ///             ordered list of polar angles \f$\phi\f$ of each spin.
    ///
    /// \return Vector of the trigonometric transformations of the
    ///         angles of the spins ordered as: cos(theta),
    ///         cos(theta_left), cos(theta_up), sin(theta),
    ///         sin(theta_left), sin(theta_up), cos(phi),
    ///         cos(phi_left), cos(phi_up), sin(phi), sin(phi_left),
    ///         sin(phi_up)
    ///////////////////////////////////////////////////////////////////////////
    std::vector<std::valarray<double> > trig_left_up(const std::valarray<double> data);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the cosine and sine of the spin angles and then shifts
    ///        them left, right, up and down.
    ///
    /// \param data Reference to the valarray where the input data is stored.
    ///             The first half of the array is the ordered list of
    ///             azimuthal angles \f$\theta\f$ and the second half is the
    ///             ordered list of polar angles \f$\phi\f$ of each spin.
    ///
    /// \return Vector of the trigonometric transformations of the
    ///         angles of the spins ordered as:
    ///         cos(theta): (asis, left, right, up, down),
    ///         sin(theta): (asis, left, right, up, down),
    ///         cos(phi): (asis, left, right, up, down),
    ///         sin(phi): (asis, left, right, up, down)
    ///////////////////////////////////////////////////////////////////////////
    std::vector<std::valarray<double> > trig_lrud(const std::valarray<double> data);
}

#endif
