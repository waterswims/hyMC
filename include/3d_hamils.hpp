#ifndef HAMIL_3D
#define HAMIL_3D

#include <valarray>
#include <vector>

namespace hmc
{
    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the cosine and sine of the spin angles and then shifts
    ///        them left, up and forward.
    ///
    /// \param data Reference to the valarray where the input data is stored.
    ///             The first half of the array is the ordered list of
    ///             azimuthal angles \f$\theta\f$ and the second half is the
    ///             ordered list of polar angles \f$\phi\f$ of each spin.
    ///
    /// \return Vector of the trigonometric transformations of the
    ///         angles of the spins ordered as:
    ///         cos(theta): (asis, left, up, forward),
    ///         sin(theta): (asis, left, up, forward),
    ///         cos(phi): (asis, left, up, forward),
    ///         sin(phi): (asis, left, up, forward)
    ///////////////////////////////////////////////////////////////////////////
    std::vector<std::valarray<double> > trig_luf(const std::valarray<double> data);

    ///////////////////////////////////////////////////////////////////////////
    /// \brief Calculates the cosine and sine of the spin angles and then shifts
    ///        them left, right, up, down, forward and backward.
    ///
    /// \param data Reference to the valarray where the input data is stored.
    ///             The first half of the array is the ordered list of
    ///             azimuthal angles \f$\theta\f$ and the second half is the
    ///             ordered list of polar angles \f$\phi\f$ of each spin.
    ///
    /// \return Vector of the trigonometric transformations of the
    ///         angles of the spins ordered as:
    ///         cos(theta): (asis, left, right, up, down, forward, backward),
    ///         sin(theta): (asis, left, right, up, down, forward, backward),
    ///         cos(phi): (asis, left, right, up, down, forward, backward),
    ///         sin(phi): (asis, left, right, up, down, forward, backward)
    ///////////////////////////////////////////////////////////////////////////
    std::vector<std::valarray<double> > trig_lrudfb(const std::valarray<double> data);
}

#endif
