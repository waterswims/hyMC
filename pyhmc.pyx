# distutils: language=c++
import cython
import numpy as np
cimport numpy as np

from libcpp.vector cimport vector

# valarray
cdef extern from "<valarray>" namespace "std" nogil:
    cdef cppclass dvarray "std::valarray<double>":
        dvarray() except+
        dvarray( size_t ) except+
        double& operator[](size_t)

# Options struct
cdef extern from "hmc.hpp" namespace "hmc":
    struct HamiltonianOptions:
        double J
        double H

# declare the Heisenberg model function
cdef extern from "hmc.hpp" namespace "hmc":
    vector[dvarray] heisenberg_model(
        dvarray &energy,
        dvarray &magnetisation,
        const vector[size_t] dims,
        const HamiltonianOptions options,
        const double beta,
        const double leapfrog_eps,
        const size_t nsamples,
        const long initial_state_seed )

# Wrap function
cpdef simulate(
    J, H, KB, T, np.ndarray[double, ndim=1, mode='c'] dimensions,
    nsamples, lf_eps, init_seed=1001):

    cdef dvarray c_energy = dvarray(nsamples)
    cdef dvarray c_magnetisation = dvarray(nsamples)

    cdef vector[size_t] c_dims
    for dim in dimensions:
        c_dims.push_back( dim )

    cdef HamiltonianOptions options
    options.J = J
    options.H = H
    cdef double beta = 1.0 / (KB * T)

    heisenberg_model( c_energy, c_magnetisation, c_dims, options, beta,
                      lf_eps, nsamples, init_seed )

    energy = np.array([c_energy[i] for i in range(nsamples)])
    magnetisation = np.array([c_magnetisation[i] for i in range(nsamples)])


    return {
        'energy': energy,
        'magnetisation': magnetisation
    }
