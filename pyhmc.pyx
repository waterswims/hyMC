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
        const vector[int] dims,
        const HamiltonianOptions options,
        const double beta,
        const double leapfrog_eps,
        const int nsamples,
        const int initial_state_seed )

# Wrap function
cpdef simulate(
    double J, double H, double KB, double T,
    long [:] dimensions,
    int nsamples, double lf_eps, int init_seed=1001):

    cdef dvarray c_energy = dvarray(nsamples)
    cdef dvarray c_magnetisation = dvarray(nsamples)

    cdef vector[int] c_dims
    for dim in dimensions:
        c_dims.push_back( dim )

    cdef HamiltonianOptions options
    options.J = J
    options.H = H
    cdef double beta = 1.0 / (KB * T)
    cdef double c_eps = lf_eps
    cdef int c_samp = nsamples
    cdef int c_seed = init_seed

    heisenberg_model( c_energy, c_magnetisation, c_dims, options, beta,
                      c_eps, c_samp, c_seed )

    energy = np.array([c_energy[i] for i in range(nsamples)])
    magnetisation = np.array([c_magnetisation[i] for i in range(nsamples)])

    return {
        'energy': energy,
        'magnetisation': magnetisation
    }
