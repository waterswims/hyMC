from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import numpy

import os
MKLROOT = os.environ['MKLROOT']

setup(
    name='pyhmc',
    version='0.1dev',
    ext_modules= cythonize(Extension(
        name='pyhmc',
        sources=["pyhmc.pyx"],
        language='c++',
        extra_compile_args=[
            "-std=c++11", '-pthread', '-O3','-fopenmp', '-simd', '-qopenmp', '-xHost',
            '-DUSEMKL', '-DMKL_ILP64', '-I{}/include'.format(MKLROOT)],
        extra_link_args=[
            '-std=c++11', '-pthread', '-fopenmp', '-use-intel-optimized-headers',
            '-Wl,--start-group', '{}/lib/intel64/libmkl_intel_ilp64.a'.format(MKLROOT),
            '{}/lib/intel64/libmkl_sequential.a'.format(MKLROOT),
            '{}/lib/intel64/libmkl_core.a'.format(MKLROOT), '-Wl,--end-group',
            '-lpthread', '-lm', '-ldl'],
        libraries=['hymc', 'm', 'dl', 'pthread'],
        library_dirs=['.', '{}/lib/intel64'.format(MKLROOT)],
        include_dirs=[numpy.get_include(), './include']
    ))
)
