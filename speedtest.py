import pyhmc
import numpy as np

J = 1.
H = 0
kb = 1.
T = 1.0
dims = np.array([10, 10, 10], dtype=np.double)
Nsamp = 10000
eps = 0.01

res = pyhmc.simulate(J, H, kb, T, dims, Nsamp, eps)
