import numpy as np
import runmc

x = np.zeros(3*31)

mciter = 1000000
stepsize = .01
Emax = 0.
radius = 2.5
xnew = runmc.mc_cython(x, mciter, stepsize, Emax, radius,)

print xnew