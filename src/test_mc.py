import numpy as np
import runmc

x = np.zeros(3*4)

mciter = 100
stepsize = .01
Emax = 0.
radius = 2.5
xnew = runmc.mc_cython(x, mciter, stepsize, Emax, radius,)

print xnew