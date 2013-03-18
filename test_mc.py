import numpy as np
from src.runmc import mc_cython

from lj_run import LJClusterNew, MonteCarloCompiled
from nested_sampling import MonteCarloChain
from pygmin.takestep import RandomDisplacement
from pygmin.utils.xyz import write_xyz

system = LJClusterNew(31)


x = system.get_random_configuration()

with open("test.xyz", "w") as fout:
    write_xyz(fout, x)


mciter = 100000
stepsize = .01
Emax = 1e20
radius = 2.5

mcc = MonteCarloCompiled(system, radius)
mcc(x.copy(), mciter, stepsize, Emax)
print mcc.naccept, mcc.nsteps, mcc.energy


takestep = RandomDisplacement(stepsize=stepsize)
mc = MonteCarloChain(system.get_potential(), x.copy(), takestep, Emax, system.get_config_tests())
for i in xrange(mciter):
    mc.step()
print mc.naccept, mc.nsteps, mc.energy


