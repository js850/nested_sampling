import argparse

from hparticle import HarParticle, HarRunner
from nested_sampling_runner import run_nested_sampling
from nested_sampling import NestedSampling


def main():
    parser = argparse.ArgumentParser(description="do nested sampling with basin sampling for lennard jones clusters")
#    parser.add_argument("--db", type=str, nargs=1, help="database filename",
#                        default="otp.db")
    parser.add_argument("-K", "--nreplicas", type=int, help="number of replicas", default=300)
    parser.add_argument("-A", "--ndof", type=int, help="number of degrees of freedom", default=4)
    parser.add_argument("-P", "--nproc", type=int, help="number of precessors", default=1)
    parser.add_argument("-p", "--trivparal", type=bool, help="set whether to do trivial parallelisation, by default True",default=True)
    parser.add_argument("-q", action="store_true", help="turn off vebose printing of information at every step")
    args = parser.parse_args()


    system = HarParticle(args.ndof, Emax_init=1000.)
    mcrunner = HarRunner(system)
    ns = NestedSampling(system, args.nreplicas, mcrunner, nproc=args.nproc, 
                        triv_paral=args.trivparal, verbose=not args.q)
    print "harmonic particle ndof", args.ndof
    run_nested_sampling(ns, label="hparticle", etol=.000001)
    

if __name__ == "__main__":
    main()