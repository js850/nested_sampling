import argparse

from nested_sampling.models.harmonic import Harmonic
from nested_sampling.models.harmonic_nowalk import HarmonicSampler
from nested_sampling._nested_sampling import NestedSampling
from nested_sampling.nested_sampling_runner import run_nested_sampling


def main():
    parser = argparse.ArgumentParser(description="do nested sampling on a p[article in a n-dimensional Harmonic well")
#    parser.add_argument("--db", type=str, nargs=1, help="database filename",
#                        default="otp.db")
    parser.add_argument("-K", "--nreplicas", type=int, help="number of replicas", default=300)
    parser.add_argument("-A", "--ndof", type=int, help="number of degrees of freedom", default=4)
    parser.add_argument("-P", "--nproc", type=int, help="number of processors", default=1)
    parser.add_argument("-q", action="store_true", help="turn off verbose printing of information at every step")
    args = parser.parse_args()


    system = Harmonic(args.ndof)
    mcrunner = HarmonicSampler(system, args.ndof)
    ns = NestedSampling(system, args.nreplicas, mcrunner, nproc=args.nproc, 
                        verbose=not args.q)
    print "harmonic particle ndof", args.ndof
    run_nested_sampling(ns, label="hparticle", etol=1e-5)
    

if __name__ == "__main__":
    main()