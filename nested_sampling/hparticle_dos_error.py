import argparse
import numpy as np
from scipy.stats import beta

def run_hparticle_ea(K,ndof,P,E_init=1000,E_fin=0.000001,ranerr=False):
    """
    generates a list of energies analytically for the harmonic oscillator provided an initial energy
    and a final energy
    """
    
    expo = float(2)/ndof
    
    if P == 1:
        alpha = float(K)/(K+1)
    else:
        alpha = 1 - float(P)/(K+1)
        
    E_list = []
    
    E = E_init
    E_list.append(E)
    sig_beta = np.sqrt(float(P)*(1+K-P) / ((K+2)*(K+1)**2)) #generalised expression of the variance of compresison for P processors 
    print 'sig_beta',sig_beta
    while E > E_fin:
        if ranerr == True:
            alpha = np.random.beta(K-P+1,P)
            E = float(E) * ((alpha)**expo)
        else:
            E = float(E) * (alpha**expo)
        E_list.append(E)
    
    return np.array(E_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculate the analytical energies given a specific number of replicas for the harmonic oscillator", 
                                     epilog="")
    parser.add_argument("K", type=int, help="number of replicas")
    parser.add_argument("ndof", type=int, help="number of degrees of freedom")
    parser.add_argument("-P", type=int, help="number of cores for parallel run", default=1)
    parser.add_argument("--Einit", type=int, help="initial energy (default=1000)", default=1000)
    parser.add_argument("--Efin", type=int, help="initial energy (default=1E-8)", default=.00000001)
    parser.add_argument("--ranerr", action='store_true', help="random error flag (default False)", default=False)
    args = parser.parse_args()
    print args
    
    E_list = run_hparticle_ea(args.K,args.ndof,args.P,E_init=args.Einit,E_fin=args.Efin,ranerr=args.ranerr)
    
    with open("hparticle_dos_ea.energies", "w") as fout:
        fout.write("#E\n")
        for vals in E_list:
            fout.write("%g\n" % vals)
            