"""
a routine to run a nested sampling class
"""


def print_replicas(fout, replicas):
    fout.write("#energy niter from_random\n")
    for r in replicas:
        fout.write("%g %d %d\n" % (r.energy, r.niter, int(r.from_random))) 
    fout.write("\n\n")

def write_energies(fout, max_energies, isave=0):
    fout.write( "\n".join([ str(e) for e in max_energies[isave:]]) )
    fout.write("\n")


def run_nested_sampling(ns, label="ns_out", etol=0.01, maxiter=None,
                           iprint_replicas=1000):
    isave = 0
    i = 0
    
    print "nreplicas", len(ns.replicas)
    
    fout_replicas = open(label+".replicas", "w")
    fout_energies = open(label+".energies", "w")

    while True:
        ediff = ns.replicas[-1].energy - ns.replicas[0].energy

        # save max energies to a file
        if i != 0 and i % 100 == 0:
            write_energies(fout_energies, ns.max_energies, isave=isave)
            isave = len(ns.max_energies)

        # write the current replicas to a file
        if i % 1000 == 0:
            print_replicas(fout_replicas, ns.replicas)

        if ediff < etol: break
        if maxiter is not None and i >= maxiter: break  
        ns.one_iteration()
        i += 1
    
    write_energies(fout_energies, ns.max_energies, isave=isave)
    fout_energies.close()

    print_replicas(fout_replicas, ns.replicas)
    fout_replicas.close()

    # save final replica energies to a file
    # save them with highest energy first
    with open(label+".replicas_final", "w") as fout:
        write_energies(fout, [r.energy for r in reversed(ns.replicas)]) 
    

    print "min replica energy", ns.replicas[0].energy
    print "max replica energy", ns.replicas[-1].energy
    return ns

