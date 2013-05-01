"""
a routine to run a nested sampling class
"""



def run_nested_sampling(ns, label="nsout", etol=0.01, maxiter=1e100):
    isave = 0
    i = 0

    with open(label+".energies", "w") as fout:
        while True:
            ediff = ns.replicas[-1].energy - ns.replicas[0].energy
            if i != 0 and i % 100 == 0:
                fout.write( "\n".join([ str(e) for e in ns.max_energies[isave:i]]) )
                fout.write("\n")
                fout.flush()
                isave = i

            if i % 1000 == 0:
                if i == 0: openas = "w"
                else: openas = "a"
                with open(label+".replicas", openas) as xout:
                    xout.write("#energy niter from_random\n")
                    for r in ns.replicas:
                        xout.write("%g %d %d\n" % (r.energy, r.niter, int(r.from_random))) 
                    xout.write("\n\n")

            if ediff < etol: break
            if maxiter is not None and i >= maxiter: break  
            ns.one_iteration()
            i += 1
        fout.write( "\n".join([ str(e) for e in ns.max_energies[isave:i]]) )
        fout.write("\n")
        fout.flush()
#    print ns.max_energies
    print "min replica energy", ns.replicas[0].energy
    print "max replica energy", ns.replicas[-1].energy
    ns.finish()
    return ns

