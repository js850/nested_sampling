"""
a routine to run a nested sampling class
"""
import cPickle as pickle
import copy

def save_replicas_to_binary(fout, ns):
    checkpoint = {}
    checkpoint['replicas'] = ns.replicas
    checkpoint['iter_number'] = ns.iter_number
    checkpoint['failed_mc_walks'] = ns.failed_mc_walks
    checkpoint['_mc_niter'] = ns._mc_niter
    pickle.dump( checkpoint, open(fout,"wb"))
    #Emax can then be obtained by get_new_Emax in _nested_sampling
    
def load_replicas_from_binary(fin):
    checkpoint = pickle.load(open(fin, "rb"))
    return checkpoint

def load_checkpoint(fin, ns):
    checkpoint = load_replicas_from_binary(fin)
    ns.replicas = copy.deepcopy(checkpoint['replicas'])
    ns.nreplicas = copy.deepcopy(len(ns.replicas))
    ns.iter_number = copy.deepcopy(checkpoint['iter_number'])
    ns.failed_mc_walks = copy.deepcopy(checkpoint['failed_mc_walks'])
    ns._mc_niter = copy.deepcopy(checkpoint['_mc_niter'])

def remove_energies(fout, Emax):
    temp_file = open(fout,'r')
    temp_energies = []
    for i, energy in enumerate(temp_file):
        if energy > Emax:
            temp_energies.append(energy)
        else:
            break
    temp_file.close()
    temp_file = open(fout, 'w')
    temp_file.writelines(temp_energies)
    temp_file.close()

def write_energies(fout, max_energies, isave=0):
    fout.write( "\n".join([ str(e) for e in max_energies[isave:]]) )
    fout.write("\n")


def run_nested_sampling(ns, label="ns_out", etol=0.01, maxiter=None,
                           iprint_replicas=1000):
    isave = 0
    i = 0
    
    print "nreplicas", len(ns.replicas)
    
    if ns.cpfile == None:
        fout_replicas_name = label + ".replicas.p"
    else:
        fout_replicas_name = ns.cpfile
    
    fout_energies_file = label+".energies"
    
    if ns.cpstart == True:
        fout_energies = open(fout_energies_file, "a")
    else:
        fout_energies = open(fout_energies_file, "w")
    
    while True:
        #start from checkpoint binary file?
        if ns.cpstart == True and i == 0 :
            load_checkpoint(fout_replicas_name, ns)
            remove_energies(fout_energies_file, ns._get_new_Emax())
        
        ediff = ns.replicas[-1].energy - ns.replicas[0].energy

        # save max energies to a file
        if i != 0 and i % 100 == 0:
            write_energies(fout_energies, ns.max_energies, isave=isave)
            isave = len(ns.max_energies)
        
        #pickle current replicas and write them to a file current replicas to a file
        cpfreq = ns.cpfreq
        if i %  cpfreq == 0:
            save_replicas_to_binary(fout_replicas_name, ns)

        if ediff < etol: break
        if maxiter is not None and i >= maxiter: break  
        ns.one_iteration()
        i += 1
    
    write_energies(fout_energies, ns.max_energies, isave=isave)
    fout_energies.close()

    #print_replicas(fout_replicas, ns.replicas)
    #fout_replicas.close()

    # save final replica energies to a file
    # save them with highest energy first
    with open(label+".replicas_final", "w") as fout:
        write_energies(fout, [r.energy for r in reversed(ns.replicas)]) 
    

    print "min replica energy", ns.replicas[0].energy
    print "max replica energy", ns.replicas[-1].energy
    
    return ns

