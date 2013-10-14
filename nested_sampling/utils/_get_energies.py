import numpy as np

def get_energies(fnames, block=False, live=False):
    """read energies from a list of file names, append the live replicas and return as one list
    """
    if len(fnames) == 1:
        return np.genfromtxt(fnames[0])
    else:
        eall = []
        i = 0
        for fname in fnames:
            e = np.genfromtxt(fname)
            if block is False:
                eall += e.tolist()
            else:
                if live is False:
                    eall.append(e.tolist())
                else:
                    if i%2 is 0:
                        eall.append(e.tolist())
                    else:
                        j = i/2
                        eall[j] += e.tolist()
            i+=1
        if block is False:
            eall.sort(key=lambda x: -x)
            eall = np.array(eall).flatten()
        return np.array(eall)

# old version of get energies, it does not make use of 
#def get_energies(fnames, block=False):
#    """read energies from a list of file names and return as one list
#    """
#    if len(fnames) == 1:
#        return np.genfromtxt(fnames[0])
#    else:
#        eall = []
#        for fname in fnames:
#            e = np.genfromtxt(fname)
#            if block is False:
#                eall += e.tolist()
#            else:
#                eall.append(e.tolist())
#        if block is False:
#            eall.sort(key=lambda x: -x)
#            eall = np.array(eall).flatten()
#        return eall
    