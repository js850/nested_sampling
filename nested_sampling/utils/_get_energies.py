import numpy as np

def file_len(fname):
    with open(fname) as f:
        for i, val in enumerate(f):
            pass
    return i + 1

def read_line(fname, eall):
    with open(fname) as f:
        for i, val in enumerate(f):
            eall[i] = val
            
def getfromtext(fname):
    print 'linecounting...'
    linecount = file_len(fname)
    print 'pre-allocating memory...'
    e = np.zeros(linecount)
    print 'reading energies'
    read_line(fname,e)
    return e

def get_energies(fnames, block=False, live=False):
    """read energies from a list of file names, append the live replicas and return as one list
    """
    if len(fnames) == 1:
        eall = getfromtext(fnames[0])
        return eall
    else:
        eall = []
        i = 0
        for fname in fnames:
            e = getfromtext(fname)
            if block is False:
                eall = np.hstack((eall,e))
            else:
                if live is False:
                    eall.append(e)
                else:
                    if i%2 is 0:
                        eall.append(e)
                        print eall
                    else:
                        j = i/2
                        print e
                        eall[j] = np.hstack((eall[j],e))
            i+=1
        if block is False:
            eall.sort(axis=0)
            eall = eall[::-1]
            print eall
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
    