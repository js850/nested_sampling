from nested_sampling.src.cv_trapezoidal import compute_cv_c



def compute_heat_capacity(energies, nreplicas, npar=1, ndof=0, Tmin=.1, Tmax=1., nT=100, live_replicas=False):
    """compute the heat capacity and average energies from the energy levels of a nested sampling run
    
    Parameters
    ----------
    energies: nparray of floats
        the array of energy levels which is the main result of a nested sampling run.  These are stored in
        NestedSampling.max_energies
    nreplicas : int
    npar : int
        the number of processors uses in a the nested sampling run.  This changes the phase space compression factor
        between the energies.
    ndof : int
        number of degrees of freedom in the system being studied.
        This is used to add in the contribution from the momentum degrees of freedom.  
        If this is not passed the heat capacity will be off by an additive constant.
    Tmin, Tmax : float
        minimum and maximum temperatures
    nT : int
        number of temperatures at which to compute the heat capacity
    live_replicas : bool
        if True then the last nreplcas energies are the sorted energies of the live replicas
        at the end of the nested sampling run.
    
    Returns
    -------
    T : array of temperatures
    Cv : array of heat capacities (including the momentum degrees of freedom)
    U : the mean energy
    U : the mean  of the (energy squared)
    """
    T, Cv, U, U2 = compute_cv_c(energies, float(npar), 
                                float(nreplicas), float(Tmin), float(Tmax), int(nT), 
                                float(ndof), bool(live_replicas))
    
    return T, Cv, U, U2
    