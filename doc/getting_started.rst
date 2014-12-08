Getting started
===============

Basic elements
--------------

Nested sampling consists of three main elements: 

1. the :class:`NestedSampling` class,
#. the :class:`MonteCarloWalker` class and 
#. the :func:`run_nested_sampling` method. 

The :class:`NestedSampling` class fundamentally requires only three input parameters: 

1. a list of initial samples, which we refer to as *replicas* 
#. a callable of type :class:`MonteCarloWalker` and 
#. the number of cores on which to run the simulation. 

Other optional arguments can be useful but we shall ignore them for the time being. The :class:`MonteCarloWalker`
takes care of walking the replicas under some energy (log-likelihood) constraint, to perform uniform sampling 
(rejection sampling) in phase space. :class:`MonteCarloWalker` takes as parameters:

1. a **potential**, which is a class that through a *get_energy()* function,
   returns an energy for a given a set of coordinates. 
#. a **takestep** method which is some callable that takes a set of coordinates and a step-size,
   displaces them (takes a step) and returns the coordinates, see for instance :func:`random_displace`. 
#. a list of accept_tests, which are additional configurational tests. For instance
   one might want the walkers to stay within a certain spherical container, as well as obeying the energy constraint.
#. one can optionally pass to the :class:`MonteCarloWalker` a list of *events* which is a set of operations
   that one want to do on the new configuration.

Integration with the `MCpele <https://pele-python.github.io/mcpele/>`_ walkers is possible 
(they provide c++ performance through a Python interface!), but this will be addressed in another tutorial.
The :func:`run_nested_sampling` method primarily takes as parameters:

1. an object of type :class:`NestedSampling`
#. a string to label the output
#. a tolerance for the termination of the Nested Sampling iteration. By deafault we terminate the
   run when the difference in energy between the highest and lowest energy replicas is less than etol
#. the maximum number of iterations

Writing a nested sampling runner
--------------------------------

Here we shall cover the basic steps necessary for the implementation of a nested sampling runner to
perform our calculations. In the example/harmonic folder we provide two basic implementation of a 
nested sampling runner with thread based parallelisation (single node) and distributed parallelisation (multiple nodes).
Let us start with the most simple case.

First of all we need to import the modules that we need::

    from nested_sampling import NestedSampling, MonteCarloWalker, Harmonic, run_nested_sampling, Replica

and choose a reasonable set of parameters (note that in the actual example we use the argparse module instead)::

    ndof = 3        #no. of degrees of freedom
    nproc = 4       #no. of cores
    nsteps = 1e3    #no. of Monte Carlo steps per walk 
    nreplicas = 1e3 #no. of initial samples (replicas)
    stepsize = 0.1  #stepsize
    etol = 0.01     #stepsize
    
For the potential we choose the most simple function form, that is a *ndof* dimensional harmonic well,
which we take from the models modules. As we mentioned before each potential needs to have a
*get_energy function* that returns the energy for a given set of coordinates::

    from nested_sampling.utils.rotations import vector_random_uniform_hypersphere

    class Harmonic(object):
        def __init__(self, ndof):
            self.ndim = ndof
        
        def get_energy(self, x):
            assert len(x) == self.ndof
            return 0.5 * x.dot(x)
        
        def get_random_configuration(self, radius=10.):
            """ return a random vector sampled uniformly from within a hypersphere of dimensions self.ndim"""
            x = vector_random_uniform_hypersphere(self.ndof) * radius
            return x
            
Now that we have a potential, we need to construct a potential object and the Monte Carlo runner::
    
    #construct potential (cost function)
    potential = Harmonic(ndof)
    
    #construct Monte Carlo walker
    mc_runner = MonteCarloWalker(potential, mciter=nsteps)
    
We the need to initialise the *nreplicas* samples, we do so by uniformly sampling a set of configurations,
and construct the :class:`NestedSampling` class object::

    #initialise replicas (initial uniformly samples set of configurations)
    replicas = []
    for _ in xrange(nreplicas):
        x = potential.get_random_configuration()
        replicas.append(Replica(x, potential.get_energy(x)))
    
    #construct Nested Sampling object
    ns = NestedSampling(replicas, mc_runner, stepsize=stepsize, nproc=nproc, max_stepsize=10)
    
Finally we can run nested sampling doing::
    
    run_nested_sampling(ns, label="run_hparticle", etol=etol)

which will perform nested sampling on *nproc* cores, on a simple node and with output:

* label.energies (one for each iteration)
* label.replicas_final (live replica energies when NS terminates)

Single node nested sampling
+++++++++++++++++++++++++++

In practice if one were to run the example provided, which makes use of argparse, he would have
to use the following terminal command-line::

    $python examples/harmonic/run_hparticle.py --nreplicas 1e3 --ndof 3 --nprocs 4 --nsteps 1e3 --stepsize 0.1 --etol 0.01

Writing a nested sampling worker
--------------------------------
The nested sampling package allows to run the algorithm on distributed architectures making use of 
the `Pyro4 <https://pythonhosted.org/Pyro4/>`_ library. First of all we need to install Pyro4 and 
add the environment variable::

    $export PYRO_SERIALIZERS_ACCEPTED=serpent,json,marshal,pickle
    
or add it to the .bashrc file if we intend to use it frequently. This parallelisation makes use of a 
**dispatcher** (the middle man) that takes care of dispatching the jobs to the **workers** 
assigned to it by the :func:`run_nested_sampling` function. Workers are very similar to the
nested sampling runners above. First we import the modules that we need::

    from nested_sampling import pyro_worker
    from nested_sampling import MonteCarloWalker, Harmonic

We then construct the potential and the MOnte Carlo objects as above::

    nsteps = 1e3    #no. of Monte Carlo steps per walk 
    ndof = 3        #no. of degrees of freedom    
    potential = Harmonic(ndof)
    mc_runner = MonteCarloWalker(potential, mciter=nsteps)

and inizialise the Pyro worker::

    dispatcher_URI = "PYRO:obj###@17###0:3###8" #address of the dispatcher
    worker_name = None                          #name of worker, when None it's chosen automatically
    host = None                                 #host name, when None found automatically
    port = 0                                    #port number
    server_type = "multiplex"                   #type of server
    
    worker = pyro_worker(dispatcher_URI, mc_runner, worker_name=worker_name, host=host, port=port, server_type=server_type)
    worker._start_worker()

In the actual example we use the argparse package to launch the worker from the command line. In practice the user 
only needs to replace the potential and the Monte Carlo runner to fit his needs.

Distributed nested sampling
+++++++++++++++++++++++++++
We start by initialising a dispatcher by running::

    $python scripts/start_dispatcher.py
    
which will choose a random *dispatcher URI* and default *port* from where it will listen for incoming
communications. One can alternatively specify the server name, the host address, the port number and
the server type (multiplex or threaded), we use the multiplex server by default. From the Pyro4 
documentation we note that 

*"a connected proxy that is unused takes up resources on the server. In the case of the threadpool server type, it locks
up a single thread. If you have too many connected proxies at the same time, the server may run out 
of threads and stops responding. (The multiplex server doesn’t have this particular issue)."*

The dispatcher will also print its URI to a default file name *dispatcher_uri.dat* from where we can
read the address of URI (as well as printing it on the termianl). Let us assume that the randomly
allocated dispatcher URI is::

    $PYRO:obj_fbe65d26b5ed49d7bf3a590bea419a63@888.88.888.888:77777
    
we can then start a worker by doing::

    $python scripts/start_worker.py 3 PYRO:obj_fbe65d26b5ed49d7bf3a590bea419a63@888.88.888.888:77777 -n 1000

where the first positional argument is *ndof*, the second positional argument is the dispatcher URI and the optional
argument *-n 1000* is the number of Monte Carlo steps to perform at each call. We can start as many workers as we like,
although we expect that the dispatcher efficiency decreases as the number of workers increases. 







