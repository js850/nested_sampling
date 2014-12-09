.. nested_sampling documentation master file, created by
   sphinx-quickstart on Sun Nov 30 12:26:21 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

nested sampling : the Python parallel nested sampling algorithm
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Source code: https://github.com/js850/nested_sampling

Documentation: http://js850.github.io/nested_sampling

Flexible and efficient Python implementation of the nested sampling algorithm.
This implementation is geared towards allowing statistical physicists to use this
method for thermodynamic analysis but is also being used by astrophysicists.

This implementation uses the language of statistical mechanics (partition function, 
phase space, configurations, energy, density of states) rather than the language 
of Bayesian sampling (likelihood, prior, evidence). This is simply for convenience, 
the method is the same.

The package goes beyond the bare implementation of the method providing:

* built-in parallelisation on single computing node (max total number of cpu threads on a single machine)

* built-in Pyro4-based parallelisation by distributed computing, ideal to run calculations on a cluster or across a network

*  ability to save and restart from checkpoint binary files, ideal for very long calculations

* scripts to compute heat capacities and perform error analysis

* integration with the `MCpele  <https://pele-python.github.io/mcpele/>`_ package to implement efficient Monte Carlo walkers.

.. figure:: nested_sampling.jpg
  :align: left
  :scale: 90%
  
  Figure 1: Snapshot of a Nested Sampling iteration on a multimodal surface. Each contour
  line corresponds to a past maximum energy (log-likelihood) constraint. The inner most
  contour line corresponds to the current constraint and the corresponding sample is removed
  (crossed in red) and replaced by walking a copy of a randomly selected replica (walk trajectory
  is the dashed line and the red point is the new configuration that satisfies the tighter constraint).
  The bottom panel shows the fraction of phase space corresponding to each sample and the
  iterative contractions (corresponding to the contours) are shown by the horizontal lines.
  Note that the nested sampling contraction of phase space is constant in the log of phase space
  volume. The full animation can be run from the example folder.

nested sampling has been authored by Stefano Martiniani Jacob D. Stevenson at the University of Cambridge.
The project is publicly available under the GNU general public licence.

Tutorials
-----------
.. toctree::
   :maxdepth: 3

   getting_started
   
Reference
---------

.. toctree::
   :maxdepth: 2
	
   NestedSampling
	
.. toctree::
   :maxdepth: 2
	
   MonteCarloWalker

.. toctree::
   :maxdepth: 2
   
   Models

.. toctree::
   :maxdepth: 2
   
   PyroParallelisation

.. toctree::
   :maxdepth: 2
   
   Utils

.. toctree::
   :maxdepth: 2
   
   Scripts

Modules
+++++++
.. toctree::
   :maxdepth: 1

   nested_sampling

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

