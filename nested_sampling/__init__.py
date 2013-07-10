"""
.. currentmodule:: nested_sampling

.. autosummary::
   :toctree: generated/

    NestedSampling
    
    run_nested_sampling


"""

from scipy.optimize import Result
from mc_walker import MonteCarloChain
from _nested_sampling import NestedSampling, MCRunner
from nested_sampling_runner import run_nested_sampling
