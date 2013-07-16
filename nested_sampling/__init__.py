"""
.. currentmodule:: nested_sampling

.. autosummary::
   :toctree: generated/

    NestedSampling
    
    run_nested_sampling


"""

from utils.result import Result
from mc_walker import MonteCarloWalker
from _nested_sampling import NestedSampling
from nested_sampling_runner import run_nested_sampling

from models.harmonic import Harmonic
