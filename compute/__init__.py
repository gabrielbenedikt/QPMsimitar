"""
compute — QPMsimitar compute engine package.

Provides the backend abstraction for running QPM simulations either
locally (in-process) or remotely (over the network).
"""

from compute.backend import ComputeBackend, LocalBackend
