import numpy as np
from scipy.signal import periodogram

def normalize_oscillations(solution):
    """
    Normalizes the oscillations to their mean.
    """
    return solution / solution.mean(axis=0)

def compute_periodogram(solution, dt):
    """
    Computes the periodogram and returns frequencies and power.
    """
    frequencies, power = periodogram(solution[:, 0], 1 / dt)
    nonzero_indices = frequencies != 0
    return frequencies[nonzero_indices], power[nonzero_indices]
