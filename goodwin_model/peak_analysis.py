import numpy as np
from scipy.signal import find_peaks

def calculate_period_from_peaks(solution, dt):
    """
    Calculate the period of oscillation using the time interval between two peaks.
    """
    # I assume the first variable (solution[:, 0], CLK/BMAL) is used for peak analysis
    signal = solution[:, 0]
    
    # Find peaks in the signal
    peaks, _ = find_peaks(signal)
    
    # Calculate the time intervals between consecutive peaks
    peak_intervals = np.diff(peaks) * dt
    
    # Calculate the average period
    if len(peak_intervals) > 0:
        average_period = np.mean(peak_intervals)
    else:
        average_period = np.nan  # Not enough peaks to calculate period
    
    return average_period, peaks
