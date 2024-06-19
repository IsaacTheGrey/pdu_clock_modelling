import unittest
import numpy as np
from goodwin_model.peak_analysis import calculate_period_from_peaks

class TestPeakAnalysis(unittest.TestCase):
    
    def test_calculate_period_from_peaks(self):
        dt = 0.01
        # Create a synthetic signal with a known period
        t = np.arange(0.0, 10, dt)
        signal = np.sin(2 * np.pi * t / 1.0)  # period of 1 second
        solution = np.column_stack((signal, np.zeros_like(signal)))

        period, peaks = calculate_period_from_peaks(solution, dt)
        expected_period = 1.0

        self.assertAlmostEqual(period, expected_period, places=2)
        self.assertGreater(len(peaks), 1)

if __name__ == '__main__':
    unittest.main()
