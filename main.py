import numpy as np
from goodwin_model.integration import integrate_model
from goodwin_model.analysis import normalize_oscillations, compute_periodogram
from goodwin_model.plotting import plot_results
from goodwin_model.peak_analysis import calculate_period_from_peaks
import yaml

def main():
    with open('config.yaml', 'r') as file:
        pars = yaml.safe_load(file)

    dt = 0.01
    t = np.arange(0.0, 5000, dt)
    y0 = [1, 0, 0, 0, 0, 0, 0]

    sol_PFL = integrate_model(y0, t, pars)

    # Remove transient behavior, keep only asymptotic dynamics
    last_osc = 48
    t_asymp = np.arange(0, last_osc, dt)
    sol_PFL_asymp = sol_PFL[-int(last_osc / dt):, :]

    sol_PFL_asymp = normalize_oscillations(sol_PFL_asymp)

    frequencies, power = compute_periodogram(sol_PFL_asymp, dt)

    # Convert frequency to period in hours
    periods = 1 / frequencies

    # Find the dominant period from periodogram
    dominant_period_idx = np.argmax(power)
    dominant_period = periods[dominant_period_idx]

    # Calculate period using peak analysis
    peak_period, peaks = calculate_period_from_peaks(sol_PFL_asymp, dt)

    # Plot results
    plot_results(t_asymp, sol_PFL_asymp, periods, power, dominant_period)

    print(f"The dominant period from periodogram is approximately {dominant_period:.2f} hours.")
    print(f"The average period from peak analysis is approximately {peak_period:.2f} hours.")

if __name__ == "__main__":
    main()
