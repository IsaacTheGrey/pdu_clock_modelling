import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import periodogram

def goodwin_model_simplified(y0, t, parameters):
    """
    Defines ODEs of the Goodwin model
    """
    X, Y, Z, R, S, C, W = y0

    nu1  = parameters['nu1']
    nu2  = parameters['nu2']
    nu3  = parameters['nu3']
    nu4  = parameters['nu4']
    nu5  = parameters['nu5']
    nu6  = parameters['nu6']
    nu7  = parameters['nu7']
    nu8  = parameters['nu8']
    nu9  = parameters['nu9']
    nu10 = parameters['nu10']
    nu11 = parameters['nu11']
    nu12 = parameters['nu12']
    nu13 = parameters['nu13']
    nu14 = parameters['nu14']

    K1  = parameters['K1']
    K2  = parameters['K2']
    K3  = parameters['K3']
    K4  = parameters['K4']
    K5  = parameters['K5']
    K6  = parameters['K6']
    K7  = parameters['K7']
    K8  = parameters['K8']
    K9  = parameters['K9']
    K10 = parameters['K10']

    hill = parameters['hill']
    hill_S = parameters['hill_S']
    hill_W = parameters['hill_W']

    b = parameters['b']
    c = parameters['c']
    d = parameters['d']

    PFL = b + c * X + d * W  # Positive feedback loop (PFL) term

    dXdt = nu1 * ((K1**hill) / (K1**hill + Z**hill + W**hill_W + S**hill_S)) * PFL - nu2 * (X / (K2 + X))  # clk/bmal
    dYdt = nu3 * X * ((K3**hill) / (K3**hill + W**hill_W)) - nu4 * (Y / (K4 + Y))  # per/tr-cry mRNA and protein
    dZdt = nu5 * Y - nu6 * (Z / (K6 + Z))
    dRdt = nu7 * X - nu8 * (R / (K7 + R))  # rev-erb mRNA and protein
    dSdt = nu9 * R - nu10 * (S / (K8 + S))
    dCdt = nu11 * X * ((K5**hill) / (K5**hill + W**hill_W)) - nu12 * (C / (K9 + C))  # cwo mRNA and protein
    dWdt = nu13 * C - nu14 * (W / (K10 + W))

    return [dXdt, dYdt, dZdt, dRdt, dSdt, dCdt, dWdt]

def integrate_model(y0, t, parameters):
    """
    Integrates the ODEs numerically.
    """
    return odeint(goodwin_model_simplified, y0, t, args=(parameters,))

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

def plot_results(t_asymp, solution, periods, power, dominant_period):
    """
    Plots the results: time series, phase plot, and periodogram.
    """
    fig1 = plt.figure(figsize=(18, 8))

    ax1 = fig1.add_subplot(131)
    ax1.plot(t_asymp, solution[:, 0], label='$bmal$')
    ax1.plot(t_asymp, solution[:, 2], label='PER')
    ax1.plot(t_asymp, solution[:, 4], label='REV-ERB')
    ax1.plot(t_asymp, solution[:, 5], label='$cwo$')
    ax1.plot(t_asymp, solution[:, 6], label='CWO')
    ax1.legend(fontsize=12, loc='upper right')
    ax1.set_xlabel('time [h]', fontsize=14)
    ax1.set_ylabel('concentration [a.u.]', fontsize=14)
    ax1.set_xticks(np.arange(0, t_asymp[-1], step=24))
    ax1.tick_params(axis='both', which='major', labelsize=12)

    ax2 = fig1.add_subplot(132)
    ax2.plot(solution[:, 0], solution[:, 2], color='grey')
    ax2.set_xlabel('CLK/BMAL concentration [a.u.]', fontsize=14)
    ax2.set_ylabel('PER/tr-CRY concentration [a.u.]', fontsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=12)

    ax3 = fig1.add_subplot(133)
    ax3.plot(periods, power)
    ax3.set_xlim(0, 36)
    ax3.set_xlabel('Period [h]', fontsize=14)
    ax3.set_ylabel('Power', fontsize=14)
    ax3.tick_params(axis='both', which='major', labelsize=12)
    ax3.set_xticks(np.arange(0, 36, step=4))

    fig1.subplots_adjust(
        top=0.88,
        bottom=0.22,
        left=0.085,
        right=0.98,
        hspace=1,
        wspace=.2
    )

    plt.show()

    print(f"The dominant period is approximately {dominant_period:.2f} hours.")

# Main function to run the simulation
def main():
    dt = 0.01
    t = np.arange(0.0, 5000, dt)
    y0 = [1, 0, 0, 0, 0, 0, 0]

    # Define dictionary of model parameters
    pars = {
        'nu1': 0.7, 'nu2': 0.5, 'nu3': 0.40, 'nu4': 0.3, 'nu5': 0.70,
        'nu6': 0.35, 'nu7': 0.3, 'nu8': 0.2, 'nu9': 0.1, 'nu10': 0.2,
        'nu11': 0.03, 'nu12': 0.05, 'nu13': 0.3, 'nu14': 0.2, 'K1': 1.00,
        'K2': 1.00, 'K3': 1.00, 'K4': 1.00, 'K5': 0.80, 'K6': 1.00,
        'K7': 1.00, 'K8': 1.00, 'K9': 1.00, 'K10': 1.00, 'hill': 3,
        'hill_S': 1.00, 'hill_W': 1.00, 'b': 1.0, 'c': 1.0, 'd': 0.3
    }

    sol_PFL = integrate_model(y0, t, pars)

    # Remove transient behavior, keep only asymptotic dynamics
    last_osc = 48
    t_asymp = np.arange(0, last_osc, dt)
    sol_PFL_asymp = sol_PFL[-int(last_osc / dt):, :]

    sol_PFL_asymp = normalize_oscillations(sol_PFL_asymp)

    frequencies, power = compute_periodogram(sol_PFL_asymp, dt)

    # Convert frequency to period in hours
    periods = 1 / frequencies

    # Find the dominant period
    dominant_period_idx = np.argmax(power)
    dominant_period = periods[dominant_period_idx]

    plot_results(t_asymp, sol_PFL_asymp, periods, power, dominant_period)

if __name__ == "__main__":
    main()

