import matplotlib.pyplot as plt
import numpy as np

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
