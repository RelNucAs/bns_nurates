## Python script to compare output of integration routines

## Import modules
import numpy as np
import math
import matplotlib.pyplot as plt

## Define functions for comparison
def compare_arrays(arr1,arr2):
    arr3 = abs(arr1-arr2)/np.minimum(abs(arr1),abs(arr2))
    arr3[np.where((arr1 == np.nan) | (arr2 == np.nan))] = np.nan
    arr3[np.where((arr1 == 0.) & (arr2 == 0.))] = 0.
    return arr3

n_array = [8,16,24,32]

fig_lag_0, axs_lag_0 = plt.subplots(2, 1, sharex='col', sharey='row')
fig_lag_5, axs_lag_5 = plt.subplots(2, 1, sharex='col', sharey='row')
fig_leg,   axs_leg   = plt.subplots(2, 1, sharex='col', sharey='row')

fig_lag_0.suptitle(r'NR Fermi-Dirac integral $F_5(\eta)$, Gauss-Laguerre ($\alpha=0$)')
fig_lag_5.suptitle(r'NR Fermi-Dirac integral $F_5(\eta)$, Gauss-Laguerre ($\alpha=5$)')
fig_leg.suptitle(r'NR Fermi-Dirac integral $F_5(\eta)$, Gauss-Legendre ($t=\eta$)')

figs = [fig_lag_0, fig_lag_5, fig_leg]
ax  = [axs_lag_0, axs_lag_5, axs_leg]

for axs in ax:
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')

    axs[1].set_xscale('log')
    axs[1].set_yscale('log')

    axs[1].set_xlabel(r'$\eta$')

    axs[0].set_ylabel(r'$F_5(\eta)$')
    axs[1].set_ylabel(r'$|\Delta I|/I$')
    
    axs[0].axvline(x=80, color='k', ls='--')
    axs[1].axvline(x=80, color='k', ls='--')

for fig in figs:
    fig.subplots_adjust(hspace=0,wspace=0)

MS = 2.

for n in n_array:

    ## Import data
    data = np.loadtxt("../output/FD_integral_alpha_0_n_"+str(n)+".txt",comments='#',unpack=True)
    eta = data[0]
    result = data[1]
    GSL_lag_0 = data[2]
    NR_lag_0  = data[3]
    NR_leg    = data[4:7]

    alpha = 5
    data = np.loadtxt("../output/FD_integral_alpha_5_n_"+str(n)+".txt",comments='#',unpack=True)
    if (n == 32):
        NR_lag_5 = data[3]
    else:
        NR_lag_5 = data[2]
    # legend =  0: eta,   1: exact result,  2: GSL_leg, 3: NR_lag,  4-7: NR_leg

    dataplot = [NR_lag_0,NR_lag_5,NR_leg[1]]
    ## Make plot
    for idx,axs in enumerate(ax):
        if (n == 8):
            axs[0].plot(eta,result,marker='.',ms=MS,label="'Exact'")
        axs[0].plot(eta,dataplot[idx],marker='.',ms=MS,label=r"n="+str(n))
        axs[1].plot(eta,compare_arrays(result,dataplot[idx]),marker='.',ms=MS,label=r"n="+str(n))

for axs in ax:
    axs[0].legend() #loc='upper right')
    axs[1].legend() #loc='upper right')

for fig in figs:
    fig.subplots_adjust(hspace=0,wspace=0)

fig_lag_0.savefig("../output/plots/FD_integral_lag_0.png", dpi=200)
fig_lag_5.savefig("../output/plots/FD_integral_lag_5.png", dpi=200)
fig_leg.savefig("../output/plots/FD_integral_leg.png", dpi=200)

