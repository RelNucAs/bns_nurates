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
fig_leg_t, axs_leg_t = plt.subplots(2, 1, sharex='col', sharey='row')

fig_lag_0.suptitle(r'NR Fermi-Dirac integral $F_5(\eta)$, Gauss-Laguerre ($\alpha=0$)')
fig_lag_5.suptitle(r'NR Fermi-Dirac integral $F_5(\eta)$, Gauss-Laguerre ($\alpha=5$)')
fig_leg.suptitle(r'NR Fermi-Dirac integral $F_5(\eta)$, Gauss-Legendre ($t=\eta$)')
fig_leg_t.suptitle(r'NR Fermi-Dirac integral $F_5(\eta)$, Gauss-Legendre ($t=f_{\rm max}$)')

figs = [fig_lag_0, fig_lag_5, fig_leg, fig_leg_t]
ax  = [axs_lag_0, axs_lag_5, axs_leg, axs_leg_t]

for axs in ax:
    #axs[0].set_xscale('log')
    axs[0].set_yscale('log')

    #axs[1].set_xscale('log')
    axs[1].set_yscale('log')

    axs[1].set_xlabel(r'$\eta$')

    axs[0].set_ylabel(r'$F_5(\eta)$')
    axs[1].set_ylabel(r'$|\Delta I|/I$')
    
    axs[0].axvline(x=-20, color='k', ls='--')
    axs[1].axvline(x=-20, color='k', ls='--')
    
    axs[0].axvline(x=80, color='k', ls='--')
    axs[1].axvline(x=80, color='k', ls='--')

for fig in figs:
    fig.subplots_adjust(hspace=0,wspace=0)

MS = 2.

for n in n_array:

    ## Import data
    data = np.loadtxt("../output/FD_integral_n_"+str(n)+".txt",comments='#',unpack=True)
    eta = data[0]
    t   = data[1]
    result = data[2]
    GSL_lag = data[3:5]
    NR_lag  = data[5:7]
    NR_leg  = data[7:11]

    # legend =  0: eta, 1: t_split, 2: exact result,    3-4: GSL_leg, 5-6: NR_lag,  7-10: NR_leg

    pos = np.where(eta > 0.)

    dataplot = [NR_lag[0],NR_lag[1],NR_leg[1],NR_leg[3]]
    ## Make plot
    for idx,axs in enumerate(ax):
        axs[0].plot(eta,dataplot[idx],marker='.',ms=MS,label=r"n="+str(n))
        if (idx == 2):
            axs[1].plot(eta[pos],compare_arrays(result[pos],dataplot[idx][pos]),marker='.',ms=MS,label=r"n="+str(n))
        else:
            axs[1].plot(eta,compare_arrays(result,dataplot[idx]),marker='.',ms=MS,label=r"n="+str(n))
        if (n == 32):
            axs[0].plot(eta,result,marker='.',ms=MS,label="'Exact'")


data_100 = np.loadtxt("../output/FD_integral_n_100.txt",comments='#',unpack=True)
eta = data_100[0]
t   = data_100[1]
result = data_100[2]
NR_leg_100 = data_100[3:7]
axs_leg[0].plot(eta[pos],NR_leg_100[1][pos],marker='.',ms=MS,label=r"n=100")
axs_leg[1].plot(eta[pos],compare_arrays(result[pos],NR_leg_100[1][pos]),marker='.',ms=MS,label=r"n=100")
axs_leg_t[0].plot(eta,NR_leg_100[3],marker='.',ms=MS,label=r"n=100")
axs_leg_t[1].plot(eta,compare_arrays(result,NR_leg_100[3]),marker='.',ms=MS,label=r"n=100")

for axs in ax:
    axs[0].legend() #loc='upper right')
    axs[1].legend() #loc='upper right')

fig_lag_0.savefig("../output/plots/FD_integral_lag_0.png", dpi=200)
fig_lag_5.savefig("../output/plots/FD_integral_lag_5.png", dpi=200)
fig_leg.savefig("../output/plots/FD_integral_leg.png", dpi=200)
fig_leg_t.savefig("../output/plots/FD_integral_leg_t.png", dpi=200)

