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

n = 32

## Import data
data = np.loadtxt("../output/FD_integral_n_"+str(n)+".txt",comments='#',unpack=True)
eta = data[0]
t   = data[1]
result  = data[2]
GSL_lag = data[3:5]
NR_lag  = data[5:7]
NR_leg  = data[7:11]

# legend =  0: eta, 1: t_split, 2: exact result,    3-4: GSL_leg, 5-6: NR_lag,  7-10: NR_leg

pos = np.where(eta>0.)

## Make plot
fig, axs = plt.subplots(2, 1, sharex='col', sharey='row')

fig.suptitle(r'NR Fermi-Dirac integral $F_5(\eta)$')

#axs[0].set_xscale('log')
axs[0].set_yscale('log')

#axs[1].set_xscale('log')
axs[1].set_yscale('log')

axs[1].set_xlabel(r'$\eta$')

axs[0].set_ylabel(r'$F_5(\eta)$')
axs[1].set_ylabel(r'$|\Delta I|/I$')

axs[0].set_xlim((-110,+110))

# Remove (horizontal) space between axes
fig.subplots_adjust(hspace=0,wspace=0)

MS = 2.

# Plot each graph, and manually set the y tick values
axs[0].plot(eta,NR_lag[0],marker='.',ms=MS,label=r"Lag ($\alpha=0$)")
axs[0].plot(eta,NR_lag[1],marker='.',ms=MS,label=r"Lag ($\alpha=5$)")
axs[0].plot(eta[pos],NR_leg[1][pos], marker='.',ms=MS,label=r"Leg ($t=\eta$)")
axs[0].plot(eta,NR_leg[3], marker='.',ms=MS,label=r"Leg ($t=f_{\rm max}$)")
axs[0].plot(eta,result, marker='.',ms=MS,label=r"'Exact'")
axs[0].axvline(x=-20, color='k', ls='--')
axs[0].axvline(x=80, color='k', ls='--')
axs[0].legend() #loc='upper right')


axs[1].plot(eta,compare_arrays(result,NR_lag[0]),marker='.',ms=MS,label=r"Lag ($\alpha=0$)")
axs[1].plot(eta,compare_arrays(result,NR_lag[1]),marker='.',ms=MS,label=r"Lag ($\alpha=5$)")
#axs[1].plot(eta,compare_arrays(result,NR_leg[0]),marker='.',ms=MS,label=r"Leg ($t=\eta/2$)")
axs[1].plot(eta[pos],compare_arrays(result[pos],NR_leg[1][pos]),marker='.',ms=MS,label=r"Leg ($t=\eta$)")
#axs[1].plot(eta,compare_arrays(result,NR_leg[2]),marker='.',ms=MS,label=r"Leg ($t=2\eta$)")
axs[1].plot(eta,compare_arrays(result,NR_leg[3]),marker='.',ms=MS,label=r"Leg ($t=f_{\rm max}$)")
axs[1].axvline(x=-20, color='k', ls='--')
axs[1].axvline(x=80, color='k', ls='--')
axs[1].legend()

#print(compare_arrays(result,NR_lag_0))
print(np.amax(abs(GSL_lag[0]-NR_lag[0])/np.minimum(GSL_lag[0],NR_lag[0])))
print(np.amax(abs(GSL_lag[1]-NR_lag[1])/np.minimum(GSL_lag[1],NR_lag[1])))

plt.savefig("../output/plots/FD_integral_n_"+str(n)+".png", dpi=200)

