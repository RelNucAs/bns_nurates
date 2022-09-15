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

n = 8
alpha = 0

## Import data
data = np.loadtxt("../output/FD_integral_alpha_"+str(alpha)+"_n_"+str(n)+".txt",comments='#',unpack=True)
eta = data[0]
result = data[1]
GSL_lag = data[2]
NR_lag_0  = data[3]
NR_leg  = data[4:7]

alpha = 5
data = np.loadtxt("../output/FD_integral_alpha_"+str(alpha)+"_n_"+str(n)+".txt",comments='#',unpack=True)
if (n == 32):
    NR_lag_5 = data[3]
else:
    NR_lag_5 = data[2]
# legend =  0: eta,   1: exact result,  2: GSL_leg, 3: NR_lag,  4-7: NR_leg

## Make plot

fig, axs = plt.subplots(2, 1, sharex='col', sharey='row')

fig.suptitle(r'NR Fermi-Dirac integral $F_5(\eta)$')

axs[0].set_xscale('log')
axs[0].set_yscale('log')

axs[1].set_xscale('log')
axs[1].set_yscale('log')

axs[1].set_xlabel(r'$\eta$')

axs[0].set_ylabel(r'$F_5(\eta)$')
axs[1].set_ylabel(r'$|\Delta I|/I$')

#axs[1][0].set_ylim((-0.0005, 0.0005))

# Remove (horizontal) space between axes
fig.subplots_adjust(hspace=0,wspace=0)

MS = 2.

# Plot each graph, and manually set the y tick values
axs[0].plot(eta,NR_lag_0,marker='.',ms=MS,label=r"Lag ($\alpha=0$)")
axs[0].plot(eta,NR_lag_5, marker='.',ms=MS,label=r"Lag ($\alpha=5$)")
axs[0].plot(eta,NR_leg[1], marker='.',ms=MS,label=r"Leg ($t=\eta$)")
axs[0].plot(eta,result, marker='.',ms=MS,label=r"'Exact'")
axs[0].axvline(x=80, color='k', ls='--')
axs[0].legend() #loc='upper right')


axs[1].plot(eta,compare_arrays(result,NR_lag_0),marker='.',ms=MS,label=r"Lag ($\alpha=0$)")
axs[1].plot(eta,compare_arrays(result,NR_lag_5),marker='.',ms=MS,label=r"Lag ($\alpha=5$)")
axs[1].plot(eta,compare_arrays(result,NR_leg[0]),marker='.',ms=MS,label=r"Leg ($t=\eta/2$)")
axs[1].plot(eta,compare_arrays(result,NR_leg[1]),marker='.',ms=MS,label=r"Leg ($t=\eta$)")
axs[1].plot(eta,compare_arrays(result,NR_leg[2]),marker='.',ms=MS,label=r"Leg ($t=2\eta$)")
axs[1].axvline(x=80, color='k', ls='--')
axs[1].legend()

#print(compare_arrays(result,NR_lag_0))
print(np.amax(abs(GSL_lag-NR_lag_0)/np.minimum(GSL_lag,NR_lag_0)))

plt.savefig("../output/plots/FD_integral_n_"+str(n)+".png", dpi=200)

