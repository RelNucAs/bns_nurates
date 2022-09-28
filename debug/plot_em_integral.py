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

dU = ""
dU = "_dU"

WM = ""
WM = "_WM"

## Import data
data = np.loadtxt("../output/j_integral_1"+str(WM)+str(dU)+".txt",comments='#',unpack=True)
r    = data[0]
T    = data[1]
mu_e = data[2]
NR_leg = data[3:8]

#data = np.loadtxt("../output/j_integral_2"+str(WM)+str(dU)+".txt",comments='#',unpack=True)
#NR_leg = data[3:8]

# legend =  0: r [cm], 1: T [MeV], 2: mu_e [MeV],    3-7: NR_leg (n=8,16,24,32,100)

## Make plot
fig, axs = plt.subplots(2, 1, sharex='col', sharey='row')

fig.suptitle(r'Integrated emissivity ($\eta = \int d\nu\,4\pi\nu^2 j_{\nu}$)')

axs[0].set_xscale('log')
axs[0].set_yscale('log')

axs[1].set_yscale('log')

axs[1].set_xlabel(r'$r$ [cm]')
axs[0].set_ylabel(r'$\eta$ [MeV$^4$/s]')
axs[1].set_ylabel(r'$|\Delta I|/I$')

# Remove (horizontal) space between axes
fig.subplots_adjust(hspace=0,wspace=0)

MS = 2.

# Plot each graph, and manually set the y tick values
axs[0].plot(r,NR_leg[0],label=r"$n=8$")
axs[0].plot(r,NR_leg[1],label=r"$n=16$")
axs[0].plot(r,NR_leg[2],label=r"$n=24$")
axs[0].plot(r,NR_leg[3],label=r"$n=32$")
axs[0].plot(r,NR_leg[4],label=r"$n=100$")
#axs[0].axvline(x=80, color='k', ls='--')
axs[0].legend() #loc='upper right')

axs[1].plot(r,compare_arrays(NR_leg[0],NR_leg[4]),label=r"$n=8$")
axs[1].plot(r,compare_arrays(NR_leg[1],NR_leg[4]),label=r"$n=16$")
axs[1].plot(r,compare_arrays(NR_leg[2],NR_leg[4]),label=r"$n=24$")
axs[1].plot(r,compare_arrays(NR_leg[3],NR_leg[4]),label=r"$n=32$")
axs[1].legend() #loc='upper right')

plt.savefig("../output/plots/j_integral_1"+str(WM)+str(dU)+".png", dpi=200)

