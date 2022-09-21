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
data = np.loadtxt("../output/FD_integral_legendre_n_"+str(n)+".txt",comments='#',delimiter=",",unpack=True)
eta = data[0]
t = data[1]
I0 = data[2]
I1 = data[3]

#data_res = np.loadtxt("../output/FD_integral_alpha_0_n_32.txt",comments='#',unpack=True)
#eta = data_res[0]
#result = data_res[1]

## Make plot
fig, axs = plt.subplots(2, 1, sharex='col', sharey='row')

fig.suptitle(r'NR Fermi-Dirac integral $F_5(\eta)$, Gauss_legendre ($t=f_{\rm max}$, $n='+str(n)+'$)')

#axs[0].set_xscale('log')
axs[0].set_yscale('log')

#axs[1].set_xscale('log')
#axs[1].set_yscale('log')

axs[1].set_xlabel(r'$\eta$')

axs[0].set_ylabel(r'$F_5(\eta)$')
axs[1].set_ylabel(r'$I_i/(I_0+I_1)$')

# Remove (horizontal) space between axes
fig.subplots_adjust(hspace=0,wspace=0)

MS = 2.

# Plot each graph, and manually set the y tick values
axs[0].plot(eta,I0,marker='.',ms=MS,label=r"$I_0$")
axs[0].plot(eta,I1,marker='.',ms=MS,label=r"$I_1$")
axs[0].plot(eta,I0+I1,marker='.',ms=MS,label=r"Leg ($t=\eta$)")
#axs[0].plot(eta,result, marker='.',ms=MS,label=r"'Exact'")
axs[0].axvline(x=-20, color='k', ls='--')
axs[0].axvline(x=80, color='k', ls='--')
axs[0].legend() #loc='upper right')


axs[1].plot(eta,I0/(I0+I1),marker='.',ms=MS,label=r"$I_0/I$")
axs[1].plot(eta,I1/(I0+I1),marker='.',ms=MS,label=r"$I_1/I$")
#axs[1].plot(eta,compare_arrays(result,NR_leg[0]),marker='.',ms=MS,label=r"Leg ($t=\eta/2$)")
axs[1].axvline(x=80, color='k', ls='--')
axs[1].legend()

plt.savefig("../output/plots/FD_integral_legendre_n_"+str(n)+".png", dpi=200)

