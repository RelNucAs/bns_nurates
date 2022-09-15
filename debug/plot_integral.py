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


## Import data
data = np.loadtxt("../output/j_integral.txt",comments='#',unpack=True)
r = data[0]
GSL_lag = data[1]
NR_lag  = data[2]
NR_leg  = data[3:8]

# legend =  0: GSL_lag,   1: NR_lag,  2-6: NR_leg

print("Difference between NR and GSL Laguerre:")
tmp = compare_arrays(GSL_lag,NR_lag)
print(tmp)
print("Max difference: %.5e\n" %np.amax(np.nan_to_num(tmp)))

print("Difference between NR Legendre and Laguerre:")
diff = compare_arrays(NR_leg[2],NR_lag)
print(diff)
print("Max difference: %.5e\n" %np.amax(np.nan_to_num(diff)))


## Make plot
fig, ax = plt.subplots(1, sharex='col', sharey='row')

fig.suptitle('Comparison for emissivity integration')

ax.set_xscale('log')
#ax.set_yscale('log')

ax.set_xlabel(r'$r$ [cm]')

ax.set_ylabel(r'$|\Delta I| / I$')

#ax.set_ylim((-0.0005, 0.0005))

ax.plot(r,compare_arrays(NR_leg[0],NR_lag),marker='.',label=r'$a=0.1\eta$')
ax.plot(r,compare_arrays(NR_leg[1],NR_lag),marker='.',label=r'$a=0.5\eta$')
ax.plot(r,compare_arrays(NR_leg[2],NR_lag),marker='.',label=r'$a=\eta$')
ax.plot(r,compare_arrays(NR_leg[3],NR_lag),marker='.',label=r'$a=5\eta$')
ax.plot(r,compare_arrays(NR_leg[4],NR_lag),marker='.',label=r'$a=10\eta$')
ax.legend() #loc='upper right')

plt.savefig("../output/plots/j_integral.png", dpi=200)
