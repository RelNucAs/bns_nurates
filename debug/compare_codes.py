## Python script to compare emissivities and opacities from Fortran and C++ codes

import numpy as np
import parameters as pars
import units
import math
import matplotlib.pyplot as plt

# Print 
if pars.WM:
    print("Comparing rates corrected by WM + recoil")
else:
    print("comparing rates without WM correction")

# Import Fortran data
f_data = np.nan_to_num(np.loadtxt(pars.f_file,comments='#',unpack=True))
r,d,T,ye  = f_data[1:5]
f_rates = f_data[11:17]
#print(f_rates)
# Import C++ data
C_data = np.nan_to_num(np.loadtxt(pars.C_file,comments='#',unpack=True))
C_rates = C_data[8:14]
#print(C_rates)
C_rates[4] = 2*math.pi*C_rates[4] #possibily 2pi factor more in Fortran code 
C_rates[5] = 2*math.pi*C_rates[5]

rel_diff = (f_rates-C_rates)/np.minimum(abs(f_rates),abs(C_rates))
rel_diff = np.where((C_rates == 0.) & (f_rates == 0.), 0., rel_diff)
print(C_rates[4])
print(C_rates[5])
rel_diff = np.nan_to_num(rel_diff)
#print([np.amax(rel_diff[i]) for i in range(rel_diff.shape[0])])

# Make plots for comparison
y1lab = [r'$j_{\nu_e}$ [1/s]', r'$\chi_{\nu_e}$ [1/cm]', r'$j_{\bar{\nu}_e}$ [1/s]', r'$\chi_{\bar{\nu}_e}$ [1/cm]',
         r'$B^{\rm IS}_{\nu_e}$ [1/cm]', r'$B^{\rm IS}_{\bar{\nu}_e}$ [1/cm]']
y2lab = [r'$|\Delta j|/j$', r'$|\Delta \chi|/\chi$', r'$|\Delta j|/j$', r'$|\Delta \chi|/\chi$',
         r'$\Delta|R^{\rm IS}|/R^{\rm IS}$', r'$\Delta|R^{\rm IS}|/R^{\rm IS}$']
tit = [r'$e^-+p\rightarrow\nu_e+n$:     emissivity', r'$\nu_e+n\rightarrow e^-+p$:  opacity', 
        r'$e^++n\rightarrow\bar{\nu}_e+p$:   emissivity', r'$\bar{\nu}_e+n\rightarrow e^++n$:   opacity', 
        r'$\nu+N\rightarrow\nu+N$:  source term', r'$\bar{\nu}+N\rightarrow\bar{\nu}+N$:  source term'] 

for i in range(rel_diff.shape[0]):
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')

    fig.suptitle(tit[i])

    for j in range(2):
        axs[0][j].set_xscale('log')

    axs[0][0].set_yscale('log')

    axs[1][0].set_xlabel(r'$r$ [cm]')
    axs[1][1].set_xlabel(r'$\rho$ [g/cm$^3$]')

    axs[0][0].set_ylabel(y1lab[i])
    axs[1][0].set_ylabel(y2lab[i])

    axs[1][0].set_ylim((-0.0005, 0.0005))

    # Remove (horizontal) space between axes
    fig.subplots_adjust(hspace=0,wspace=0)

    # Plot each graph, and manually set the y tick values
    axs[0][0].plot(r,abs(f_rates[i]),marker='.',label='Fortran')
    axs[0][0].plot(r,abs(C_rates[i]),marker='.',label='C++')

    axs[0][1].plot(d,abs(f_rates[i]),marker='.')
    axs[0][1].plot(d,abs(C_rates[i]),marker='.')

    axs[1][0].plot(r,rel_diff[i],marker='.')
    axs[1][0].axhline(y=0., color='k', linestyle='--')

    axs[1][1].plot(d,rel_diff[i],marker='.')
    axs[1][1].axhline(y=0., color='k', linestyle='--')

    fig.legend(loc='upper right')

    plt.show()
