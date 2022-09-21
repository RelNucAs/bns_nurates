import numpy as np
import matplotlib.pyplot as plt
import py_parameters as pars

## Import data
data = np.loadtxt("../output/j_integral.txt",comments='#',unpack=True)
r = data[0]
T = data[1]
mu = data[2]

eta = mu/T

## Make plot
#create subplot
fig, ax = plt.subplots(1, sharex='col', sharey='row')

fig.suptitle('CCSN 1D radial profile')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$r$ [cm]')
ax.set_ylabel(r'[MeV]')

MS=2

#plot each graph
ax.plot(r,T,marker='.',label=r"$T$")
ax.plot(r,mu,marker='.',label=r"$\mu_e$")
ax.axvline(r[0],color='k',ls='--',lw=0.8,label='idx=[0,20,48]')
ax.axvline(r[20],color='k',ls='--',lw=0.8)
ax.axvline(r[48],color='k',ls='--',lw=0.8)
ax2 = ax.twinx()
ax2.plot(r,eta,marker='.',color='tab:green',label=r"$\eta_e$")
ax2.set_yscale('log')

ax.legend()
ax2.legend(loc=(0.855,0.70))

plt.savefig(pars.plotfolder+'eta_radial_profile.png',dpi=200)



