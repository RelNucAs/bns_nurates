import numpy as np
import matplotlib.pyplot as plt
import py_parameters as pars

## Define constants
Q = 1.2935 #MeV
GF = 8.957e-44 #MeV*cm^3
pi = 3.1415926535898
h = 4.1356943e-21 #Mev*s
MeV = 1.602176634e-6 #conversion MeV to CGS
c = 2.997924562e+10 #cm/s
me = 0.510998928 #MeV
mb = 1.674e-24 #g
gA = 1.23
gV = 1.

#r, d, mu_eq, mu_1, mu_2 = np.loadtxt("../mu_test.txt",delimiter=",",unpack=True)
r, d, T, T_1, T_2 = np.loadtxt("../t_test.txt",delimiter=",",unpack=True)

#create subplot
fig, ax = plt.subplots(1, sharex='col', sharey='row')

fig.suptitle(r'Neutrino chemical potential')

ax.set_xscale('log')
#ax.set_yscale('log')

ax.set_xlabel(r'$r$ [cm]')
ax.set_ylabel(r'$\mu_\nu$ [MeV]')

MS = 2

ax.plot(d,T,marker='.',ms=MS,color='tab:orange',label=r'$\mu_{\rm eq}$')
ax.plot(d,1./T_1,marker='.',ms=MS,color='tab:green',label=r'$\mu_1$')
ax.plot(d,1./T_2,marker='.',ms=MS,color='tab:blue',label=r'$\mu_2$')
#ax.plot(d,T-T_1,marker='.',ms=MS,color='tab:purple',label=r'$\mu_2$') 
ax.legend()
ax.set_ylim((0., 40.))
ax.set_xlim((5.e11, 5.e14))
#ax.hlines(0.,5.e11, 5.e14)
#print(np.amax(abs(mu_1-mu_2)/np.minimum(abs(mu_1),abs(mu_2))))
plt.show()
plt.savefig(pars.plotfolder+'mu_test.png',dpi=200)

