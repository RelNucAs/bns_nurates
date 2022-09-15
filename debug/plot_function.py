import numpy as np
import matplotlib.pyplot as plt
import py_parameters as pars

# Profile id 
idx = 0

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

c1 = (GF*GF/pi) / (h*c/(2*pi))**4.

c2 = gV*gV+3*gA*gA

print(c*c1*c2)

f_data = np.loadtxt("./nue_em_data.txt") #unpack=True)

def nue_emissivity(E, eta_pn, T, mu_e):
    x = E + Q
    return c * c1 * c2 * eta_pn * x**2. * (1.-(me/x)**2.)**0.5 / (1. + np.exp((x-mu_e)/T))

eta_pn = f_data[idx][0]
T      = f_data[idx][1]
mu_e   = f_data[idx][2]


## Plot function as a function of energy

#define energy array
E = np.logspace(np.log10(1.e-6),np.log10(1.e3),num=1000)

#create subplot
fig, ax = plt.subplots(1, sharex='col', sharey='row')

fig.suptitle('Emissivity of electron neutrinos')

#for j in range(2):
    #axs[j].set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E$ [MeV]')
ax.set_ylabel(r'$j_{\nu}\,(\nu_e)$')

#plot each graph
ax.plot(E,4*pi*E**3. *nue_emissivity(E,eta_pn,T,mu_e),marker='.') #,label='Fortran')
#plt.legend()

plt.savefig(pars.plotfolder+'nue_em_vs_energy_idx_'+str(idx)+'.png',dpi=200)



