import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

complete_file = '../eos_table/electrons/eos_electrons_complete_leo.txt'
ne = np.loadtxt(complete_file,skiprows=0,max_rows=1,dtype=int)
nt  = np.loadtxt(complete_file,skiprows=1,max_rows=1,dtype=int)
ne_array = np.loadtxt(complete_file,skiprows=2,max_rows=1,dtype=float)
t_array  = np.loadtxt(complete_file,skiprows=3,max_rows=1,dtype=float)
s   = np.loadtxt(complete_file,skiprows=4+700*5,max_rows=700,unpack=True,dtype=float)

pi = 3.1415926535898
MeV = 1.602176634e-6
mLep = 0.510998928
h = 4.1356943e-21
c = 2.997924562e+10
kB = 8.617333262e-11
K = 8.*np.sqrt(2.)*pi*(mLep/(h*c))**3.

theta = t_array/mLep
for i in range(s.shape[1]):
    s[:,i] = s[:,i]/(K*kB*MeV*theta**1.5)
fig = plt.figure(1)
plt.suptitle('Electron entropy on-the-fly')
plt.pcolormesh(ne_array,t_array,s,norm=colors.LogNorm(vmin=np.amin(s), vmax=np.amax(s)),shading='gouraud')
plt.xlabel(r'$n_e$ [fm$^{-3}$]')
plt.ylabel(r'$T$ [MeV]')
plt.xscale('log')
plt.yscale('log')
plt.colorbar()
plt.savefig("../output/plots/eos/electrons/plot_eos_entropy.png")
plt.close()
