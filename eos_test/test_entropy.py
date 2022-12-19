import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

input_table = '../eos_table/eos_electrons_complete_leo.txt'
input_file  = 'test_entropy.txt'

#grid parameters
nne = np.loadtxt(input_table,skiprows=0,max_rows=1,dtype=int)
nt  = np.loadtxt(input_table,skiprows=1,max_rows=1,dtype=int)
ne_array = np.loadtxt(input_table,skiprows=2,max_rows=1,dtype=float)
t_array  = np.loadtxt(input_table,skiprows=3,max_rows=1,dtype=float)

[X, Y] = np.meshgrid(ne_array,t_array)

#import data
s  = np.loadtxt(input_file,unpack=True)

plot_folder = '../output/plots/eos/'

#plot figure
fig = plt.figure(1)
plt.suptitle('Electron entropy: sensitivity')
#cx = ax.contourf(X,Y,s,norm=colors.LogNorm(vmin=np.amin(s),vmax=np.amax(s)),levels=10**np.arange(math.floor(np.log10(np.amin(ref[i]))),math.floor(np.log10(np.amax(ref[i]))),2,dtype=np.float),extend='both')
cx = plt.pcolormesh(X,Y,s,norm=colors.LogNorm(vmin=np.amin(s),vmax=np.amax(s)),shading='gouraud')
plt.colorbar(cx)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$n_e$ [fm$^{-3}$]')
plt.ylabel(r'$T$ [MeV]')

plt.savefig(plot_folder+"test_entropy.png",dpi=200,bbox_inches='tight')
plt.close()


