import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import argparse

HR = False

if HR:
    res = '_HR'
else:
    res = ''

species = 1

if (species==1):
    sp = 'electrons'
elif (species==2):
    sp = 'muons'

input_table = '../eos_table/'+sp+'/eos_'+sp+'_leo'+res+'.txt'
input_file  = sp+'/test_NR_steps_'+sp+res+'.txt'
#grid parameters
nne = np.loadtxt(input_table,skiprows=0,max_rows=1,dtype=int)
nt  = np.loadtxt(input_table,skiprows=1,max_rows=1,dtype=int)

#number density array
ne_edges = np.loadtxt(input_table,skiprows=2,max_rows=1,dtype=float)
ne_array = np.array([0.5*(ne_edges[i]+ne_edges[i+1]) for i in range(ne_edges.size-1)])

#temperature array
t_edges  = np.loadtxt(input_table,skiprows=3,max_rows=1,dtype=float)
t_array = np.array([0.5*(t_edges[i]+t_edges[i+1]) for i in range(t_edges.size-1)])

#import data
n_step = np.loadtxt(input_file,unpack=True)

#reshape arrays
[n_s1, n_s2] = np.reshape(n_step,(2,nne-1,nt-1))
n_s1 = np.transpose(n_s1)
n_s2 = np.transpose(n_s2)

plot_folder = '../output/plots/eos/'

#plot figure
fig, axs = plt.subplots(1, 3, sharex='all', sharey='all', figsize=(16,4))
plt.suptitle('Number of iterations in Newton-Raphson')

axs[0].set_title('Analytic guess')
cx1 = axs[0].pcolormesh(ne_edges,t_edges,n_s1,shading='flat')
plt.colorbar(cx1,ax=axs[0])

axs[1].set_title('Interpolated guess')
cx2 = axs[1].pcolormesh(ne_edges,t_edges,n_s2,shading='flat')
plt.colorbar(cx2,ax=axs[1])

axs[2].set_title('Ratio guess')
cx3 = axs[2].pcolormesh(ne_edges,t_edges,n_s2/n_s1,shading='flat')
plt.colorbar(cx3,ax=axs[2])

  #  ax = axs[int(i/3)][i%3]
    #cx = ax.pcolormesh(ne_rand,t_rand,ref[i],norm=colors.LogNorm(vmin=1.e-10,vmax=np.amax(ref[i])),shading='gouraud')
        #cx = ax.pcolormesh(ne_rand,t_rand,ref[i],norm=colors.LogNorm(vmin=np.amin(ref[i]),vmax=np.amax(ref[i])),shading='gouraud')
   #     cx = ax.contourf(ne_rand,t_rand,ref[i],norm=colors.LogNorm(vmin=np.amin(ref[i]),vmax=np.amax(ref[i])),levels=10**np.arange(math.floor(np.log10(np.amin(ref[i]))),math.floor(np.log10(np.amax(ref[i]))),2,dtype=np.float),extend='both')
        #cx = ax.pcolormesh(ne_rand,t_rand,ref[i],norm=colors.LogNorm(vmin=1.e-10,vmax=np.amax(ref[i])),shading='gouraud')
    #    cx = ax.contourf(ne_rand,t_rand,ref[i],norm=colors.LogNorm(vmin=1.e-10,vmax=np.amax(ref[i])),levels=10**np.arange(-10,math.floor(np.log10(np.amax(ref[i]))),2,dtype=np.float),extend='both')

axs[0].set_xscale('log')
axs[0].set_yscale('log')

for i in range(3):
    axs[i].set_xlabel(r'$n_L$ [fm$^{-3}$]')
axs[0].set_ylabel(r'$T$ [MeV]')

plt.subplots_adjust(top=0.85,bottom=0.1,right=0.8,left=0.1,wspace=0.12)
plt.savefig(plot_folder+"/"+sp+"/test_NR_steps_"+sp+res+".png",dpi=200,bbox_inches='tight')
plt.close()

