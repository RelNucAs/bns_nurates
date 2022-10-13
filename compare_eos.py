import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#grid parameters
nne = 700
nt = 150

amu_g = 1.66054e-24
yemin = 5.e-3
yemax = 5.5e-1
denmin = 5.e+3
denmax = 2.e+16
ne_min = (yemin*denmin)/(amu_g*1.e39)
ne_max = (yemax*denmax)/(amu_g*1.e39)
t_min = 5.e-02
t_max = 1.05e+2

#number density array
ne_edges = np.logspace(np.log10(ne_min),np.log10(ne_max),num=nne+1)
ne_array = np.array([0.5*(ne_edges[i]+ne_edges[i+1]) for i in range(nne)])

#temperature array
t_edges = np.logspace(np.log10(t_min) ,np.log10(t_max),num=nt+1)
t_array = np.array([0.5*(t_edges[i]+t_edges[i+1]) for i in range(nt)])

#import data
ne_rand, t_rand = np.loadtxt("input_entries.txt",unpack=True)
Y1, Y2, Y3 = np.loadtxt("electron_energy.txt",unpack=True)

#reshape arrays
Y1 = np.reshape(Y1,(nt-1,nne-1))
Y2 = np.reshape(Y2,(nt-1,nne-1))
Y3 = np.reshape(Y3,(nt-1,nne-1))
ne_rand = np.reshape(ne_rand,(nt-1,nne-1))
t_rand  = np.reshape(t_rand,(nt-1,nne-1))

x1 = Y1
x2 = Y2

Y = abs(x1-x2)/np.minimum(x1,x2)

nne_cut = 70
nt_cut  = 75

#define mask
mask_ne = np.append(np.arange(0,nne,nne/nne_cut,dtype=int),nne-1)
mask_t  = np.append(np.arange(0,nt ,nt/nt_cut  ,dtype=int),nt-1)

#apply mask
ne_cut = ne_edges[mask_ne]
t_cut  = t_edges[mask_t]

Y_cut = np.zeros((nt_cut,nne_cut))

for i in range(nt_cut):
    for j in range(nne_cut):
        idx = int(0.5*(mask_t[i]+mask_t[i+1]))
        jdx = int(0.5*(mask_ne[j]+mask_ne[j+1]))
        Y_cut[i,j] = Y[idx,jdx]

#plot figure
fig = plt.figure()
#plt.pcolormesh(ne_cut,t_cut,Y_cut)
plt.pcolormesh(ne_rand,t_rand,Y,norm=colors.LogNorm(vmin=1.e-5, vmax=Y.max()),shading='gouraud')
plt.xlabel(r'$n_e$ [fm$^{-3}$]')
plt.ylabel(r'$T$ [MeV]')
plt.xscale('log')
plt.yscale('log')
plt.colorbar()
#plt.show()
plt.savefig("test_eos_tris.png")
