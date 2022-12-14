import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

HR = True

#grid parameters
nne = 700
nt = 150

if HR:
    res = '_HR'
    nne = nne*2
    nt  = nt*2
else:
    res = ''

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
ne_rand, t_rand = np.loadtxt("input_entries"+res+".txt",unpack=True)
first  = np.loadtxt("first_method"+res+".txt",unpack=True)
second = np.loadtxt("second_method"+res+".txt",unpack=True)
third  = np.loadtxt("third_method"+res+".txt",unpack=True)

#reshape arrays
first  = [np.reshape(therm,(nne-1,nt-1)) for therm in first]
second = [np.reshape(therm,(nne-1,nt-1)) for therm in second]
third  = [np.reshape(therm,(nne-1,nt-1)) for therm in third]

[mu1, p1, e1, s1, a_p1, a_e1, a_s1] = first
[mu2, p2, e2, s2, a_p2, a_e2, a_s2] = second
[mu3, p3, e3, s3, a_p3, a_e3, a_s3] = third

ne_rand = np.reshape(ne_rand,(nt-1,nne-1))
t_rand  = np.reshape(t_rand,(nt-1,nne-1))

titles = [r'Chemical potential $\mu_e$', r'Electron pressure $P_{e^-}$', r'Electron internal energy $e_{e^-}$', r'Electron entropy $s_{e^-}$', r'Positron pressure $P_{e^+}$', r'Positron internal energy $e^+$', r'Positron entropy $s_{e^+}$']

total = [first, second, third]

#plot figure
for i in range(1,len(first)):
    fig = plt.figure(i-1)
    plt.suptitle(titles[i])
    plt.pcolormesh(ne_rand,t_rand,first[i],norm=colors.LogNorm(vmin=1.e-5, vmax=np.amax(first[i])),shading='gouraud')
    plt.xlabel(r'$n_e$ [fm$^{-3}$]')
    plt.ylabel(r'$T$ [MeV]')
    plt.xscale('log')
    plt.yscale('log')
    plt.colorbar()
    plt.savefig("plot_eos_first_"+res+"%d.png" %i)
    plt.close()
