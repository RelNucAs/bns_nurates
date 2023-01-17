import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

kB = 8.617333262e-11 #MeV/K
me = 9.1093837e-28 #g
plot_folder = '../output/plots/eos/'

complete_file = '../eos_table/eos_electrons_complete_leo.txt'
ne = np.loadtxt(complete_file,skiprows=0,max_rows=1,dtype=int)
nt  = np.loadtxt(complete_file,skiprows=1,max_rows=1,dtype=int)
ne_array = np.loadtxt(complete_file,skiprows=2,max_rows=1,dtype=float)
t_array  = np.loadtxt(complete_file,skiprows=3,max_rows=1,dtype=float)
mu  = np.loadtxt(complete_file,skiprows=4,max_rows=700,unpack=True,dtype=float)
P   = np.loadtxt(complete_file,skiprows=4+700*1,max_rows=700,unpack=True,dtype=float)
a_P = np.loadtxt(complete_file,skiprows=4+700*2,max_rows=700,unpack=True,dtype=float)
e   = np.loadtxt(complete_file,skiprows=4+700*3,max_rows=700,unpack=True,dtype=float)
a_e = np.loadtxt(complete_file,skiprows=4+700*4,max_rows=700,unpack=True,dtype=float)
s   = np.loadtxt(complete_file,skiprows=4+700*5,max_rows=700,unpack=True,dtype=float)
a_s = np.loadtxt(complete_file,skiprows=4+700*6,max_rows=700,unpack=True,dtype=float)

therm_leo = np.array([P, e, s, a_P, a_e, a_s])
therm_leo = np.where(therm_leo>0.,therm_leo,1.e-50)

d_array = ne_array*me*1.e39*1.e4
T = 10**np.linspace(8,10,num=3)

fig = plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.xlim((1.e-4,1.e12))
plt.ylim((1.e8,1.e31))
for i in range(T.size):
    print(find_nearest(t_array,T[i]*kB))
    plt.plot(d_array,P[find_nearest(t_array,T[i]*kB),:],label=r'$T=%.0e$ GK' %T[i])
plt.legend()
plt.show()


