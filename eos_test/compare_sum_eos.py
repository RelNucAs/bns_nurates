import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import argparse

parser = argparse.ArgumentParser(description="Choosing methods to compare")
parser.add_argument("-a", dest="method_a", required=True, type=int)
parser.add_argument("-b", dest="method_b", required=True, type=int)
args = parser.parse_args() #pass arguments

if ((args.method_a != 1) and (args.method_a != 2) and (args.method_a != 3)):
    print("-a must be 1, 2 or 3")
    exit()

if ((args.method_b != 1) and (args.method_b != 2) and (args.method_b != 3)):
    print("-b must be 1, 2 or 3")
    exit()

if (args.method_a == args.method_b):
    print("-a and -b must be different")
    exit()

HR = False

if HR:
    res = '_HR'
else:
    res = ''

species = 2
if (species==1):
    sp = 'electrons'
elif (species==2):
    sp = 'muons'

input_table = '../eos_table/'+sp+'/eos_'+sp+'_complete_leo'+res+'.txt'

#grid parameters
nne = np.loadtxt(input_table,skiprows=0,max_rows=1,dtype=int)
nt  = np.loadtxt(input_table,skiprows=1,max_rows=1,dtype=int)

#number density array
ne_edges = np.loadtxt(input_table,skiprows=2,max_rows=1,dtype=float)
ne_array = np.array([0.5*(ne_edges[i]+ne_edges[i+1]) for i in range(ne_edges.size-1)])

#temperature array
t_edges  = np.loadtxt(input_table,skiprows=3,max_rows=1,dtype=float)
t_array = np.array([0.5*(t_edges[i]+t_edges[i+1]) for i in range(t_edges.size-1)])

#eta_file = '../eos_table/eos_electrons_v2.txt'
#eta_ele = np.loadtxt(eta_file,skiprows=4,max_rows=700,unpack=True,dtype=float)

[X, Y] = np.meshgrid(ne_array,t_array)

#import data
ne_rand, t_rand = np.loadtxt(sp+"/input_entries_"+sp+res+".txt",unpack=True)
first  = np.loadtxt(sp+"/first_method_"+sp+res+".txt",unpack=True)
second = np.loadtxt(sp+"/second_method_"+sp+res+".txt",unpack=True)
third  = np.loadtxt(sp+"/third_method_"+sp+res+".txt",unpack=True)
total = np.array([first, second, third])

#reshape arrays
total = np.reshape(total,(3,7,nne-1,nt-1))
ne_rand = np.reshape(ne_rand,(nne-1,nt-1))
t_rand  = np.reshape(t_rand,(nne-1,nt-1))

[mu1, p1, e1, s1, a_p1, a_e1, a_s1] = total[0]
[mu2, p2, e2, s2, a_p2, a_e2, a_s2] = total[1]
[mu3, p3, e3, s3, a_p3, a_e3, a_s3] = total[2]

titles = [r'Chemical potential $\mu_L$', r'Pressure $P_{L}$', r'Int. energy $e_{L}$', r'Entropy $s_{L}$', r'Anti Pressure $P_{\bar{L}}$', r'Anti Int. energy $e_{\bar{L}}$', r'Anti entropy $s_{\bar{L}}$']

id_a = args.method_a
id_b = args.method_b

ref = total[0][1:]
Y1  = total[id_a-1][1:]
Y2  = total[id_b-1][1:]
titles = titles[1:]

ref = np.array([ref[i]+ref[i+3] for i in range(3)])
Y1  = np.array([Y1[i]+Y1[i+3]   for i in range(3)])
Y2  = np.array([Y2[i]+Y2[i+3]   for i in range(3)])

thres = 1.e-10
Y = np.where(Y1==Y2, 1.e-10, abs(Y1-Y2)/Y1)
#Y = np.where(Y1>0., Y, 1.e-10) 
#Y = np.where((Y1>thres) | (Y2>thres), Y, 1.e-10)
#ref = np.where(ref>0.,ref,1.e-50)

for i in range(Y[0].size):
    if (Y[0].flatten()[i]==0.):
        print("%.3e, %.3e, %.3e" %(Y[0].flatten()[i], Y1[0].flatten()[i], Y2[0].flatten()[i]))
plot_folder = '../output/plots/eos/'

#plot figure
fig, axs = plt.subplots(2, 3, sharex='all', sharey='all', figsize=(16,8))
plt.suptitle('Method %d vs method %d' %(id_a,id_b))
for i in range(3):
    ax = axs[0][i]
    cx1 = ax.contourf(ne_rand,t_rand,ref[i],norm=colors.LogNorm(vmin=np.amin(ref[i]),vmax=np.amax(ref[i])),levels=10**np.arange(math.floor(np.log10(np.amin(ref[i]))),math.floor(np.log10(np.amax(ref[i]))),2,dtype=np.float),extend='both')
    plt.colorbar(cx1,ax=ax)
    ax.set_title(titles[i])
    ax = axs[1][i]
    cx2 = ax.contourf(ne_rand,t_rand,Y[i],norm=colors.LogNorm(vmin=1.e-5,vmax=1.e0),levels=np.logspace(-5,0,num=6),extend='both')
    plt.colorbar(cx2,ax=ax)
    ax.set_title('Max diff: %.2e' %np.amax(Y[i]))
    
axs[0][0].set_xscale('log')
axs[0][0].set_yscale('log')
for i in range(3):
    axs[1][i].set_xlabel(r'$n_L$ [fm$^{-3}$]')
for i in range(2):
    axs[i][0].set_ylabel(r'$T$ [MeV]')

plt.subplots_adjust(top=0.85,bottom=0.1,right=0.8,left=0.1,wspace=0.12)
plt.savefig(plot_folder+sp+"/test_eos_sum_"+sp+"_%d_vs_%d" %(id_a,id_b)+res+".png",dpi=200,bbox_inches='tight')
plt.close()

