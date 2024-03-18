import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser 

## Plot M1 source coefficients as a function of the
## number of quadrature points

## Parse reaction type(s)
parser = ArgumentParser(description="Thermodynamic point index && reaction type")
parser.add_argument('--reac', dest='reac', required=True)
parser.add_argument('--id'  , dest='id'  , required=True)
args = parser.parse_args()

if (args.reac not in ["beta", "pair", "brem", "iso", "inel", "all", "noinel"]):
    print("--reac must be one between 'beta', 'pair', 'brem', 'iso', 'inel', 'all', 'noinel'")
    exit()

def read_data(n_leg, id_file):
    folder = "../tests/tests_quadrature_convergence/output/integrated/"
    filename = "bns_nurates_data_" + args.reac + "_quad_" + str(n_leg) + "_" + str(id_file) + ".txt"
    
    data = np.loadtxt(folder + filename, comments="#", unpack=True)

    return data


n_list = [10, 20, 40, 60, 80, 100]
n_max = n_list[-1]

data = []
for n in n_list:
    data.append(read_data(n, args.id))


fig_name = "Emissivity comparison for n = " + str(n_max) + " points (" + args.reac + ")"

fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]}, num=fig_name)

# Suptitle
plt.suptitle(fig_name)

# Ylabels
axs[0][0].set_ylabel(r"$\eta~[{\rm s}^{-1}]$")
axs[0][0].set_xscale("log")

for idx in range(3):
    axs[0][idx].set_yscale("log")

plt.subplots_adjust(hspace=0.)
    
# Title
axs[0][0].set_title(r"$\nu_{\rm e}$")
axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
axs[0][2].set_title(r"$\nu_{\rm x}$")

# Ylabel
axs[1][0].set_ylabel("Ratio")
axs[1][0].set_xlabel("Neutrino energy [MeV]")
   
idx_0 = 6
if (args.reac == "iso"):
    idx_0 = 18

data_ref = data[-1]
for id_file in range(4):
    try:
        data_plot = read_data(n_max, id_file)
        #data_plot[3,:] = 2. * data_plot[3,:] + 2. * data_plot[4,:]
        ratio = data_plot[idx_0:idx_0+3] / data_ref[idx_0:idx_0+3]
    
        for idx in range(3):
            axs[0][idx].plot(data_plot[0], data_plot[idx_0+idx], marker="", label=r"Data " + str(id_file))
            axs[1][idx].plot(data_plot[0], ratio[idx], marker="", label=r"Data " + str(id_file) + " / Data " + str(args.id))

    except Exception as e:
            print("Exception for id_file = %d: " %id_file)
            print(e)
            continue

    axs[0][0].legend()
    axs[1][0].legend()

        


ratio = data / data[-1]

if (args.reac != "iso"):

#########################################################################

    fig_name = "Emissivity as a function of quadrature points (" + args.reac + ")"

    fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]}, num=fig_name)

    # Suptitle
    plt.suptitle(fig_name)

    # Ylabels
    axs[0][0].set_ylabel(r"$\eta~[{\rm s}^{-1}]$")
    axs[0][0].set_xscale("log")

    for idx in range(3):
        axs[0][idx].set_yscale("log")

    plt.subplots_adjust(hspace=0.)
    
    # Title
    axs[0][0].set_title(r"$\nu_{\rm e}$")
    axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
    axs[0][2].set_title(r"$\nu_{\rm x}$")

    # Ylabel
    axs[1][0].set_ylabel("Ratio")
    axs[1][0].set_xlabel("Radius [cm]")

    for id_n,n in enumerate(n_list):
        data_plot = data[id_n]
    
        for idx in range(3):
            axs[0][idx].plot(data_plot[0], data_plot[idx+6], marker="", label=r"$n=%d$" %n)
            axs[1][idx].plot(data_plot[0], ratio[id_n][idx+6], marker="", label=r"$n=%d$" %n)

    axs[0][0].legend()
    axs[1][0].legend()
    
    plt.show()
    exit()
    
#########################################################################

    fig_name = "Absoprtion opacity as a function of quadrature points (" + args.reac + ")"
    fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]}, num=fig_name)

    # Suptitle
    plt.suptitle(fig_name)

    # Ylabels
    axs[0][0].set_ylabel(r"$\lambda^{-1}~[{\rm cm}^{-1}]$")
    axs[0][0].set_xscale("log")

    for idx in range(3):
        axs[0][idx].set_yscale("log")

    plt.subplots_adjust(hspace=0.)
    
    # Title
    axs[0][0].set_title(r"$\nu_{\rm e}$")
    axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
    axs[0][2].set_title(r"$\nu_{\rm x}$")

    # Ylabel
    axs[1][0].set_ylabel("Ratio")
    axs[1][0].set_xlabel("Radius [cm]")

    for id_n,n in enumerate(n_list):
        data_plot = data[id_n]
    
        for idx in range(3):
            axs[0][idx].plot(data_plot[0], data_plot[idx+14], marker="", label=r"$n=%d$" %n)
            axs[1][idx].plot(data_plot[0], ratio[id_n][idx+14], marker="", label=r"$n=%d$" %n)

    axs[0][0].legend()
    axs[1][0].legend()

#########################################################################

if ((args.reac == "iso") | (args.reac == "all")):

    fig_name = "Scattering opacity as a function of quadrature points (" + args.reac + ")"
    fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]}, num=fig_name)

    # Suptitle
    plt.suptitle(fig_name)

    # Ylabels
    axs[0][0].set_ylabel(r"$\lambda^{-1}~[{\rm cm}^{-1}]$")
    axs[0][0].set_xscale("log")

    for idx in range(3):
        axs[0][idx].set_yscale("log")

    plt.subplots_adjust(hspace=0.)
    
    # Title
    axs[0][0].set_title(r"$\nu_{\rm e}$")
    axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
    axs[0][2].set_title(r"$\nu_{\rm x}$")

    # Ylabel
    axs[1][0].set_ylabel("Ratio")
    axs[1][0].set_xlabel("Radius [cm]")

    for id_n,n in enumerate(n_list):
        data_plot = data[id_n]
    
        for idx in range(3):
            axs[0][idx].plot(data_plot[0], data_plot[idx+18], marker="", label=r"$n=%d$" %n)
            axs[1][idx].plot(data_plot[0], ratio[id_n][idx+18], marker="", label=r"$n=%d$" %n)

    axs[0][0].legend()
    axs[1][0].legend()

plt.show()
