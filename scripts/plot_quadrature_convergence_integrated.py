import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser 

## Plot M1 source coefficients as a function of the
## number of quadrature points

## Parse reaction type(s)
parser = ArgumentParser(description="Thermodynamic point index && reaction type")
parser.add_argument('--reac', dest='reac', required=True)
args = parser.parse_args()

if (args.reac not in ["beta", "pair", "brem", "iso", "inel", "all", "noinel"]):
    print("--reac must be one between 'beta', 'pair', 'brem', 'iso', 'inel', 'all', 'noinel'")
    exit()

def read_data(n_leg):
    folder = "../tests/tests_quadrature_convergence/output/integrated/"
    filename = "bns_nurates_data_" + args.reac + "_quad_" + str(n_leg) + ".txt"
    
    data = np.loadtxt(folder + filename, comments="#", unpack=True)

    return data


n_list = [10, 20, 40, 60, 80, 100]

data = []
for n in n_list:
    data.append(read_data(n))

try:
    folder = "../tests/tests_quadrature_convergence/output/integrated/"
    n_max = n_list[-1]
    filename = "bns_nurates_data_" + args.reac + "_quad_" + str(n_max) + "_backup.txt"
    data_1 = data[-1]
    data_2 = np.loadtxt(folder + filename, comments="#", unpack=True)
    ratio = data_1[6:9] / data_2[6:9]
    print(ratio)

    fig_name = "Ratio for n = " + str(n_max)

    fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]}, num=fig_name)

    # Suptitle
    plt.suptitle(fig_name)

    # Ylabels
    axs[0][0].set_ylabel(r"$\eta~[{\rm s}^{-1}]$")
    axs[0][0].set_xscale("log")

    for idx in range(3):
        axs[0][idx].set_yscale("log")
        
        axs[0][idx].plot(data_1[0], data_1[6+idx], marker="", label=r"Data 1")
        axs[0][idx].plot(data_2[0], data_2[6+idx], marker="", label=r"Data 2")
        axs[1][idx].plot(data_1[0], ratio[idx], marker="", label=r"$n=%d$" %n)     

    plt.subplots_adjust(hspace=0.)
    
    # Title
    axs[0][0].set_title(r"$\nu_{\rm e}$")
    axs[0][1].set_title(r"$\bar{\nu}_{\rm e}$")
    axs[0][2].set_title(r"$\nu_{\rm x}$")

    # Ylabel
    axs[1][0].set_ylabel("Ratio")
    axs[1][0].set_xlabel("Radius [cm]")
                

    axs[0][0].legend()
    axs[1][0].legend()
except Exception as e:
    print(e)
    print("Backup data not found...")

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
