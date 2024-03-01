import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser 

## Plot spectral rates as a function of the
## number of quadrature points

## Parse reaction type(s)
parser = ArgumentParser(description="Thermodynamic point index && reaction type")
parser.add_argument('-p', dest='id_point', type=int, required=True)
parser.add_argument('--reac', dest='reac', required=True)
args = parser.parse_args()

if (args.reac not in ["beta", "pair", "brem", "iso", "inel", "all", "noinel"]):
    print("--reac must be one between 'beta', 'pair', 'brem', 'iso', 'inel', 'all', 'noinel'")
    exit()

def read_data(n_leg):
    folder = "../tests/tests_quadrature_convergence/output/spectral/"
    filename = "bns_nurates_data_" + args.reac + "_point_" + str(args.id_point) + "_quad_" + str(n_leg) + ".txt"
    
    data = np.loadtxt(folder + filename, comments="#", unpack=True)

    return data


n_list = [10, 20, 40, 60, 80, 100]

data = []
for n in n_list:
    data.append(read_data(n))
    
ratio = data / data[-1]

if (args.reac != "iso"):

#########################################################################

    fig_name = "Emissivity as a function of quadrature points"

    fig, axs = plt.subplots(2, 3, sharex="all", gridspec_kw={'height_ratios': [3, 1]}, num=fig_name)

    # Suptitle
    plt.suptitle(fig_name)

    # Ylabels
    axs[0][0].set_ylabel(r"$j~[{\rm s}^{-1}]$")
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

    for id_n,n in enumerate(n_list):
        data_plot = data[id_n]
    
        for idx in range(3):
            axs[0][idx].plot(data_plot[0], data_plot[idx+1], marker="", label=r"$n=%d$" %n)
            axs[1][idx].plot(data_plot[0], ratio[id_n][idx+1], marker="", label=r"$n=%d$" %n)

    axs[0][0].legend()
    axs[1][0].legend()

#########################################################################

    fig_name = "Absoprtion opacity as a function of quadrature points"
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
    axs[1][0].set_xlabel("Neutrino energy [MeV]")

    for id_n,n in enumerate(n_list):
        data_plot = data[id_n]
    
        for idx in range(3):
            axs[0][idx].plot(data_plot[0], data_plot[idx+5], marker="", label=r"$n=%d$" %n)
            axs[1][idx].plot(data_plot[0], ratio[id_n][idx+5], marker="", label=r"$n=%d$" %n)

    axs[0][0].legend()
    axs[1][0].legend()

#########################################################################

if ((args.reac == "iso") | (args.reac == "all")):

    fig_name = "Scattering opacity as a function of quadrature points"
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
    axs[1][0].set_xlabel("Neutrino energy [MeV]")

    for id_n,n in enumerate(n_list):
        data_plot = data[id_n]
    
        for idx in range(3):
            axs[0][idx].plot(data_plot[0], data_plot[idx+13], marker="", label=r"$n=%d$" %n)
            axs[1][idx].plot(data_plot[0], ratio[id_n][idx+13], marker="", label=r"$n=%d$" %n)

    axs[0][0].legend()
    axs[1][0].legend()

plt.show()