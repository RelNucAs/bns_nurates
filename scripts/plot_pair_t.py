import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

###############################################
# Plot PairT function from Pons et al. (1998) #
###############################################

## T_l(alpha) = \sum_{n=1}^\inf  (-1)^{n+1} * exp(-n*alpha) / n^l 

## Find index function
def find_index(array, value):
    return np.argmin(abs(array - value))

## Read data
filename = "../inputs/pair_t_tabulated_more_digits_log.txt"
data = np.loadtxt(filename, unpack=True)

alpha  = data[0]
pair_t = data[1:7] # T_k(alpha), k = 1,2,3,4,5,6

alpha_max = 700. # restrict alpha range since last numbers are just zeros
id_max = find_index(alpha, alpha_max)

alpha  = alpha[:id_max+1]
pair_t = pair_t[:,:id_max+1]

nt = pair_t.shape[0] # 6

def fit_function(alpha_fix, n):
    id_fix = find_index(alpha, alpha_fix)

    fit_array = []

    for idx in range(nt):
        #z1 = np.polyfit(alpha[id_fix+1:], np.log(pair_t[idx,id_fix+1:]), 1)
        #tmp_1 = np.exp(np.poly1d(z1)(alpha))
        
        tmp_1 = np.exp(-alpha) 

        ratio = pair_t[idx] / tmp_1

        #z2 = np.polyfit(alpha[:id_fix+1], ratio[:id_fix+1], n)
        z2 = np.polyfit(np.exp(-alpha[:id_fix+1]), ratio[:id_fix+1], n)
        print(z2)
        #z2 = np.flip([(-1)**(nn+1) / nn**(idx+1) for nn in range(1,n+2)])

        print(z2)
        print("")

        #tmp_2 = np.poly1d(z2)(alpha) * tmp_1
        tmp_2 = np.poly1d(z2)(np.exp(-alpha)) * tmp_1

        #print("{", end="")
        #for z in z2:
        #    print(str(z) +", ",end="")
        #print("},")

        fit = np.zeros(alpha.shape)

        fit[:id_fix+1] = tmp_2[:id_fix+1]
        fit[id_fix+1:] = tmp_1[id_fix+1:]
   
        fit_array.append(fit)

    return np.array(fit_array)

def residuals(fit):
    #return np.sqrt( np.sum( ((fit - pair_t) / pair_t)**2. ) )
    return np.sum( abs((fit - pair_t) / pair_t) )

alpha_fix = 25.
n = 12

fit = fit_function(alpha_fix, n)

print(residuals(fit))

#print("Residuals:")
#for n in range(5,21):
#    print("n = %d: %.5e" %(n, residuals(fit_function_pol(alpha_fix,n))))

#print("")

#n = 12

#for alpha_fix in range(5,21):
#     print("alpha_fix = %d: %.5e" %(alpha_fix, residuals(fit_function_pol(alpha_fix,n))))   

## Plot data
fig, axs = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.)

for idx in range(6):
    axs[0].plot(alpha, pair_t[idx], label=r"$T_{%d}$" %(idx+1))
    axs[0].plot(alpha, fit[idx], ls = "--", label="Fit")
    #axs[1].plot(alpha, (pair_t/fit)[idx], label=r"$T_{%d}$" %(idx+1))
    axs[1].plot(alpha, (abs(pair_t - fit)/pair_t)[idx], label=r"$T_{%d}$" %(idx+1))

axs[1].axhline(1., color="k", ls="--")
axs[1].axhline(1.0-1.0E-06, color="k", ls="--", alpha=0.5)
axs[1].axhline(1.0+1.0E-06, color="k", ls="--", alpha=0.5)

xlim_min = 0.
xlim_max = 700.

ylim_min = np.amin(pair_t[:,:find_index(alpha, xlim_max)+1])
ylim_max = np.amax(pair_t[:,:find_index(alpha, xlim_max)+1])

axs[0].set_xlim((xlim_min, xlim_max))
axs[0].set_ylim((ylim_min, ylim_max))

axs[0].set_xscale("log")
axs[0].set_xlim((1.0E-05, 7.0E+02))

axs[0].set_yscale("log")
axs[1].set_yscale("log")

axs[1].set_xlabel(r"$\alpha$")
axs[0].set_ylabel(r"$T(\alpha)$")
axs[1].set_ylabel(r"$T$-to-fit ratio")

axs[0].legend()
axs[1].legend()

plt.show()