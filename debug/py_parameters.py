## Parameters for Python debugging routine

# Weak magnetism
#WM_abs   = False #for (anti)neutrino absorption on nucleons
#WM_scatt = False #for (anti)neutrino scattering on nucleons

# dU correction
use_dU = True

# Energy of the comparison
E = "3.000e+02"

# Fortran data file
if use_dU:
    f_file = '../input/nurates_'+E.replace("e","E")+'_dU.txt'
else:
    f_file = '../input/nurates_'+E.replace("e","E")+'.txt'

# C++ data file
if use_dU:
    C_file = '../input/newrates_'+E+'_dU.txt'
else:
    C_file = '../input/newrates_'+E+'.txt'

# Plot folder
plotfolder = '../output/plots/'

# Max number of rows to be read
nrowmax = 104
