## Parameters for Python debugging routine

# Weak magnetism
#WM_abs   = False #for (anti)neutrino absorption on nucleons
#WM_scatt = False #for (anti)neutrino scattering on nucleons
use_dU = False

# Fortran data file
if use_dU:
    f_file = '../input/nurates_1.008E+01_dU.txt'
else:
    f_file = '../input/nurates_1.008E+01.txt'


# C++ data file
if use_dU:
    C_file = '../input/newrates_dU.txt'
else:
    C_file = '../input/newrates.txt'
