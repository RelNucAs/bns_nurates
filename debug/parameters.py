## Parameters for Python debugging routine

# Weak magnetism
#WM_abs   = False #for (anti)neutrino absorption on nucleons
#WM_scatt = False #for (anti)neutrino scattering on nucleons
WM = False

# Fortran data file
if WM:
    f_file = '../input/nurates_1.008E+01_WM.txt'
else:
    f_file = '../input/nurates_1.008E+01.txt'


# C++ data file
if WM:
    C_file = '../input/newrates_WM.txt'
else:
    C_file = '../input/newrates.txt'
