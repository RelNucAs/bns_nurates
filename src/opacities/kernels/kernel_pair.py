import numpy as np
import sympy as sp
import sys

print("Script used to generate the expressions for the optimized pair process"
      "kernel functions, from:")
print("  Pons et al., Legendre expansion of the nu nubar <--> eplus eminus kernel: Influence of high order terms, Astron. Astrophys. Suppl. Ser. 129, 343-351 (1998)")
print("Run it twice, once with argument y (meaning neutrino energy < antineutrino energy),")
print("once with argument z (meaning neutrino energy > antineutrino energy)")
print("")


# Sympy variables
eta = sp.symbols("eta", real=True) # Electron degeneracy
y = sp.symbols("y", real=True, positive=True) # Neutrino energy rescaled with temperature
z = sp.symbols("z", real=True, positive=True) # Antieutrino energy rescaled with temperature

# A generic (from sympy point of view) function, representing a complete Fermi-Dirac integral
F = sp.Function("F")

# Get from the user which energy is the smaller/greater
if sys.argv[1] == "y":
    m = y
    M = z
elif sys.argv[1] == "z":
    m = z
    M = y
else:
    exit("Pass either y or z")

# A dictionary useful to manipulate the way the Fermi-Dirac integrals are handled and printed to screen
FDI_names = {}
for i in range(6):
    for arg in (eta, eta - y, eta - z, eta + y, eta + z, eta - y - z, eta + y + z):
        name = f"FDI_{'p' if i > 0 else '':s}{i:d}_e"
        if arg == eta - y:
            name += "my"
        elif arg == eta - z:
            name += "mz"
        elif arg == eta + y:
            name += "py"
        elif arg == eta + z:
            name += "pz"
        elif arg == eta - y - z:
            name += "myz"
        elif arg == eta + y + z:
            name += "pyz"
        FDI_names[F(i, arg)] = sp.symbols(name, real=True, positive=True)


# Incomplete Fermi-Dirac integral, with integration endpoints called u and v.
# This function expresses the incomplete Fermi integrals in terms of the
# complete ones. Note: only makes sense when n is an integer!
def F_inc(n, eta, u, v):

    k = sp.symbols("k", integer=True)

    if v == +sp.oo:
        expr = u**(n - k) * F(k, eta - u)
    else:
        expr = u**(n - k) * F(k, eta - u) - v**(n - k) * F(k, eta - v)

    return sp.Sum(sp.binomial(n, k) * expr, (k, 0, n)).doit().subs(FDI_names)


# G function defined in eq. 12 of Pons et al.
def G(n, eta, u, v):

    return F_inc(n, eta, u, v) - F_inc(n, eta + y + z, u, v)


# Coefficients a, c and d, only for l=1, defined in appendix A of Pons et al.
def a(n, u, v):
    if n == 3:
        return 8/(3*u**2)
    elif n == 4:
        return -4/(3*u**2*v)
    elif n == 5:
        return 4/(15*u**2*v**2)


def c(n, u, v):
    if n == 0:
        return 4*u/v**2*(2*v**2/3+u*v+2*u**2/5)
    elif n == 1:
        return -4*u/(3*v**2)*(3*u+4*v)
    elif n == 2:
        return 8*u/(3*v**2)


def d(n, u, v):
    if n == 0:
        return 4*v**3/(15*u**2)
    elif n == 1:
        return -4*v**2/(3*u**2)
    elif n == 2:
        return 8*v/(3*u**2)


# Function Psi defined in eq. 11 of Pons et al. (only for l=0)
def Psi(eta, u, v):

    r = 0
    for n in range(3):
        r += c(n, u, v) * G(n, eta, u, u + v) + d(n, u, v) * G(n, eta, v, u + v)
    for n in range(3, 6):
        r += a(n, u, v) * (G(n, eta, 0, m) - G(n, eta, M, u + v))

    return r


# for n in range(6):
    # print(f"G_{n:d}(y, y + z):  ", sp.horner(sp.simplify(sp.expand(G(n, eta, y, y + z))), wrt=y))
    # print(f"G_{n:d}(z, y + z):  ", sp.horner(sp.simplify(sp.expand(G(n, eta, z, y + z))), wrt=z))
    # print(f"G_{n:d}(0, min(y, z)):  ", sp.horner(sp.simplify(sp.expand(G(n, eta, 0, m))), wrt=m))
    # print("")

psiyz = sp.horner(sp.simplify(sp.expand(Psi(eta, y, z) * (15*y**2*z**2))), wrt=y)
psizy = sp.horner(sp.simplify(sp.expand(Psi(eta, z, y) * (15*y**2*z**2))), wrt=z)

print("Psi(y, z) * (15*y**2*z**2):\n", psiyz)
print("")
print("Psi(z, y) * (15*y**2*z**2):\n", psizy)
print("")

print("Fermi integrals actually appearing in the expressions for Psi(y,z) and Psi(z,y):  ")
i = 0
for fn, f in FDI_names.items():

    if len(psiyz.find(f) | psizy.find(f)) == 0:
        continue
    else:
        fn = str(fn).replace("F", "FDI")
        for n in range(6):
            fn = fn.replace(f"({n:d}, ", f"_p{n:d}(")
        print(i, "  ", f, "=", fn)
        i += 1
print("")

replacements, reduced_exprs = sp.cse([psiyz, psizy])

print("Shorthands suggested by common subexpression elimination algorithm:")
for r in replacements:
    print("   ", r[0], "=", r[1])
print("")


print("Psi(y, z) * (15*y**2*z**2) (but now as a function of the shorthands):\n", reduced_exprs[0])
print("")
print("Psi(z, y) * (15*y**2*z**2) (but now as a function of the shorthands):\n", reduced_exprs[1])
print("")
