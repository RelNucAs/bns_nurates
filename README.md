# NuO-BNSlib
A C++ library of Neutrino opacities in the context of BNS mergers

The Opacities are based on Burrows(2004), Bruenn (1985) and references there in. 

We implemented the following reactions's opacities:

1. Neutrino absorption on protons
Eq:C13 from Bruenn(1985) with blocking factor and Weak magnetism from Burrows(2004) included.

2. Neutrino absorption on neutrons
Eq:C19 from Bruenn(1985) with blocking factor and Weak magnetism from Burrows(2004) included.

Code description:

main.cpp contains all the key routines for calculating the opacities. 

To create an exe:

g++ -o nu_rates main.cpp -I. 

You can also do make in the main folder that contains the makefile which creates a exe with the name nu_rates.
