# mixed-basis-be
mixed-basis calculations based originally on QuEmb repository from oimeitei, 
and now later versions from quemb from troyvvgroup.

Allows for BE calculations to perfomed from a mixed-basis reference. 
Rather than using a single large-basis HF reference wavefunction, 
mixed-basis-be generates a set of BE(n)-in-BE(m) HF reference calculations
with a dense-in-minimal basis set for a given BE(m) region. 
For sufficiently large systems, this reduces the time and memory costs for 
both the HF reference steps and the integral transformations to each integral 
space while preserving BE(n) accuracy. This allows for dense basis, large accuracy
BE calculations to be performed for extended chemical systems.

To run the code, first install quemb using guide (https://github.com/troyvvgroup/quemb).
Then, simply provide the geometry and charge of the system interest; size of big basis
region, using the BE(x) prescription based on x connected coordination shells; chosen
basis for the big (BE(m)) and small basis set regions; and BE(n) choice for BE fragment
size.
We provide an example shell script that passes this information to the python script
interfacing with quemb. This script returns the correlation energy of the system, an
estimation of the mixed HF energy, nuclear and core energies, and number of electrons.
A crude chemical potential optimizer is added to optionally correct the number of
electrons in the mixed calculations.

TODO: further guide and documentation

This work, nomenclature, and related results are given in the paper on arXiv: 
Multiscale Embedding for Quantum Computing, arXiv:2409.06813
