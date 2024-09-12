# mixed-basis-be
mixed-basis calculations based on QuEmb repository from oimeitei.

Allows for BE calculations to perfomed from a mixed-basis reference. 
Rather than using a single large-basis HF reference wavefunction, 
mixed-basis-be generates a set of BE(n)-in-BE(m) HF reference calculations
with a dense-in-minimal basis set for a given BE(m) region. 
For sufficiently large systems, this reduces the time and memory costs for 
both the HF reference steps and the integral transformations to each integral 
space while preserving BE(n) accuracy. This allows for dense basis, large accuracy
BE calculations to be performed for extended chemical systems.

TODO: add documentation and guide to run the code.

This work, nomenclature, and related results are given in the paper on arXiv: 
Multiscale Embedding for Quantum Computing, arXiv:2409.06813
