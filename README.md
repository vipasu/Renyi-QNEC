# Renyi Quantum Null Energy Condition

The Quantum Null Energy condition is a statement about lower bounds on energy density in quantum field theory. Its proof in general was recently shown using techniques from quantum information theory by studying the relative entropy. We considered a Renyi generalization called the Sandwiched Renyi Divergence (SRD) and prove its validity for free field theories in dimensions greater than 2 for n > 1. The arXiv preprint can be found [here](https://arxiv.org/abs/2007.15025).

This repository contains code for the numerical study of SRD in critical spin chains that are known to be dual to 1+1 Conformal Field Theories in the low energy sector.

We make use of the [ITensors.jl](https://github.com/ITensor/ITensors.jl) library for Density Matrix Renormalization Group (DMRG) methods for numerical access to these critical spin chain states.

Relevant figures are presented in the paper and all the code used to generate them are contained in this repository

## Features
- Different hamiltonians (TFIM, XXZ, WZW SU(2)_2)
- Iteratively calculated excited ground states
- Entanglement entropies of chains, relative entropies and sandwiched renyi divergences
- TODO Check the support condition for relative entropies
