# Hartree Fock Method
Computes the Hartree Fock energy for CUO2 using model parameters.

## Table of Contents:

1. [Basics](README.md#basics)

1. [System Description](README.md#system-descprition) 

1. [Versions](README.md#versions) 

### Git Directory at origin:

`/Users/pcs/Documents/LSCO` --> Note: Do not change this

### Working Directory:

Current: `/home/niladri/HF_Ara/TEST/VERSION/JH0.0` 


## Basics
Running this code need ARMADILLO. Make sure to install ARMADILLO package in your
system.  

#### To compile:

First load `module.sh`

`g++ -O2 -g -Wall -o hf main.cpp lattice.h matrix.cpp -fopenmp  -I/opt/blas+lapack/3.8.0/gcc-4.8.5/include -L/opt/blas+lapack/3.8.0/gcc-4.8.5/lib -llapack -lblas -I/opt/armadillo/9.200.7/gcc-4.8.5/include -L/opt/armadillo/9.200.7/gcc-4.8.5/lib64 -larmadillo -I/opt/intel/2017up1/compilers_and_libraries_2017.1.132/linux/mkl/include`

#### Run:

example of sample run:

`./hf <U> <H>`

where `U` is the electron correlation and `H` is external magnetic field. 


#### Plotting:

Python scripts for plotting are added.
      
       




### v1.00

OPENMP version added. 




