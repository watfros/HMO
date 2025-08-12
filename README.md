# HMO
This repository contains programs for calculating the Hückel Molecular Orbital (HMO) theory using both Fortran and Python. 
HMO theory is a fundamental method in computational chemistry for determining the molecular orbital energies and wavefunctions of conjugated hydrocarbons.

## Fortran Implementation
The Fortran program is designed to efficiently compute the molecular orbital energies and eigenvectors using the Hückel approximation. Fortran is well-suited for scientific computing due to its high performance and robust numerical capabilities 
### Usage:
The HMO-test.exe is already compiled and suitable for Windows systems.
Input file format can be referenced from 1.dat, which defines the molecular connectivity.
Prepare your input file (e.g., input.dat) with the same format as 1.dat.
Place the input file in the same directory as HMO-test.exe.
Double-click HMO-test.exe to run the program, which will automatically read the input and output the results.
Note: The input file must be correctly formatted; otherwise, the program may crash or produce incorrect results.

## Python Implementation
The Python script provides a more modern and accessible approach to solving the same problem. It leverages Python's simplicity and powerful libraries like NumPy for linear algebra operations, making it easier to understand and modify for educational purposes.
### Usage:
The HMO.py script can be run directly or used with readxyz.py, which generates input data from an .xyz file that can then be copied into HMO.py.
