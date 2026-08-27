[Github readme.txt](https://github.com/user-attachments/files/31531163/Github.readme.txt)

Passive Modeling of Multiconductor Transmission Lines in MATLAB.

MATLAB implementation of a methodology for the passive modeling and time-domain simulation of frequency-dependent multiconductor transmission lines (MTLs).

The repository supports three-phase and double-circuit three-phase transmission lines, including configurations with ground wires and, a new case of MTL can be defined.

Main Features:

- Frequency-dependent transmission-line parameter calculation.
- Kron reduction and ideal transposition.
- Clarke modal transformation.
- Rational fitting using Vector Fitting (VF) and Vector Fitting with Real Poles (VF-RP).
- Passive residue calculation using Non-Negative Least Squares (NNLS).
- Passive Foster-type RL and cascaded π-equivalent circuit synthesis.
- Time-domain simulation using:
  1) State-space.
  2) Sparse state-space.
  3) EMTP-type.
- Numerical Laplace Transform (NLT).
- Optional passivity assessment.
- ATP input-file generation for the supported of the three-phase cases.

Software:

The project is implemented in MATLAB and uses the Vector Fitting (VF) software for rational approximation. ATP is required only for the optional ATP-based validation.

Quick Start:

Open Main.m and execute the sections sequentially:

1) Parameter calculation.
2) Rational fitting with VF, VF-RP, and NNLS.
3) Time-domain simulation.
4) ATP file generation and comparison (optional).
The main stages exchange data through the generated .mat files.

Case Studies:

The repository includes three-phase and double-circuit three-phase transmission-line configurations with optional ground wires. Custom geometries can also be entered through the graphical interface.

NNLS Solvers:

Several NNLS options are available, including MATLAB lsqnonneg, Lawson-Hanson, BW-NNLS, FNNLS, and TNT-NN. Some of these routines were developed by other authors and are included as computational components of the proposed methodology.

Citations:

1) If you use this code in academic work, please cite the associated paper describing the methodology and MATLAB implementation.
2) The Vector Fitting, NNLS algorithms, and other third-party methods and software used in this repository are properly referenced in the associated paper. Please also cite the original sources when required.

Authors:

-José de Jesús Reyes Ramírez.
-Eduardo Salvador Bañuelos Cabral.
-José Alberto Gutiérrez Robles.
-José de Jesús Nuño Ayón.

Acknowledgment:

The authors gratefully acknowledge the developers of the third-party algorithms and software used in this project.
