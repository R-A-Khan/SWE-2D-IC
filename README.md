# SWE-2D-IC
Extension to 2D of data assimilation scheme used to reconstruct missing initial data for shallow water equations using existing measurements, with applications to tsunami simulation. Full description of model derivation and evaluation can be found at:

R. A. Khan & N. K.-R. Kevlahan (2022) Data assimilation for the two-dimensional shallow water equations: Optimal initial conditions for tsunami modelling, Ocean Modelling, Vol 174, DOI: //doi.org/10.1016/j.ocemod.2022.102009

## Matlab Usage
Adjust model parameters in 
* __run_data_assimil_2D.m__  
 
Utilises the data assimilation algorithm in __data_assimil_2D.m__.

Verification of numerical solvers used in DA scheme done using a kappa convergence test in 
* __Kappa_test_2D.m__

Visualize results using 
* __Plots.m__ 

Built using Matlab 2018a, including Optimization toolbox.

## Fortran Usage
Adjust model parameters in 
* __data_assimil_2D.f95__  
* Run command `gfortran -o test mod_procedures.f95 mod_SWEandAdjointSolvers.f95 mod_RungeKuttaSolvers.f95 mod_FWandBWSolvers.f95 mod_ParametersDA.f95  data_assimil_2D.f95 
./test`
 

Verification of numerical solvers used in DA scheme done using a kappa convergence test in 
* __kappa_test_2D.f95__


Built using GNU Fortran.


