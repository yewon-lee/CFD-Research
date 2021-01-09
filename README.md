# CFD-Research
Implementation of Summation-by-Parts operators and Discontinuous Galerkin methods, which are numerical methods for solving PDEs (more specifically, aerodynamic flows). @ Computational Aerodynamics Lab, UTIAS.

Details:
- Higher-order Discontinuous Galerkin methods (energy stable scheme for varying orders p>=1)
- Summation-by-Parts operators (classical, p=2 and p=4 element-type with upwind and symmetric fluxes)

The above programs include functions to march the solution in time, compute energies at each time step, and plot the final solution. The final report summarizes detailed results and analyses from a comparative study between Summation-by-Parts and Discontinuous Galerkin methods.
