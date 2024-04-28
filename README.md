# `mutil` MATLAB Utilities for Trajectory Optimization

See example implementations provided at [`https://github.com/purnanandelango/ct-scvx`](https://github.com/purnanandelango/ct-scvx)

## `scp` Successive Convexification with Continuous-Time Constraint Satisfaction

 - `scp.ctscvx_[...]_handparse_noparam`  

## `disc` Discretization and Parameterization

 - First-order hold (FOH)
 - Zero-order hold (ZOH)
 - Finite-burn pulse (FBP)
 - Impulse

## `plant` System Model

- Double integrator with drag
- 6-DoF rocket with drag

## `solvers` Conic Optimization Solver

 - `pipg` Extrapolated Proportional Integral Projected Gradient Method

### Requirements for *node-only* constrained trajectory optimization

 - [YALMIP](https://yalmip.github.io/)
 - [OSQP](https://osqp.org/) &ndash; for solving quadratic programs
 - [ECOS](https://github.com/embotech/ecos) &ndash; for solving second-order cone programs
