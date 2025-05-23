Full Euler Solver on Deep Water
=======================

![Peregrine breather](pics/PB_shot.png)

## Overview

The present collection of [MATLAB](http://mathworks.com/products/matlab/) scripts provides a high-fidelity Fourier pseudo-spectral solver for the full Euler equations with free surface on a fluid layer of infinite depth (the so-called *deep water approximation*). The time-dependent fluid domain is transformed into a fixed computational strip using the conformal mapping technique (pioneered by L.V. Ovsyannikov in 1972 and further developed by A. Dyachenko *et al.* in 1996).

### Key Features

* **Conformal Mapping**: Transforms the time-dependent physical domain into a fixed computational strip
* **Pseudo-spectral Method**: High-order spatial discretization using FFT with exponential dealiasing
* **Adaptive Time-stepping**: Embedded Cash-Karp 5(4) method with PI control for optimal time step selection
* **Integrating Factor Technique**: Exact integration of linear dispersive terms
* **Surface Velocity Computation**: Optional post-processing to compute horizontal and vertical velocities (NEW)
* **Conservation Monitoring**: Real-time tracking of mass, momentum, and total energy

## Mathematical Framework

The solver implements the following key mathematical components:

1. **Spatial Discretization**: Fourier pseudo-spectral method with N = 16384 modes
2. **Time Integration**: Cash-Karp embedded Runge-Kutta method of order 5(4)
3. **Dealiasing**: Exponential filter following Hou et al. with Λ = 15
4. **Linear Terms**: Exact integration using rotation matrices
5. **Conformal Transform**: x = ξ - H[η], where H is the Hilbert transform

## Running the Simulation

To execute the simulation, simply run:
```matlab
run_Euler
```

No command-line arguments are required as all parameters are defined within the script.

### Configuration Options

The main script `run_Euler.m` includes a toggle for velocity computation:
```matlab
compute_velocities = false;  % Set to true to compute surface velocities
```

When enabled, the solver computes:
* **U**: Horizontal surface velocity
* **V**: Vertical surface velocity  
* **U_x**: Horizontal derivative of U
* **V_x**: Horizontal derivative of V

## Code Structure

* **`run_Euler.m`**: Main simulation driver with adaptive time-stepping
* **`RHS.m`**: Computes nonlinear right-hand side for conformal Euler evolution
* **`turn.m`**: Handles linear rotation steps with spectral filtering
* **`Conf2Real.m`**: Transforms from conformal space to physical coordinates
* **`Plot.m`**: Real-time visualization of surface elevation and Fourier spectrum

## Documentation

* **`doc/VelocityField.tex`**: Mathematical derivation of surface velocities and gradients in conformal mapping for water waves. This document provides:
  - Derivation of velocity components from the complex potential
  - Transformation formulas between physical and conformal coordinates
  - Laplace equation in conformal coordinates and the relationship between velocity potential and streamfunction
  - Explicit formulas for computing surface velocity gradients
  - Authors: Francesco Fedele (Georgia Tech) and Denys Dutykh (Khalifa University)

## Physical Problem

The solver is initialized to simulate the celebrated Peregrine breather evolution in the full Euler equations. This is a fundamental nonlinear wave solution that models extreme wave events. We refer to the following publication for more details:

* L. Shemer & L. Alperovich. *Peregrine breather revisited*. Phys. Fluids, **25**, 051701, 2013

## Author

**Denys Dutykh**  
Khalifa University of Science and Technology  
Abu Dhabi, UAE  
GitHub: [https://github.com/dutykh/](https://github.com/dutykh/)

Please don't hesitate to contact the author with any questions, remarks, and bug reports. The latest contact details can be found at:

* [http://www.denys-dutykh.com](http://www.denys-dutykh.com/)

## Acknowledgements

The author would like to thank the following colleagues (in alphabetical order) who helped develop and enhance this solver:

* **[Didier Clamond](http://math.unice.fr/~didierc/)**, [University of Nice Sophia Antipolis](http://unice.fr/), [Laboratoire J.A. Dieudonné](http://math.unice.fr/), France  
  *Contributions: Conformal mapping methodology and numerical implementation*

* **[Bernard Ee](https://www.researchgate.net/profile/Bernard_Ee2)**, [Tel Aviv University](http://english.tau.ac.il/), [School of Mechanical Engineering](http://engineering.tau.ac.il/), Israel  
  *Contributions: Peregrine breather initialization and validation*

* **[Francesco Fedele](https://ce.gatech.edu/directory/person/francesco-fedele)**, Associate Professor, School of Civil & Environmental Engineering and School of Electrical & Computer Engineering, [Georgia Institute of Technology](https://gatech.edu/), Atlanta, GA, USA  
  *Contributions: Surface velocity computation module and physical derivatives transformation. Special thanks for initiating and preparing the initial draft of the mathematical documentation (VelocityField.tex)*

* **[Lev Shemer](http://www.eng.tau.ac.il/~shemer/)**, [Tel Aviv University](http://english.tau.ac.il/), [School of Mechanical Engineering](http://engineering.tau.ac.il/), Israel  
  *Contributions: Peregrine breather theory and experimental validation*

## Version History

* **v1.1** (2025): Added surface velocity computation capability with contributions from Prof. Francesco Fedele
* **v1.0** (2025): Initial release with full Euler solver using conformal mapping technique

## License

This software is distributed under the MIT License. See the LICENSE file for details.