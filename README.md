<!-- Author: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) -->

# ConformalEulerDeepWater: the full Euler equations of water waves in deep water, in conformal variables

![Focusing of a Peregrine breather computed by ConformalEulerDeepWater](docs/figures/conformal-breather.png)

The first 26 seconds of the run of `run_Euler.m`, computed in GNU Octave with the parameters of the script (deep
water, carrier period 0.75 s, wavelength 0.878 m, amplitude 0.01 m, 16384 Fourier modes). Panel (a) shows the free
surface at t = 0, a wave train whose envelope is the Peregrine breather of the nonlinear Schrödinger equation, and at
t = 25.7 s, when the modulation of the train, carried along at the group velocity, has grown into a single steep
crest of 0.076 m at x = 16.7 m, 4.1 times the highest crest at t = 0 (k₀η ≈ 0.54). Panel (b) shows the highest crest
of the free surface every 0.1 s, and panel (c) the time step chosen by the adaptive scheme, which falls from 10⁻² s to
below 10⁻⁴ s as this crest steepens.

ConformalEulerDeepWater is a set of MATLAB scripts, which also run in GNU Octave, that solve the full Euler equations
of irrotational surface gravity waves on water of infinite depth (deep water). The time-dependent fluid domain is
mapped conformally onto a fixed half-plane, a technique pioneered by L. V. Ovsyannikov and developed by Dyachenko,
Kuznetsov, Spector and Zakharov (1996), and the equations are solved by a Fourier pseudo-spectral method with an
exponential filter, exact integration of the linear terms and an adaptive embedded Runge–Kutta scheme. The example is
the evolution of the Peregrine breather, the prototype of extreme (rogue) waves, initialised as by Shemer and
Alperovich (2013). No article accompanies the code: the software itself is cited (see How to cite).

## The problem

The fluid is ideal and incompressible, the flow irrotational and two-dimensional, the depth infinite, and gravity
g = 9.81 m/s² acts alone; all quantities are in SI units. The fluid domain is mapped conformally onto the lower
half-plane of a variable ξ + iζ, so that the free surface becomes x(ξ, t) + iη(ξ, t). The unknowns are the elevation
η(ξ, t) and the velocity potential φ(ξ, t) at the surface, which obey

```text
η_t = η_ξ H[ψ_ξ / J] − x_ξ ψ_ξ / J,
φ_t = φ_ξ H[ψ_ξ / J] − (φ_ξ² − ψ_ξ²) / (2J) − g η,
x = ξ − H[η],     ψ_ξ = H[φ_ξ],     J = x_ξ² + η_ξ²,
```

where ψ is the stream function at the surface and H is the Fourier multiplier of symbol i sign(k) in ξ (a Hilbert
transform, with this sign). The linear part of these equations, whose waves obey ω² = g|k|, is integrated exactly in
Fourier space (an integrating factor), and the nonlinear part is computed on the grid.

The functions are expanded on N = 16384 Fourier modes of the periodic interval −l < ξ ⩽ l, with l = 35 λ₀ = 30.74 m.
After every step the spectrum is multiplied by the exponential filter exp(−5 (|k|/k_c)¹⁵), with k_c = 0.6 times the
largest wavenumber, a filter of the kind proposed by Hou and Li (2007). Time is advanced by the embedded Runge–Kutta
scheme of order 5(4) of Cash and Karp (1990), whose step is chosen by a PI control (Hairer and Wanner, 1996) with the
tolerance 5 × 10⁻¹⁰.

The initial wave train has the carrier period T₀ = 0.75 s, hence the frequency ω₀ = 8.378 rad/s, the wavenumber
k₀ = ω₀²/g = 7.154 m⁻¹, the wavelength λ₀ = 0.878 m and the group velocity c_g = 0.585 m/s, and the amplitude
0.01 m, that is a steepness ε = 0.051. Its envelope is the Peregrine solution of the nonlinear Schrödinger equation,
written as by Shemer and Alperovich (2013) with the shifts x0 = −17.74 m and Tmax = −30.30 of the script, and it is
cut off by a window that vanishes for |ξ| > l/2. The elevation is transferred to the conformal variable by a
fixed-point iteration (tolerance 10⁻¹⁴), and the potential is that of linear waves, φ̂ = −i sign(k) √(g/|k|) η̂.

## Quick start

From the folder `scripts/`, in MATLAB or GNU Octave:

```matlab
cd scripts
run_Euler
```

The script takes no argument; its parameters are set at its top. It draws the free surface and the Fourier spectrum
of its elevation every 0.1 s, and keeps the mass, the momentum and the energy at the same times in the variables
`Time`, `Mass`, `Moment` and `Energy`. Setting `compute_velocities = true` computes, at the same times, the velocities
`U` and `V` at the surface and their horizontal derivatives `U_x` and `V_x`, whose formulas are derived in
`doc/VelocityField.pdf`. The final time set by the script is Tf = 223.2 s (it writes 95 T₀ √g). In GNU Octave 11.1,
with the drawing replaced by the recording of the figure scripts, the computation reached t = 23 s in about 75 s; the
time step then falls below 10⁻⁴ s as the crest of the figure steepens, and after five minutes the run had reached only
t = 27.8 s.

## Contents

| Path | Content |
|---|---|
| `scripts/run_Euler.m` | the main script: parameters, initial breather, time loop, mass, momentum, energy and the optional surface velocities |
| `scripts/RHS.m` | the nonlinear right-hand side of the equations in conformal variables |
| `scripts/turn.m` | the exact integration of the linear terms, with the exponential filter |
| `scripts/Conf2Real.m` | the physical abscissa x = ξ − H[η] of the surface |
| `scripts/Plot.m` | the drawing of the free surface and of its spectrum |
| `doc/VelocityField.tex`, `doc/VelocityField.pdf` | the derivation of the surface velocities and of their gradients in conformal variables, by Francesco Fedele and Denys Dutykh |
| `pics/PB_shot.png` | the drawing of `Plot.m` at the start of the run, shown by the first version of this page |
| `docs/figures/` | the scripts of the figure of this page |

The figure is made by two scripts, run from `docs/figures/` in this order. The first, `conformal_readme_data.m`
(GNU Octave or MATLAB, about two minutes), runs `run_Euler.m` as it stands from the folder `docs/figures/recording/`,
whose `Plot.m` records the free surface and the time step instead of drawing them and stops the run at t = 26 s. The
second, `plot_readme_figure.py` (Python 3 with NumPy and Matplotlib, and the style file beside it), draws the figure.

## Requirements

MATLAB or GNU Octave (tested with Octave 11.1); no toolbox is needed.

## How to cite

No article accompanies this code, so please cite the software itself:

D. Dutykh, ConformalEulerDeepWater: the full Euler equations of water waves in deep water, in conformal variables, by
a Fourier pseudo-spectral method, version 1.1.0 (2026).

https://github.com/dutykh/ConformalEulerDeepWater

```bibtex
@software{Dutykh_ConformalEulerDeepWater,
  author = {Dutykh, Denys},
  title = {{ConformalEulerDeepWater}: the full {Euler} equations of water waves in deep water, in conformal variables, by a {Fourier} pseudo-spectral method},
  year = {2026},
  version = {1.1.0},
  url = {https://github.com/dutykh/ConformalEulerDeepWater},
  license = {MIT}
}
```

The button "Cite this repository" in the right-hand column of this page gives the same reference, read from
`CITATION.cff`. The code rests on the following works, which deserve credit wherever it is used:

| Work | Used for |
|---|---|
| A. I. Dyachenko, E. A. Kuznetsov, M. D. Spector and V. E. Zakharov, Analytical description of the free surface dynamics of an ideal fluid (canonical formalism and conformal mapping), Physics Letters A 221, 73–79 (1996), [doi:10.1016/0375-9601(96)00417-3](https://doi.org/10.1016/0375-9601(96)00417-3) | the equations in conformal variables |
| D. H. Peregrine, Water waves, nonlinear Schrödinger equations and their solutions, The Journal of the Australian Mathematical Society, Series B 25, 16–43 (1983), [doi:10.1017/S0334270000003891](https://doi.org/10.1017/S0334270000003891) | the Peregrine breather |
| L. Shemer and L. Alperovich, Peregrine breather revisited, Physics of Fluids 25, 051701 (2013), [doi:10.1063/1.4807055](https://doi.org/10.1063/1.4807055) | the initial breather |
| T. Y. Hou and R. Li, Computing nearly singular solutions using pseudo-spectral methods, Journal of Computational Physics 226, 379–397 (2007), [doi:10.1016/j.jcp.2007.04.014](https://doi.org/10.1016/j.jcp.2007.04.014) | the exponential filter |
| J. R. Cash and A. H. Karp, A variable order Runge–Kutta method for initial value problems with rapidly varying right-hand sides, ACM Transactions on Mathematical Software 16, 201–222 (1990), [doi:10.1145/79505.79507](https://doi.org/10.1145/79505.79507) | the embedded Runge–Kutta scheme |
| E. Hairer and G. Wanner, Solving Ordinary Differential Equations II. Stiff and Differential-Algebraic Problems, Springer (1996), [doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7) | the PI control of the time step |

## Version history

Version 1.0, the tag `v1.0` of May 2025, holds the solver as it stands, with the computation of the surface velocities
contributed by Francesco Fedele. Version 1.1.0 (November 2026) leaves the solver unchanged and documents it: the
equations and the parameters of the breather, the figure of its focusing with the scripts that make it, the citation
file and the credits.

## Licence

MIT licence (see `LICENSE`).

## Author and acknowledgements

Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi).

https://www.denys-dutykh.com/

Comments, questions and reports of errors are welcome. The author thanks the colleagues who helped to develop the
solver: Didier Clamond (Université Côte d'Azur, Laboratoire J.A. Dieudonné, Nice), for the conformal method and its
implementation; Bernard Ee (Tel Aviv University, School of Mechanical Engineering), for the initialisation of the
Peregrine breather and its validation; Francesco Fedele (Georgia Institute of Technology, Atlanta), for the module
that computes the surface velocities and their derivatives, and for the first draft of `doc/VelocityField.tex`; and
Lev Shemer (Tel Aviv University, School of Mechanical Engineering), for the theory of the Peregrine breather and its
experimental validation.
