#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the Peregrine breather of Shemer & Alperovich (2013) computed by run_Euler.m.

Run conformal_readme_data.m first (GNU Octave or MATLAB), then, from this folder,
    python3 plot_readme_figure.py
which writes conformal-breather.png. The style is the stamped copy manuscript_style.py beside it.

Panel (a): the free surface around the centre of the group, x - cg*t, when the centre leaves the
wavemaker and when it passes the gauges at 2.25 m, 6.25 m and 9 m, shifted by 40 mm each, as in
Figure 1 of the article. Panel (b), the counterpart of Figure 3 of the article: the highest
elevation met at each distance from the wavemaker during the run, the largest envelope of the
Peregrine solution, and the highest crests measured at the three gauges (the values given in
the text of the article). Panel (c): the changes of the energy and of the momentum since the
first output time, relative to their values then, and the change of the mass divided by zeta0
times the length of the group (35 wavelengths).
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

# Parameters of run_Euler.m (SI units)
G, T0, ZETA0, X0 = 9.81, 0.587, 0.01, 9.0
OM0 = 2 * np.pi / T0
K0 = OM0**2 / G
L0 = 2 * np.pi / K0
CG = 0.5 * np.sqrt(G / K0)
A0 = ZETA0 / np.sqrt(2)
EPS = A0 * K0
HALF = 35 * L0                                  # half-length of the periodic domain
GAUGES = np.array([2.25, 6.25, 9.0])            # wave gauges of Figure 1 of the article (m)
MEASURED = np.array([14.2, 20.1, 25.2])         # highest crests measured there (mm)

rec = np.loadtxt(HERE / "record.csv", delimiter=",")
x_u, crest_u = np.loadtxt(HERE / "maxima.csv", delimiter=",", unpack=True)
snap = np.loadtxt(HERE / "snapshots.csv", delimiter=",")
t_inv, mass, moment, energy = np.loadtxt(HERE / "invariants.csv", delimiter=",", unpack=True)
t = rec[:, 0]

computed = 1000 * np.interp(GAUGES, x_u, crest_u)   # highest crest at the gauges (mm)


def peregrine(x):
    """Largest envelope of the Peregrine solution when its centre is at the distance x."""
    eta = EPS**2 * K0 * (x - X0)
    return np.sqrt(2) * A0 * np.abs(1 - 4 * (1 - 4j * eta) / (1 + 16 * eta**2))


hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig = plt.figure(figsize=(hs.WIDTH, 4.9))
grid = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0])
top = fig.add_subplot(grid[0, :])
left = fig.add_subplot(grid[1, 0])
right = fig.add_subplot(grid[1, 1])

# (a) free surface in the frame of the group
instants = np.unique(snap[:, 0])
names = ["leaving the wavemaker", "at 2.25 m", "at 6.25 m", "at 9 m"]
side = blended_transform_factory(top.transAxes, top.transData)
for n, (when, name) in enumerate(zip(instants, names)):
    rows = snap[:, 0] == when
    s = np.mod(snap[rows, 1] - CG * when + HALF, 2 * HALF) - HALF
    order = np.argsort(s)
    top.plot(s[order], 1000 * snap[rows, 2][order] + 40 * n, color=hs.NAVY, lw=0.6)
    top.text(1.01, 40 * n, name, transform=side, ha="left", va="center",
             fontsize=plt.rcParams["legend.fontsize"], color=hs.INK)
top.text(1.01, 168, "centre of the group", transform=side, ha="left", va="center",
         fontsize=plt.rcParams["legend.fontsize"], color=hs.SLATE)
top.set_xlim(-8, 8)
top.set_ylim(-30, 175)
top.set_xlabel(r"$x - c_g t$ (m)")
top.set_ylabel(r"$\eta$ (mm), shifted by 40 mm")
hs.panel(top, "(a) free surface around the centre of the group")

# (b) maxima along the tank
reach = np.linspace(0, X0, 400)
seen = (x_u >= 0) & (x_u <= X0)
left.plot(x_u[seen], 1000 * crest_u[seen], color=hs.NAVY, lw=0.7, label="highest crest, computed")
left.plot(reach, 1000 * peregrine(reach), color=hs.SLATE, lw=0.9, linestyle=(0, (1.6, 1.6)),
          label="envelope, Peregrine")
left.plot(GAUGES, MEASURED, linestyle="none", marker="D", markersize=4.0, markerfacecolor="white",
          markeredgecolor=hs.WINE, label="highest crest, measured")
left.set_xlim(0, X0 + 0.25)
left.set_ylim(10, 32)
left.set_xlabel("distance from the wavemaker (m)")
left.set_ylabel("elevation (mm)")
hs.panel(left, "(b) growth along the tank")
hs.legend_outside(left, where="above", ncol=2)

# (c) conservation
right.semilogy(t_inv, np.abs(energy - energy[0]) / energy[0], color=hs.NAVY, lw=1.0,
               label="energy")
right.semilogy(t_inv, np.abs(moment - moment[0]) / abs(moment[0]), color=hs.WINE, lw=1.0,
               linestyle=(0, (5.0, 2.0)), label="momentum")
right.semilogy(t_inv, np.abs(mass - mass[0]) / (ZETA0 * HALF), color=hs.TEAL, lw=1.0,
               linestyle=(0, (1.6, 1.6)), label="mass")
right.set_xlim(0, t[-1])
right.set_ylim(1e-11, 1e-4)
right.set_xlabel(r"$t$ (s)")
right.set_ylabel("relative change")
hs.grid(right)
hs.panel(right, "(c) invariants")
hs.legend_outside(right, where="above", ncol=3)

out = HERE / "conformal-breather.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out)
print("highest crest at the gauges", GAUGES, "m: computed", np.round(computed, 2),
      "mm, measured", MEASURED, "mm; envelope of the Peregrine solution",
      np.round(1000 * peregrine(GAUGES), 2), "mm")
print(f"highest elevation {1000 * rec[:, 1].max():.2f} mm at t = {t[np.argmax(rec[:, 1])]:.3f} s,"
      f" x = {rec[np.argmax(rec[:, 1]), 2]:.2f} m")
print(f"largest changes: energy {np.abs(energy - energy[0]).max() / energy[0]:.2e},"
      f" momentum {np.abs(moment - moment[0]).max() / abs(moment[0]):.2e},"
      f" mass {np.abs(mass - mass[0]).max():.2e} m^2")
