#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the focusing of the Peregrine breather computed by run_Euler.m.

Run conformal_readme_data.m first (GNU Octave or MATLAB), then, from this folder,
    python3 plot_readme_figure.py
which writes conformal-breather.png. The style is the stamped copy manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

t, crest, step = np.loadtxt(HERE / "crest.csv", delimiter=",", unpack=True)
x0, eta0 = np.loadtxt(HERE / "initial.csv", delimiter=",", unpack=True)
x1, eta1, t1 = np.loadtxt(HERE / "peak.csv", delimiter=",", unpack=True)
t_peak = t1[0]
AMPLITUDE = 0.01                                  # et0 of run_Euler.m, in metres

hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig = plt.figure(figsize=(hs.WIDTH, 4.0))
grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0])
top = fig.add_subplot(grid[0, :])
left = fig.add_subplot(grid[1, 0])
right = fig.add_subplot(grid[1, 1], sharex=left)

for n, (x, eta, label) in enumerate([(x0, eta0, r"$t = 0$"),
                                      (x1, eta1, rf"$t = {t_peak:.1f}$ s, highest crest")]):
    order = np.argsort(x)
    style = hs.series(n)
    style["marker"] = None
    top.plot(x[order], eta[order], lw=0.8, label=label, **style)
top.set_xlim(min(x0.min(), x1.min()), max(x0.max(), x1.max()))
top.set_xlabel(r"$x$ (m)")
top.set_ylabel(r"$\eta$ (m)")
hs.panel(top, "(a) free surface")
hs.legend_outside(top, where="above", ncol=2)

left.plot(t, crest, color=hs.NAVY, lw=1.0)
left.axhline(AMPLITUDE, color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)))
left.set_xlim(t.min(), t.max())
left.set_xlabel(r"$t$ (s)")
left.set_ylabel(r"$\max_x \eta$ (m)")
hs.panel(left, "(b) highest crest; dotted: the amplitude of the carrier")

right.semilogy(t, step, color=hs.WINE, lw=1.0)
right.set_xlabel(r"$t$ (s)")
right.set_ylabel(r"time step (s)")
hs.grid(right)
hs.panel(right, "(c) adaptive time step")

out = HERE / "conformal-breather.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out, f"(highest crest {eta1.max():.4f} m at t = {t_peak:.2f} s)")
