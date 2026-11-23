# -*- coding: utf-8 -*-
"""Fig. 1: module / cold-plate layout and thermal network (drawn, not simulated)."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

OUT = Path(__file__).parent / "figures"
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 10})

fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.6), gridspec_kw=dict(width_ratios=[1.35, 1.0, 1.0]))

# (a) side view of module on cold plate
ax = axs[0]
for i in range(12):
    ax.add_patch(Rectangle((i * 2.8, 1.0), 2.65, 9.1, fc="#cfd8e3", ec="k", lw=0.6))
    if i < 11:
        ax.add_patch(Rectangle((i * 2.8 + 2.65, 1.0), 0.15, 9.1, fc="#f2c14e", ec="none"))
ax.add_patch(Rectangle((0, 0.8), 33.5, 0.2, fc="#9b9b9b", ec="none"))
ax.add_patch(Rectangle((0, 0.0), 33.5, 0.8, fc="#b0b7bf", ec="k", lw=0.6))
ax.text(16.7, 0.35, "cold plate (Al, 4 mm)", ha="center", va="center", fontsize=8)
ax.text(1.3, 10.6, "cell 1", ha="center", fontsize=8)
ax.text(32.1, 10.6, "cell 12", ha="center", fontsize=8)
ax.text(16.7, 5.5, "12 prismatic NMC cells\n50 Ah, 148 × 26.5 × 91 mm", ha="center", va="center", fontsize=8.5,
        bbox=dict(fc="white", ec="none", alpha=0.85))
ax.annotate("TIM 1 mm", xy=(25, 0.9), xytext=(25, -1.6), fontsize=8, ha="center",
            arrowprops=dict(arrowstyle="-", lw=0.6))
ax.annotate("pad 1.5 mm", xy=(8.2, 8.0), xytext=(11.5, 11.8), fontsize=8,
            arrowprops=dict(arrowstyle="-", lw=0.6))
ax.add_patch(FancyArrowPatch((-3.5, 0.4), (-0.2, 0.4), arrowstyle="-|>", mutation_scale=10, color="#1f5fa8"))
ax.text(-3.6, -0.6, "inlet", fontsize=8, color="#1f5fa8")
ax.set_xlim(-4.5, 35); ax.set_ylim(-2.4, 12.6); ax.axis("off")
ax.set_title("(a) Module on cold plate (side view)", fontsize=10)

# (b) plan view of the two layouts
ax = axs[1]
def chans(y0, n, col):
    for j in range(n):
        ax.add_patch(Rectangle((0.5, y0 + j * 0.42), 8.0, 0.22, fc=col, ec="none"))
ax.add_patch(Rectangle((0.3, 5.3), 8.4, 4.4, fc="#eef1f4", ec="k", lw=0.6))
chans(5.45, 10, "#7fb3e0")
ax.annotate("", xy=(8.9, 7.5), xytext=(-0.6, 7.5), arrowprops=dict(arrowstyle="-|>", color="#1f5fa8", lw=1.2))
ax.text(4.5, 10.0, "Single pass: 10 channels in parallel", ha="center", fontsize=8.5)
ax.add_patch(Rectangle((0.3, 0.0), 8.4, 4.4, fc="#eef1f4", ec="k", lw=0.6))
chans(0.15, 5, "#7fb3e0"); chans(2.25, 5, "#e3a07a")
ax.annotate("", xy=(8.6, 1.0), xytext=(-0.6, 1.0), arrowprops=dict(arrowstyle="-|>", color="#1f5fa8", lw=1.2))
ax.annotate("", xy=(-0.6, 3.2), xytext=(8.6, 3.2), arrowprops=dict(arrowstyle="-|>", color="#c0532b", lw=1.2))
ax.annotate("", xy=(8.95, 3.2), xytext=(8.95, 1.0),
            arrowprops=dict(arrowstyle="-", color="0.3", lw=1.0, connectionstyle="arc3,rad=0.6"))
ax.text(4.5, 4.7, "U-pass: 5 channels out, 5 back", ha="center", fontsize=8.5)
ax.text(4.5, -0.7, "channels 10 × 2 mm, length 336 mm per pass", ha="center", fontsize=8)
ax.set_xlim(-1.2, 10); ax.set_ylim(-1.2, 10.6); ax.axis("off")
ax.set_title("(b) Channel layouts (plan view)", fontsize=10)

# (c) thermal network of one cell
ax = axs[2]
nodes = [(2.0, 9.0, "layer NZ (top)"), (2.0, 7.0, "..."), (2.0, 5.0, "layer 1 (bottom)"), (2.0, 3.0, "plate node"),
         (2.0, 1.0, "coolant segment(s)")]
for x, y, t in nodes:
    if t != "...":
        ax.add_patch(Rectangle((x - 0.25, y - 0.25), 0.5, 0.5, fc="white", ec="k", lw=0.8))
    ax.text(x + 0.5, y, t, va="center", fontsize=8.5)
for (y1, y2, lab) in [(8.75, 7.4, r"$k_z A/\Delta z$"), (6.6, 5.25, ""), (4.75, 3.25, "TIM + half layer"),
                      (2.75, 1.25, r"$\dot{m}c_p(1-e^{-hA/\dot{m}c_p})$")]:
    ax.plot([2.0, 2.0], [y1, y2], color="k", lw=0.8)
    if lab:
        ax.text(1.75, (y1 + y2) / 2, lab.replace("e^{−hA/ṁc_p}", "e$^{-hA/\\dot{m}c_p}$"), ha="right", va="center", fontsize=8)
ax.annotate("", xy=(0.6, 9.0), xytext=(1.75, 9.0), arrowprops=dict(arrowstyle="-", lw=0.8, ls="--"))
ax.text(0.5, 9.0, "air", ha="right", va="center", fontsize=8)
ax.annotate("neighbour\ncells (pads)", xy=(2.25, 5.6), xytext=(3.5, 6.3), fontsize=8,
            arrowprops=dict(arrowstyle="<->", lw=0.6))
ax.text(3.35, 9.6, r"$I^2R/N_Z$ in every layer", fontsize=8, style="italic")
ax.set_xlim(-1.6, 6.5); ax.set_ylim(0.2, 10.4); ax.axis("off")
ax.set_title("(c) Thermal network of one cell", fontsize=10)

fig.tight_layout()
fig.savefig(OUT / "fig1_layout.png", dpi=300)
print("saved")
