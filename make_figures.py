# -*- coding: utf-8 -*-
"""Figures 2-6 from results/study.json, results/ts_sweep.json and results/traces.npz."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).parent
FIG = ROOT / "figures"
S = json.load(open(ROOT / "results" / "study.json"))
SW = json.load(open(ROOT / "results" / "ts_sweep.json"))
TR = np.load(ROOT / "results" / "traces.npz")
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
# colour-blind-safe (Okabe-Ito)
C1, C2, C3, C4, C5 = "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"
SITES = [("jodhpur", "Jodhpur"), ("bhuj", "Bhuj"), ("anantapur", "Anantapur"), ("paradip", "Paradip")]
SC = dict(jodhpur=C2, bhuj=C5, anantapur=C3, paradip=C1)


def limit(ax):
    ax.axhline(45, color="k", ls=":", lw=1, label="45 °C charge limit")


# ---------- Fig. 2: time traces, 3C, 2 L/min ----------
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.4))
for key, lab, col, ls in [("single_25", "Single pass, 25 °C", C1, "-"), ("upass_25", "U-pass, 25 °C", C3, "-"),
                          ("single_35", "Single pass, 35 °C", C1, "--"), ("upass_35", "U-pass, 35 °C", C3, "--")]:
    tr = TR[key]; tm = tr[:, 0] / 60
    axs[0].plot(tm, tr[:, 3], color=col, ls=ls, label=lab)
    axs[1].plot(tm, tr[:, 1], color=col, ls=ls)
    axs[2].plot(tm, tr[:, 4], color=col, ls=ls)
limit(axs[0])
axs[0].set_ylabel("Hot-spot temperature (°C)"); axs[1].set_ylabel("Charge current (A)")
axs[2].set_ylabel("Cell-to-cell ΔT (K)")
for ax, t in zip(axs, "abc"):
    ax.set_xlabel("Time (min)"); ax.set_title(f"({t})", loc="left", fontsize=10)
axs[0].set_ylim(24, 47)
fig.legend(*axs[0].get_legend_handles_labels(), loc="upper center", ncol=5, frameon=False, fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.9))
fig.savefig(FIG / "fig2_traces.png", dpi=300); plt.close(fig)

# ---------- Fig. 3: flow rate, 3C, Ts = 25 °C, TIM 3 ----------
D = [r for r in S["design"] if r["C"] == 3 and r["Ts"] == 25 and r["tim_k"] == 3.0]
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3))
for lay, lab, col, mk in [("single", "Single pass", C1, "o"), ("upass", "U-pass", C3, "s")]:
    rr = sorted([r for r in D if r["layout"] == lay], key=lambda r: r["Q"])
    q = [r["Q"] for r in rr]
    axs[0].plot(q, [r["Tmax"] for r in rr], marker=mk, color=col, label=lab)
    axs[1].plot(q, [r["dT_cells"] for r in rr], marker=mk, color=col)
    axs[1].plot(q, [r["dT_int"] for r in rr], marker=mk, color=col, ls="--", mfc="white")
    axs[2].plot(q, [r["P_pump"] for r in rr], marker=mk, color=col)
axs[1].plot([], [], color="0.3", label="cell-to-cell"); axs[1].plot([], [], color="0.3", ls="--", label="within cell")
axs[1].legend(frameon=False, fontsize=8.5)
axs[0].set_ylabel("Peak hot-spot temperature (°C)"); axs[1].set_ylabel("Peak temperature difference (K)")
axs[2].set_ylabel("Pump power per module (W)"); axs[2].set_yscale("log")
for ax, t in zip(axs, "abc"):
    ax.set_xscale("log", base=2); ax.set_xticks([0.5, 1, 2, 4, 8]); ax.set_xticklabels(["0.5", "1", "2", "4", "8"])
    ax.set_xlabel("Coolant flow rate (L/min)"); ax.set_title(f"({t})", loc="left", fontsize=10)
axs[0].legend(frameon=False)
fig.tight_layout(); fig.savefig(FIG / "fig3_flow.png", dpi=300); plt.close(fig)

# ---------- Fig. 4: supply temperature sweep ----------
fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.5))
for c, ls in [(3, "-"), (2, "--")]:
    for k, col in [(1.0, C2), (3.0, C1), (6.0, C3)]:
        rr = sorted([r for r in SW if r["C"] == c and r["tim_k"] == k], key=lambda r: r["Ts"])
        lab = f"{c}C, TIM {k:g} W/(m·K)"
        axs[0].plot([r["Ts"] for r in rr], [r["t_charge"] for r in rr], color=col, ls=ls, label=lab)
        axs[1].plot([r["Ts"] for r in rr], [r["Tmax"] for r in rr], color=col, ls=ls)
limit(axs[1])
axs[0].set_ylabel("Charge time, 10→80 % SOC (min)"); axs[1].set_ylabel("Peak hot-spot temperature (°C)")
for ax, t in zip(axs, "ab"):
    ax.set_xlabel("Coolant supply temperature (°C)"); ax.set_title(f"({t})", loc="left", fontsize=10)
axs[0].legend(frameon=False, fontsize=8, ncol=2)
axs[1].legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig(FIG / "fig4_supply_temp.png", dpi=300); plt.close(fig)

# ---------- Fig. 5: climate ----------
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.4))
for s, lab in SITES:
    r = S["sites"][s]; e = np.array(r["edges"]); h = np.array(r["hist"]) / r["n_hours"] * 100
    axs[0].step(e[:-1], h, where="post", color=SC[s], label=lab)
axs[0].set_xlabel("Daytime air temperature (°C)"); axs[0].set_ylabel("Share of daytime hours (%)")
axs[0].set_xlim(10, 50); axs[0].legend(frameon=False, fontsize=8.5)
for pol, lab, col in [("radiator", "Radiator (ambient + 5 K)", C2), ("chiller30", "Chiller, 30 °C supply", C5),
                      ("chiller25", "Chiller, 25 °C supply", C1)]:
    g = S["climate_grid"][pol]
    ta = [x["Ta"] for x in g]
    axs[1].plot(ta, [min(x["t_charge"], 60) for x in g], color=col, label=lab)
    if pol != "radiator":
        axs[2].plot(ta, [x["E_chiller"] for x in g], color=col, label=lab)
axs[1].set_ylim(10, 62); axs[1].text(41.0, 48, "charging stalls\n(supply ≥ 45 °C)", fontsize=8, color=C2)
axs[1].set_xlabel("Ambient temperature (°C)"); axs[1].set_ylabel("Charge time at 3C (min)")
axs[2].set_xlabel("Ambient temperature (°C)"); axs[2].set_ylabel("Chiller electricity per charge (Wh)")
axs[1].legend(frameon=False, fontsize=8.5, loc="upper left")
for ax, t in zip(axs, "abc"):
    ax.set_title(f"({t})", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig(FIG / "fig5_climate.png", dpi=300); plt.close(fig)

# ---------- Fig. 6: sensitivity ----------
NAMES = dict(R0="Cell resistance R$_0$", cp_cell="Cell heat capacity", k_z="In-plane conductivity k$_z$",
             h="Convective coefficient h", B_R="Resistance temperature coefficient", tim_k="TIM conductivity",
             pad_k="Pad conductivity")
fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.3))
for ax, ts, met, xl in [(axs[0], "25", "Tmax", "Change in peak temperature (K), 25 °C supply"),
                        (axs[1], "35", "t_charge", "Change in charge time (min), 35 °C supply")]:
    blk = S["sensitivity"][ts]; b = blk["base"][met]
    items = sorted(blk["cases"].items(), key=lambda kv: abs(kv[1]["hi"][met] - kv[1]["lo"][met]))
    for i, (k, c) in enumerate(items):
        lo, hi = c["lo"][met] - b, c["hi"][met] - b
        ax.barh(i, lo, color=C1, height=0.6, label="−20 %" if i == 0 else None)
        ax.barh(i, hi, color=C2, height=0.6, label="+20 %" if i == 0 else None)
    ax.set_yticks(range(len(items))); ax.set_yticklabels([NAMES[k] for k, _ in items], fontsize=9)
    ax.axvline(0, color="k", lw=0.8); ax.set_xlabel(xl)
axs[0].set_title("(a)", loc="left", fontsize=10); axs[1].set_title("(b)", loc="left", fontsize=10)
axs[1].legend(frameon=False, fontsize=8.5, loc="lower right")
fig.tight_layout(); fig.savefig(FIG / "fig6_sensitivity.png", dpi=300); plt.close(fig)
print("figures saved")
