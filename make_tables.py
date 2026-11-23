# -*- coding: utf-8 -*-
"""Tables 1-5 as formatted text in results/tables.xlsx (one sheet per table; the first row is the header)."""
import json
from pathlib import Path

import openpyxl

import lc_model as M

ROOT = Path(__file__).parent
S = json.load(open(ROOT / "results" / "study.json"))
P = M.P
f1 = lambda v: f"{v:.1f}"
f2 = lambda v: f"{v:.2f}"
T = {}

T["T1"] = [["Quantity", "Value", "Basis"],
    ["Cell format", "Prismatic NMC, 148 × 26.5 × 91 mm, 50 Ah", "Generic VDA PHEV2-size cell"],
    ["Module", "12 cells in series, bottom cooled", "Study design"],
    ["Cell mass / specific heat", f"{P['m_cell']} kg / {P['cp_cell']:.0f} J/(kg·K)", "Typical prismatic NMC values"],
    ["In-plane (vertical) conductivity k_z", f"{P['k_z']:.0f} W/(m·K)", "Effective value, jelly roll plus can [@chen96; @drake]"],
    ["Through-thickness conductivity", f"{P['k_t']:.1f} W/(m·K)", "Assumed, includes can walls; stack alone is lower [@drake]"],
    ["Resistance at 25 °C, R_0", f"{P['R0']*1e3:.1f} mΩ", "Effective charge resistance, assumed"],
    ["Resistance temperature constant B_R", f"{P['B_R']:.0f} K", "About −2.7 %/K near 25 °C"],
    ["TIM thickness / conductivity", f"1.0 mm / 1, 3, 6 W/(m·K)", "Gap-pad range"],
    ["Inter-cell pad", "1.5 mm, 0.15 W/(m·K)", "Compression pad"],
    ["Cold plate", "Al, 4 mm, 10 channels 10 × 2 mm", "Study design"],
    ["Coolant", "Ethylene glycol–water 50/50", "ρ 1065 kg/m³, c_p 3300 J/(kg·K), k 0.38 W/(m·K)"],
    ["Coolant viscosity", "3.4 mPa·s at 25 °C, −2.2 %/K", "Fit to tabulated EG50 data"],
    ["Flow rate", "0.5, 1, 2, 4, 8 L/min", "Design variable"],
    ["Supply temperature", "20–44 °C", "Design variable"],
    ["Charge window / C-rate", "10→80 % SOC / 1C, 2C, 3C", "Study design"],
    ["Hot-spot limit, derating band", "45 °C, 2 K", "Eq. (12)"],
    ["Pump efficiency", "0.5", "Assumed"],
]

H = [["Layout", "Flow (L/min)", "Velocity (m/s)", "Re", "Mean Nu", "Mean h (W/(m²·K))", "Δp (kPa)", "Pump power (W)"]]
for r in S["hydraulics"]:
    H.append(["Single pass" if r["layout"] == "single" else "U-pass", f"{r['Q']:g}", f"{r['V']:.3f}", f"{r['Re']:.0f}",
              f1(r["Nu"]), f"{r['h']:.0f}", f2(r["dP"] / 1e3), f"{r['P_pump']:.3f}"])
T["T3"] = H

D = S["design"]
X = [["Layout", "Flow (L/min)", "Supply T (°C)", "Peak T (°C)", "Cell-to-cell ΔT (K)", "Within-cell ΔT (K)", "Charge time (min)"]]
for lay in ("single", "upass"):
    for q in (0.5, 2.0, 8.0):
        for ts in (25, 35, 40):
            r = next(x for x in D if x["C"] == 3 and x["tim_k"] == 3.0 and x["layout"] == lay and x["Q"] == q and x["Ts"] == ts)
            X.append(["Single pass" if lay == "single" else "U-pass", f"{q:g}", str(ts), f1(r["Tmax"]), f1(r["dT_cells"]),
                      f1(r["dT_int"]), f1(r["t_charge"])])
T["T4"] = X

Y = [["Site", "Daytime T p50 / p99 / max (°C)", "Hours ≥ 40 °C (%)", "Radiator: > 15 min (%)", "Radiator: stalled (%)",
      "Chiller 30 °C: > 15 min (%)", "Chiller 25 °C: > 15 min (%)", "Chiller energy 30 / 25 °C (Wh/charge)"]]
for s, lab in [("jodhpur", "Jodhpur"), ("bhuj", "Bhuj"), ("anantapur", "Anantapur"), ("paradip", "Paradip")]:
    r = S["sites"][s]
    Y.append([lab, f"{r['Ta_p50']:.1f} / {r['Ta_p99']:.1f} / {r['Ta_max']:.1f}", f1(100 * r["share_ge40"]),
              f1(100 * r["radiator"]["share_extended"]), f1(100 * r["radiator"]["share_incomplete"]),
              f1(100 * r["chiller30"]["share_extended"]), f1(100 * r["chiller25"]["share_extended"]),
              f"{r['chiller30']['E_ch_mean']:.1f} / {r['chiller25']['E_ch_mean']:.1f}"])
T["T5"] = Y

V = S["verification"]
Z = [["Check", "Setting", "Result", "Note"]]
for n, tm, dti, cl in V["NZ"]:
    Z.append(["Layers per cell", f"NZ = {n}", f"Peak T {tm:.2f} °C; within-cell ΔT {dti:.2f} K", "dt = 2 s" + ("; used" if n == 8 else "")])
for d, tm, tc, cl in V["dt"]:
    Z.append(["Time step", f"dt = {d:g} s", f"Peak T {tm:.2f} °C; charge time {tc:.2f} min", "NZ = 8" + ("; used" if d == 2 else "")])
an = V["steady_analytic"]
for n, num, err, top in V["steady_numeric"]:
    if n in (4, 32):
        Z.append(["Steady conduction stack", f"NZ = {n}", f"Mean top-to-plate rise {num:.3f} K vs exact {an:.3f} K",
                  f"Largest cell deviation {err:.2f} K (pad exchange)"])
cl = max(abs(x[3]) for x in V["NZ"] + V["dt"])
Z.append(["Energy balance", "All runs above", f"|closure| < {max(cl, 1e-12):.0e}", "Relative to heat generated"])
T["T2"] = Z

wb = openpyxl.Workbook(); wb.remove(wb.active)
for k, rows in T.items():
    ws = wb.create_sheet(k)
    for r in rows:
        ws.append(r)
wb.save(ROOT / "results" / "tables.xlsx")
for k, rows in T.items():
    print(k); [print("  ", r) for r in rows[:30]]
