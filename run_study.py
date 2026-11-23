# -*- coding: utf-8 -*-
"""Run every simulation behind the study and write results/study.json + results/traces.npz."""
import json
import time
from itertools import product
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import lc_model as M

ROOT = Path(__file__).parent
RES = ROOT / "results"
SITES = ["jodhpur", "bhuj", "anantapur", "paradip"]
REF = dict(C_rate=3, Q_lpm=2.0, layout="upass", tim_k=3.0)


def _run(kw):
    r = M.simulate(**kw)
    r.pop("trace", None)
    return {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in r.items()}


def design_matrix(pool):
    cases = [dict(C_rate=c, Q_lpm=q, layout=l, tim_k=k, T_supply=ts)
             for c, q, l, k, ts in product([1, 2, 3], [0.5, 1, 2, 4, 8], ["single", "upass"],
                                           [1.0, 3.0, 6.0], [20, 25, 30, 35, 40])]
    return pool.map(_run, cases)


def climate_grid(pool):
    """Reference design, 3C, T0 = ambient. Policies: radiator (supply = ambient + 5 K), chiller 25 / 30 °C."""
    Ta = np.arange(15, 51, 1.0)
    out = {}
    for pol in ("radiator", "chiller25", "chiller30"):
        cases = []
        for ta in Ta:
            ts = ta + 5.0 if pol == "radiator" else (25.0 if pol == "chiller25" else 30.0)
            cases.append(dict(REF, T_supply=ts, T0=float(ta), T_amb=float(ta)))
        rs = pool.map(_run, cases)
        for ta, r in zip(Ta, rs):
            ts = r["Ts"]
            cop = M.chiller_cop(ts, ta) if pol != "radiator" else np.inf
            r["E_chiller"] = 0.0 if pol == "radiator" else max(r["E_cool"], 0.0) / cop
            r["COP"] = None if pol == "radiator" else cop
            r["Ta"] = float(ta)
        out[pol] = rs
    return out


def site_climate(grid):
    """Map every daytime hour (08:00-20:00 IST) of 2005-2023 onto the climate grid."""
    res = {}
    for s in SITES:
        d = np.load(ROOT / "data" / f"climate_{s}.npz")
        hour_ist = (d["hour_utc"] + 5.5) % 24
        day = (hour_ist >= 8) & (hour_ist < 20)
        Ta = d["T"][day]
        r = dict(n_hours=int(day.sum()), Ta_p50=float(np.percentile(Ta, 50)), Ta_p95=float(np.percentile(Ta, 95)),
                 Ta_p99=float(np.percentile(Ta, 99)), Ta_max=float(Ta.max()),
                 share_ge35=float((Ta >= 35).mean()), share_ge40=float((Ta >= 40).mean()))
        for pol, rows in grid.items():
            ta_g = np.array([x["Ta"] for x in rows])
            tc = np.interp(Ta, ta_g, [x["t_charge"] for x in rows])
            ext = np.interp(Ta, ta_g, [1.0 if x["t_charge"] > 15.0 else 0.0 for x in rows])
            qpk = np.interp(Ta, ta_g, [x["q_peak"] for x in rows])
            ech = np.interp(Ta, ta_g, [x["E_chiller"] for x in rows])
            comp = np.interp(Ta, ta_g, [1.0 if x["complete"] else 0.0 for x in rows])
            r[pol] = dict(t_p50=float(np.percentile(tc, 50)), t_p95=float(np.percentile(tc, 95)),
                          t_max=float(tc.max()), share_extended=float((ext >= 0.5).mean()), q_peak_p99=float(np.percentile(qpk, 99)),
                          share_incomplete=float((tc >= 89.9).mean()),
                          E_ch_mean=float(ech.mean()), E_ch_p95=float(np.percentile(ech, 95)))
        hist, edges = np.histogram(Ta, bins=np.arange(0, 52, 2))
        r["hist"] = hist.tolist(); r["edges"] = edges.tolist()
        res[s] = r
    return res


def sensitivity(pool):
    keys = ["k_z", "R0", "cp_cell", "tim_k", "h", "B_R", "pad_k"]
    cases = []
    for ts in (25.0, 35.0):
        cases.append(dict(REF, T_supply=ts))
        for k, f in product(keys, (0.8, 1.2)):
            cases.append(dict(REF, T_supply=ts, scale={k: f}))
    rs = pool.map(_run, cases)
    out = {}
    for j, ts in enumerate((25.0, 35.0)):
        blk = rs[j * (1 + 2 * len(keys)):(j + 1) * (1 + 2 * len(keys))]
        base = blk[0]
        out[str(int(ts))] = dict(base=base, cases={k: dict(lo=blk[1 + 2 * i], hi=blk[2 + 2 * i]) for i, k in enumerate(keys)})
    return out


def hydraulics_table():
    rows = []
    for lay, q in product(["single", "upass"], [0.5, 1, 2, 4, 8]):
        h = M.hydraulics(q, lay, 25.0)
        _, nu, hm = M.seg_hA(h)
        rows.append(dict(layout=lay, Q=q, V=h["V"], Re=h["Re"], Nu=float(nu), h=float(hm), dP=h["dP"], P_pump=h["P_pump"]))
    return rows


def R_conv_ref():
    """Plate-to-coolant resistance of one cell for the reference design (2 L/min, U-pass, 25 °C)."""
    h = M.hydraulics(2.0, "upass", 25.0)
    hA, _, _ = M.seg_hA(h)
    mcp = h["mdot"] * M.P["cp_c"]
    G = mcp * (1 - np.exp(-hA / mcp))          # one segment per pass
    return float(1.0 / (G + G[::-1]).mean())


def verification(pool):
    """Grid and time-step independence, energy closure and a steady analytic check."""
    nz = pool.map(_run, [dict(REF, T_supply=25.0, NZ=n) for n in (4, 8, 16, 32)])
    dts = pool.map(_run, [dict(REF, T_supply=25.0, dt=d) for d in (0.5, 1.0, 2.0, 5.0)])
    # steady state: constant R (B_R = 0), no air loss, no derating, 3C heat held for 3 h.
    # Exact steady solution of the conduction stack: T_top - T_plate = q (R_TIM + H / (2 k_z A)).
    p = dict(M.P, h_air=0.0, T_lim=1e3, soc1=1e9, t_cap=3 * 3600.0)
    A = p["W"] * p["t_cell"]
    q = (3 * p["cap_Ah"]) ** 2 * p["R0"]
    R_int = p["H"] / (2 * p["k_z"] * A); R_tim = p["tim_t"] / (p["tim_k"] * A)
    steady = []
    for n in (4, 8, 16, 32):
        r = M.simulate(3, 2.0, "upass", 25.0, NZ=n, dt=20.0, p=p, scale={"B_R": 0.0}, record=True)
        T = r["T_final"]; nc = 12 * n
        top = T[:nc].reshape(12, n)[:, -1]; plate = T[nc:nc + 12]
        num = top - plate
        steady.append((n, float(num.mean()), float(np.abs(num - q * (R_int + R_tim)).max()), float(top.max())))
    T_an = float(q * (R_int + R_tim))
    return dict(NZ=[(n, r["Tmax"], r["dT_int"], r["closure"]) for n, r in zip((4, 8, 16, 32), nz)],
                dt=[(d, r["Tmax"], r["t_charge"], r["closure"]) for d, r in zip((0.5, 1.0, 2.0, 5.0), dts)],
                steady_numeric=steady, steady_analytic=T_an,
                steady_terms=dict(q=q, R_int=R_int, R_tim=R_tim, R_conv_ref=R_conv_ref()))


def traces():
    tr = {}
    for lay, ts in product(["single", "upass"], [25.0, 35.0]):
        r = M.simulate(3, 2.0, lay, ts, record=True)
        tr[f"{lay}_{int(ts)}"] = r["trace"]
    for q in (0.5, 8.0):
        tr[f"upass_25_Q{q}"] = M.simulate(3, q, "upass", 25.0, record=True)["trace"]
    return tr


if __name__ == "__main__":
    t0 = time.time()
    with Pool(8) as pool:
        out = dict(params=M.P, ref=REF)
        out["design"] = design_matrix(pool); print("design", len(out["design"]), round(time.time() - t0))
        grid = climate_grid(pool); out["climate_grid"] = grid; print("grid", round(time.time() - t0))
        out["sites"] = site_climate(grid)
        out["sensitivity"] = sensitivity(pool); print("sens", round(time.time() - t0))
        out["hydraulics"] = hydraulics_table()
        out["verification"] = verification(pool); print("verif", round(time.time() - t0))
    np.savez_compressed(RES / "traces.npz", **traces())
    json.dump(out, open(RES / "study.json", "w"), indent=1, default=float)
    print("done", round(time.time() - t0), "s")
