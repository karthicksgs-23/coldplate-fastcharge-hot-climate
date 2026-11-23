# -*- coding: utf-8 -*-
"""Cold-plate liquid cooling of a 12-cell prismatic NMC module during DC fast charging.

Thermal network (backward Euler, coolant segments quasi-steady):
  * every cell is split into NZ horizontal layers (vertical conduction, uniform Joule heat),
  * the bottom layer is linked through the thermal interface material (TIM) to one plate node per cell,
  * plate nodes conduct along the module and give heat to the coolant segment(s) under the cell,
  * each coolant segment uses the exact constant-wall-temperature solution
      T_out = a T_in + (1 - a) T_plate,  a = exp(-hA / (m_dot cp)).
Layouts: 'single' = 10 parallel channels, one pass under cells 1..12;
         'upass'  = 5 channels out under cells 1..12 and 5 channels back (U-turn), same total area.
Heat transfer: local laminar thermal-entry Nusselt number (Shah & London, H boundary) shifted to the
rectangular-duct asymptote; Gnielinski above Re = 2300.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")          # small dense solves are much faster single-threaded
import numpy as np

P = dict(
    # cell (generic VDA PHEV2-size prismatic NMC, not a measured product)
    n_cells=12, cap_Ah=50.0, m_cell=0.90, cp_cell=1050.0,
    H=0.091, W=0.148, t_cell=0.0265, pitch=0.028,
    k_z=25.0, k_t=1.0,                 # in-plane (vertical) / through-thickness conductivity, W/(m K)
    pad_t=0.0015, pad_k=0.15,          # compression pad between cells
    R0=1.2e-3, B_R=2400.0, T_ref=298.15,
    # interface and cold plate
    tim_t=1.0e-3, tim_k=3.0,
    plate_k=167.0, plate_rho=2700.0, plate_cp=900.0, plate_base=0.002, plate_t=0.004,
    n_ch=10, ch_w=0.010, ch_h=0.002,
    # coolant: ethylene glycol / water 50/50
    rho_c=1065.0, cp_c=3300.0, k_c=0.38, mu25=3.4e-3, mu_slope=0.022,
    Nu_fd=5.39, fRe=76.3,              # rectangular duct (aspect 0.2), one wide wall heated; Darcy f*Re
    K_manifold=1.5, K_dev=1.2, K_bend=1.5, eta_pump=0.5,
    h_air=5.0,
    # charging
    soc0=0.10, soc1=0.80, T_lim=45.0, derate_band=2.0, t_cap=90 * 60.0,
)


def mu_c(T, p=P):
    return p["mu25"] * np.exp(-p["mu_slope"] * (T - 25.0))


def nu_circ_H(xs):
    """Local Nusselt number, thermally developing laminar flow, circular duct, H boundary (Shah & London)."""
    xs = np.maximum(xs, 1e-7)
    return np.where(xs <= 5e-5, 3.302 * xs ** (-1 / 3) - 1.0,
           np.where(xs <= 1.5e-3, 1.302 * xs ** (-1 / 3) - 0.5,
                    4.364 + 8.68 * (1e3 * xs) ** (-0.506) * np.exp(-41.0 * xs)))


def hydraulics(Q_lpm, layout, T_mean=25.0, p=P):
    """Per-channel velocity, Re, pressure drop and pump power for one layout and flow rate."""
    n_par = p["n_ch"] if layout == "single" else p["n_ch"] // 2
    A = p["ch_w"] * p["ch_h"]
    Dh = 2 * p["ch_w"] * p["ch_h"] / (p["ch_w"] + p["ch_h"])
    Q = Q_lpm / 60000.0
    V = Q / (n_par * A)
    mu = mu_c(T_mean, p)
    Re = p["rho_c"] * V * Dh / mu
    Pr = p["cp_c"] * mu / p["k_c"]
    L_pass = p["n_cells"] * p["pitch"]
    n_pass = 1 if layout == "single" else 2
    fD = p["fRe"] / Re if Re < 2300 else (0.79 * np.log(Re) - 1.64) ** -2
    K = p["K_manifold"] + n_pass * p["K_dev"] + (p["K_bend"] if n_pass == 2 else 0.0)
    dP = (fD * n_pass * L_pass / Dh + K) * 0.5 * p["rho_c"] * V ** 2
    return dict(V=V, Re=Re, Pr=Pr, Dh=Dh, dP=dP, P_pump=dP * Q / p["eta_pump"], n_par=n_par,
                n_pass=n_pass, L_pass=L_pass, mdot=Q * p["rho_c"])


def seg_hA(hyd, p=P):
    """Segment hA (W/K) for each cell position along one pass, plus mean Nu."""
    x_mid = (np.arange(p["n_cells"]) + 0.5) * p["pitch"]
    Re, Pr, Dh = hyd["Re"], hyd["Pr"], hyd["Dh"]
    nu_lam = p["Nu_fd"] + (nu_circ_H(x_mid / (Dh * Re * Pr)) - 4.364)
    if Re >= 2300:
        f = (0.79 * np.log(Re) - 1.64) ** -2
        nu_t = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * np.sqrt(f / 8) * (Pr ** (2 / 3) - 1))
        nu = np.maximum(nu_lam, nu_t)
    else:
        nu = nu_lam
    h = nu * p["k_c"] / Dh
    perim = p["ch_w"] + 2 * p["ch_h"]          # top wall + two webs (fin efficiency ~1 for Al webs)
    return h * hyd["n_par"] * perim * p["pitch"], nu.mean(), h.mean()


def simulate(C_rate, Q_lpm, layout, T_supply, T0=None, T_amb=None, tim_k=None, NZ=8, dt=2.0,
             p=P, scale=None, record=False):
    """Fast charge SOC 0.10 -> 0.80 at C_rate with temperature derating. Returns summary (and traces)."""
    p = dict(p)
    if tim_k is not None:
        p["tim_k"] = tim_k
    sc = dict(k_z=1, R0=1, cp_cell=1, tim_k=1, h=1, B_R=1, pad_k=1)
    if scale:
        sc.update(scale)
    for k in ("k_z", "R0", "cp_cell", "tim_k", "B_R", "pad_k"):
        p[k] = p[k] * sc[k]
    T0 = T_supply if T0 is None else T0
    T_amb = T0 if T_amb is None else T_amb
    n = p["n_cells"]
    A_b = p["W"] * p["t_cell"]
    dz = p["H"] / NZ
    Cl = p["m_cell"] * p["cp_cell"] / NZ
    Gv = p["k_z"] * A_b / dz
    Gb = 1.0 / (dz / 2 / (p["k_z"] * A_b) + p["tim_t"] / (p["tim_k"] * A_b))
    A_face = p["W"] * p["H"] / NZ
    Glat = 1.0 / (p["t_cell"] / (p["k_t"] * A_face) + p["pad_t"] / (p["pad_k"] * A_face))
    Gair_top = p["h_air"] * A_b
    Gair_end = p["h_air"] * A_face
    plate_vol = p["W"] * p["pitch"] * p["plate_t"] - p["n_ch"] * p["ch_w"] * p["ch_h"] * p["pitch"]
    Cp_pl = plate_vol * p["plate_rho"] * p["plate_cp"]
    A_pl = p["W"] * p["plate_base"] + (p["W"] - p["n_ch"] * p["ch_w"]) * p["ch_h"]
    Gpl = p["plate_k"] * A_pl / p["pitch"]

    hyd = hydraulics(Q_lpm, layout, T_supply, p)
    n_pass = hyd["n_pass"]
    # unknown indexing
    nc = n * NZ
    ip = nc                       # plate nodes
    ic = nc + n                   # coolant segment outlets
    nseg = n * n_pass
    N = nc + n + nseg
    seg_cell = np.arange(n) if n_pass == 1 else np.r_[np.arange(n), np.arange(n)[::-1]]
    cidx = lambda i, k: i * NZ + k

    # constant part of the conduction matrix (K T = loads), capacitance vector
    Kc = np.zeros((N, N)); cap = np.zeros(N); bair = np.zeros(N)
    def link(a, b, G):
        Kc[a, a] += G; Kc[b, b] += G; Kc[a, b] -= G; Kc[b, a] -= G
    for i in range(n):
        for k in range(NZ):
            cap[cidx(i, k)] = Cl
            if k < NZ - 1: link(cidx(i, k), cidx(i, k + 1), Gv)
            if i < n - 1: link(cidx(i, k), cidx(i + 1, k), Glat)
        link(cidx(i, 0), ip + i, Gb)
        Kc[cidx(i, NZ - 1), cidx(i, NZ - 1)] += Gair_top; bair[cidx(i, NZ - 1)] += Gair_top
        cap[ip + i] = Cp_pl
        if i < n - 1: link(ip + i, ip + i + 1, Gpl)
    for i in (0, n - 1):
        for k in range(NZ):
            Kc[cidx(i, k), cidx(i, k)] += Gair_end; bair[cidx(i, k)] += Gair_end

    T = np.full(N, float(T0)); T[ic:] = T_supply
    soc, t = p["soc0"], 0.0
    I_set = C_rate * p["cap_Ah"]
    Tm = T_supply
    E_gen = E_cool = E_air = 0.0
    E0 = (cap * T).sum()
    rec = []
    Tmax = T0; dT_cells = 0.0; dT_int = 0.0; derated = False; t_derate = 0.0; q_peak = 0.0
    f_tim = (p["tim_t"] / (p["tim_k"] * A_b)) / (dz / 2 / (p["k_z"] * A_b) + p["tim_t"] / (p["tim_k"] * A_b))
    while soc < p["soc1"] and t < p["t_cap"]:
        # hydraulics/heat transfer with lagged mean coolant temperature
        hyd = hydraulics(Q_lpm, layout, Tm, p)
        hA, _, _ = seg_hA(hyd, p)
        hA = hA * sc["h"]
        mcp = hyd["mdot"] * p["cp_c"]
        a = np.exp(-np.r_[hA] / mcp) if n_pass == 1 else np.exp(-np.r_[hA, hA] / mcp)
        # current from previous-step hot spot
        Thot = T[:nc].max()
        I = I_set * float(np.clip((p["T_lim"] - Thot) / p["derate_band"], 0.0, 1.0))
        if I < I_set: derated = True; t_derate += dt
        Tcell = T[:nc].reshape(n, NZ)
        R = p["R0"] * np.exp(p["B_R"] * (1 / (Tcell + 273.15) - 1 / p["T_ref"]))
        q = np.zeros(N); q[:nc] = (I ** 2 * R / NZ).ravel()
        # assemble implicit system
        M = Kc.copy(); rhs = q + bair * T_amb
        M[np.arange(N), np.arange(N)] += cap / dt
        rhs = rhs + cap / dt * T
        for s in range(nseg):
            i = seg_cell[s]; row = ic + s
            M[row, :] = 0.0
            M[row, row] = 1.0
            M[row, ip + i] = -(1 - a[s])
            if s == 0:
                rhs[row] = a[s] * T_supply
            else:
                M[row, row - 1] = -a[s]; rhs[row] = 0.0
            # plate heat loss = mcp (T_out - T_in) = mcp (1-a)(Tp - T_in)
            M[ip + i, ip + i] += mcp * (1 - a[s])
            if s == 0:
                rhs[ip + i] += mcp * (1 - a[s]) * T_supply
            else:
                M[ip + i, row - 1] -= mcp * (1 - a[s])
        Tn = np.linalg.solve(M, rhs)
        T_in = np.r_[T_supply, Tn[ic:ic + nseg - 1]]
        q_cool = mcp * (Tn[ic:ic + nseg] - T_in)
        E_gen += q.sum() * dt; E_cool += q_cool.sum() * dt
        E_air += (bair * (Tn - T_amb)).sum() * dt
        T = Tn; t += dt
        soc += I * dt / 3600.0 / p["cap_Ah"]
        Tm = 0.5 * (T_supply + T[ic + nseg - 1])
        Tc3 = T[:nc].reshape(n, NZ)
        Tmax = max(Tmax, Tc3.max())
        dT_cells = max(dT_cells, Tc3.mean(1).max() - Tc3.mean(1).min())
        T_face = T[ip:ip + n] + (Tc3[:, 0] - T[ip:ip + n]) * f_tim     # cell bottom-face temperature
        dT_int = max(dT_int, (Tc3.max(1) - T_face).max())
        q_peak = max(q_peak, q_cool.sum())
        if record:
            rec.append((t, I, soc, Tc3.max(), Tc3.mean(1).max() - Tc3.mean(1).min(), T[ic + nseg - 1],
                        *Tc3.max(1), q_cool.sum(), q.sum()))
    closure = (E_gen - E_cool - E_air - ((cap * T).sum() - E0)) / max(E_gen, 1e-9)
    hyd = hydraulics(Q_lpm, layout, Tm, p)
    hA, nu_mean, h_mean = seg_hA(hyd, p)
    out = dict(C=C_rate, Q=Q_lpm, layout=layout, Ts=T_supply, T0=T0, tim_k=p["tim_k"],
               t_charge=t / 60.0, complete=bool(soc >= p["soc1"]), soc_end=soc,
               Tmax=Tmax, dT_cells=dT_cells, dT_int=dT_int, derated=derated, t_derate=t_derate / 60.0,
               E_gen=E_gen / 3600.0, E_cool=E_cool / 3600.0, closure=closure,
               Re=hyd["Re"], V=hyd["V"], dP=hyd["dP"], P_pump=hyd["P_pump"], Nu=nu_mean, h=h_mean,
               E_pump=hyd["P_pump"] * t / 3600.0, q_peak=q_peak)
    if record:
        out["trace"] = np.array(rec); out["T_final"] = T.copy(); out["NZ"] = NZ
    return out


def chiller_cop(T_supply, T_amb, eta=0.45, cap=6.0):
    """Vapour-compression COP: fraction of Carnot between evaporator (supply - 5 K) and condenser (ambient + 10 K)."""
    Te = T_supply - 5.0 + 273.15
    Tc = T_amb + 10.0 + 273.15
    return min(cap, eta * Te / max(Tc - Te, 1.0))


if __name__ == "__main__":
    import time
    t0 = time.time()
    for lay in ("single", "upass"):
        r = simulate(3, 2.0, lay, 25.0)
        print(lay, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    print("elapsed", round(time.time() - t0, 2))
