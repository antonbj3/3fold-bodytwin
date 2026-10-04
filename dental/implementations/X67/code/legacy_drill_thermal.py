"DENT-PROC-DRILL-THERMAL — transient heat in bone during incremental implant drilling (K3, 2026-09-23).\n\nFree standing cell. Axisymmetric (r,z) finit-volym, explicit tidsstegning (numba), med\n  * Moving heat source at the drill point: P_ben = eta * u(f_rev) * v * A_cut * phi(z)       [W]\n      u(f_rev) = U0 * (f_rev/F0)^(-M_SIZE)   Specific cutting energy with size effect (Kienzle-form)\n      A_cut    = pi/4 (D^2 - d_prev^2)       ring surface cut (stegvis borr) ; phi = 1 kortikalt, BV/TV spongious\n  * material removal: cells inside the drill contour (cylinder + 118°-spetskon) removed -> their heat goes with the shavings,\n    The hole wall becomes Robin-rand (spolning h_irr eller luft h_air) — the change in topology is visible in the heat problem\n  * Pennes-perfusion w_b*rho_b*c_b*(T_a - T) (in vivo), ingen perfusion in vitro\n  * layer: cortical crown layer t_c over spongious leg (Effective properties: bone + marrow, Hashin–Shtrikman boundaries)\n  * utdata: T(t) i sonder, T_max-fhigh, CEM43- field (Sapareto & Dewey 1984, PMID 6547421), time over 47 °C\nEnheter: indata i mm, s, rpm, °C; internt SI (m, s, W).\nParameters and sources: see PARAMS below and results/K3_drilling/PREREG.md.\nIngen kod kopierad. Formen P = cutting -> heat source follows COMPUTE_CELL_INVENTORY_cs_engines C35 (Merchant/Boothroyd,\ncad-to-simulation-I cnc_cutting_process.py @4fd03eb) but the Constitution is replaced by measured bentorque (Ganeyev 2025).\nThe heat line core is verified against analytical Green function (se selftest()).\n"
import numpy as np
import numba as nb
PARAMS = {'k_cort': 0.58, 'rho_cort': 1800.0, 'cp_cort': 1260.0, 'k_marrow': 0.3, 'rho_marrow': 1029.0, 'cp_marrow': 2666.0, 'bvtv': 0.31, 'U0_J_mm3': 0.683, 'F0_mm_rev': 0.03, 'M_SIZE': 0.44, 'eta': 0.5, 'h_irr': 4000.0, 'h_air': 25.0, 'perf_cort_ml_min_100g': 3.71, 'perf_canc_ml_min_100g': 10.0, 'rho_blood': 1050.0, 'cp_blood': 3617.0, 'point_half_angle_deg': 59.0, 'chi_wall': 0.5}

def hs_bounds(k1, k2, f1):
    "Hashin–Shtrikman isotropic double-phase conductivity limits (k1 > k2, volymandel f1 av fas 1)."
    f2 = 1 - f1
    upper = k1 + f2 / (1 / (k2 - k1) + f1 / (3 * k1))
    lower = k2 + f1 / (1 / (k1 - k2) + f2 / (3 * k2))
    return (lower, upper)

def canc_props(p, bvtv=None, bound='upper'):
    f = p['bvtv'] if bvtv is None else bvtv
    (lo, up) = hs_bounds(p['k_cort'], p['k_marrow'], f)
    k = up if bound == 'upper' else lo if bound == 'lower' else 0.5 * (lo + up)
    rc = f * p['rho_cort'] * p['cp_cort'] + (1 - f) * p['rho_marrow'] * p['cp_marrow']
    return (k, rc)

def u_spec(f_rev_mm, p):
    "Specific cutting energy J/mm3 (= N/mm2) vid matning per varv f_rev [mm/varv]."
    return p['U0_J_mm3'] * (np.maximum(f_rev_mm, 0.0001) / p['F0_mm_rev']) ** (-p['M_SIZE'])

def feed_from_force(F_N, D_mm, rpm, p=None):
    """Matning [mm/s] ur axialkraft via Ganeyev 2025 MONO: F = 16,91 (D/4) (f_rev/0,03)^0,41 N (kortikal, hel borr)."""
    f_rev = 0.03 * (F_N / (16.91 * D_mm / 4.0)) ** (1 / 0.41)
    return f_rev * rpm / 60.0

@nb.njit(cache=True)
def _step(T, Tn, solid, kc, rcV, Gr, Gz, perf, Ta, Q, hwall, Tcool, T0, dt, Nr, Nz, dr, dz, htop, rface, rc):
    for i in range(Nr):
        for j in range(Nz):
            if not solid[i, j]:
                Tn[i, j] = T[i, j]
                continue
            s = Q[i, j] + perf[i, j] * (Ta - T[i, j])
            if i + 1 < Nr:
                if solid[i + 1, j]:
                    s += Gr[i, j] * (T[i + 1, j] - T[i, j])
                else:
                    A = rface[i + 1] * dz
                    G = 1.0 / (1.0 / (hwall * A) + 0.5 * dr / (kc[i, j] * A))
                    s += G * (Tcool - T[i, j])
            else:
                A = rface[Nr] * dz
                s += kc[i, j] * A / dr * (T0 - T[i, j])
            if i > 0:
                if solid[i - 1, j]:
                    s += Gr[i - 1, j] * (T[i - 1, j] - T[i, j])
                else:
                    A = rface[i] * dz
                    G = 1.0 / (1.0 / (hwall * A) + 0.5 * dr / (kc[i, j] * A))
                    s += G * (Tcool - T[i, j])
            A = rc[i] * dr
            if j + 1 < Nz:
                if solid[i, j + 1]:
                    s += Gz[i, j] * (T[i, j + 1] - T[i, j])
                else:
                    G = 1.0 / (1.0 / (hwall * A) + 0.5 * dz / (kc[i, j] * A))
                    s += G * (Tcool - T[i, j])
            else:
                s += kc[i, j] * A / dz * (T0 - T[i, j])
            if j > 0:
                if solid[i, j - 1]:
                    s += Gz[i, j - 1] * (T[i, j - 1] - T[i, j])
                else:
                    G = 1.0 / (1.0 / (hwall * A) + 0.5 * dz / (kc[i, j] * A))
                    s += G * (Tcool - T[i, j])
            else:
                G = 1.0 / (1.0 / (htop * A) + 0.5 * dz / (kc[i, j] * A))
                s += G * (Tcool - T[i, j])
            Tn[i, j] = T[i, j] + dt * s / rcV[i, j]

@nb.njit(cache=True)
def _cem(T, solid, cem, tabove, Tmax, dt, Nr, Nz):
    for i in range(Nr):
        for j in range(Nz):
            if solid[i, j]:
                t = T[i, j]
                if t > Tmax[i, j]:
                    Tmax[i, j] = t
                R = 0.5 if t >= 43.0 else 0.25
                cem[i, j] += dt * R ** (43.0 - t) / 60.0
                if t >= 47.0:
                    tabove[i, j] += dt

@nb.njit(cache=True)
def _run(nsteps, T, Tn, solid, kc, rcV, Gr, Gz, perf, Ta, Q, hwall, Tcool, T0, dt, Nr, Nz, dr, dz, htop, rface, rc, cem, tabove, Tmax):
    """nsteps explicita steg + CEM43/Tmax/tid>47 °C-uppdatering; returnerar aktuell array (T eller Tn)."""
    ln05 = np.log(0.5)
    ln025 = np.log(0.25)
    for _ in range(nsteps):
        _step(T, Tn, solid, kc, rcV, Gr, Gz, perf, Ta, Q, hwall, Tcool, T0, dt, Nr, Nz, dr, dz, htop, rface, rc)
        tmp = T
        T = Tn
        Tn = tmp
        for i in range(Nr):
            for j in range(Nz):
                if solid[i, j]:
                    t = T[i, j]
                    if t > Tmax[i, j]:
                        Tmax[i, j] = t
                    lr = ln05 if t >= 43.0 else ln025
                    cem[i, j] += dt * np.exp(lr * (43.0 - t)) / 60.0
                    if t >= 47.0:
                        tabove[i, j] += dt
    return (T, Tn)

def simulate(steps, t_cort_mm=1.5, T0=37.0, perfusion=True, h_mm=0.1, p=None, probes=(), cool_after_s=60.0, R_extra_mm=10.0, Z_extra_mm=8.0, canc_bound='upper', bvtv=None, eta=None, t_cort_bottom_mm=0.0, record_every_s=0.05, T_art=None, chi=None, wall_irrig_during=True):
    "steps: list of dict(D, d_prev, depth, rpm, feed, irrig(bool), Tcool, pause) [mm, mm, mm, rpm, mm/s, -, °C, s].\n    depth = the depth of the tip under the crown (including the length of the tip). probes: (name, r_mm , z_mm ).\n    Returns dict with probe stories, T_max -, CEM43 fields, etc. (fields on grid , mm )."
    p = dict(PARAMS if p is None else p)
    if eta is not None:
        p['eta'] = eta
    if chi is not None:
        p['chi_wall'] = chi
    Dmax = max((s['D'] for s in steps))
    Hmax = max((s['depth'] for s in steps))
    R = Dmax / 2 + R_extra_mm
    Z = Hmax + Z_extra_mm
    dr = dz = h_mm * 0.001
    Nr = int(round(R / h_mm))
    Nz = int(round(Z / h_mm))
    rc = (np.arange(Nr) + 0.5) * dr
    zc = (np.arange(Nz) + 0.5) * dz
    rface = np.arange(Nr + 1) * dr
    (RR, ZZ) = np.meshgrid(rc, zc, indexing='ij')
    (kc_c, rc_c) = (p['k_cort'], p['rho_cort'] * p['cp_cort'])
    (kc_s, rc_s) = canc_props(p, bvtv, canc_bound)
    cort = ZZ < t_cort_mm * 0.001
    if t_cort_bottom_mm > 0:
        cort |= ZZ > (Z - t_cort_bottom_mm) * 0.001
    kc = np.where(cort, kc_c, kc_s)
    rcv = np.where(cort, rc_c, rc_s)
    V = rc[:, None] * dr * dz * np.ones((1, Nz))
    rcV = rcv * V
    Gr = np.zeros((Nr, Nz))
    Gz = np.zeros((Nr, Nz))
    kh = 2 * kc[:-1] * kc[1:] / (kc[:-1] + kc[1:])
    Gr[:-1] = kh * (rface[1:-1, None] * dz) / dr
    kz = 2 * kc[:, :-1] * kc[:, 1:] / (kc[:, :-1] + kc[:, 1:])
    Gz[:, :-1] = kz * (rc[:, None] * dr) / dz
    wb = np.where(cort, p['perf_cort_ml_min_100g'] * p['rho_cort'], p['perf_canc_ml_min_100g'] * (bvtv or p['bvtv']) * p['rho_cort'] + p['perf_canc_ml_min_100g'] * (1 - (bvtv or p['bvtv'])) * p['rho_marrow']) / 60.0 / 100.0 / 1000.0
    perf = wb * p['rho_blood'] * p['cp_blood'] * V if perfusion else np.zeros((Nr, Nz))
    Ta = T0 if T_art is None else T_art
    phi_mat = np.where(cort, 1.0, bvtv or p['bvtv'])
    Gb_r = kc * (rface[1:, None] * dz) / dr
    Gb_z = kc * (rc[:, None] * dr) / dz
    Gsum = 2 * Gb_r + 2 * Gb_z + perf
    dt = 0.95 * float(np.min(rcV / Gsum))
    T = np.full((Nr, Nz), float(T0))
    Tn = T.copy()
    solid = np.ones((Nr, Nz), bool)
    cem = np.zeros((Nr, Nz))
    tab = np.zeros((Nr, Nz))
    Tmax = T.copy()
    cotb = 1.0 / np.tan(np.radians(p['point_half_angle_deg']))
    pr_idx = [(n, min(max(int(round(r / h_mm - 0.5)), 0), Nr - 1), min(max(int(round(z / h_mm - 0.5)), 0), Nz - 1)) for (n, r, z) in probes]
    hist = {n: [] for (n, _, _) in pr_idx}
    th = []
    t = 0.0
    nextrec = 0.0
    energy_in = 0.0
    step_log = []

    def advance(duration, Q, hwall, htop, Tcool):
        nonlocal t, nextrec, T, Tn
        n = max(1, int(np.ceil(duration / dt)))
        d = duration / n
        done = 0
        while done < n:
            k = n - done
            if nextrec > t:
                k = min(k, max(1, int(np.ceil((nextrec - t) / d))))
            (T, Tn) = _run(k, T, Tn, solid, kc, rcV, Gr, Gz, perf, Ta, Q, hwall, Tcool, T0, d, Nr, Nz, dr, dz, htop, rface, rc, cem, tab, Tmax)
            done += k
            t += k * d
            if t >= nextrec - 1e-12:
                th.append(t)
                for (nm, i, j) in pr_idx:
                    hist[nm].append(T[i, j] if solid[i, j] else np.nan)
                nextrec += record_every_s
    Q0 = np.zeros((Nr, Nz))
    for (k, s) in enumerate(steps):
        irr = bool(s.get('irrig', True))
        Tcool = float(s.get('Tcool', T0))
        hwall = p['h_irr'] if irr else p['h_air']
        htop = hwall
        hwall_d = hwall if wall_irrig_during else p['h_air']
        if s.get('pause', 0) > 0:
            advance(s['pause'], Q0, hwall, htop, Tcool)
        t_start = t
        v = s['feed']
        f_rev = v / (s['rpm'] / 60.0)
        u = u_spec(f_rev, p)
        (Dm, dm, H) = (s['D'], s.get('d_prev', 0.0), s['depth'])
        t_drill = H / v
        nsub = max(20, int(np.ceil(t_drill / max(dt, 0.02))))
        dts = t_drill / nsub
        E_step = 0.0
        for q in range(nsub):
            ztip = (q + 1) * dts * v
            zcut = ztip - RR * 1000.0 * cotb
            inside = (RR * 1000.0 < Dm / 2) & (ZZ * 1000.0 < zcut)
            newly = inside & solid
            if newly.any():
                solid[newly] = False
            band = solid & (RR * 1000.0 < Dm / 2 + 0.5 * h_mm) & (RR * 1000.0 >= dm / 2 - 0.5 * h_mm) & (ZZ * 1000.0 >= zcut) & (ZZ * 1000.0 < zcut + h_mm)
            if not band.any():
                advance(dts, Q0, hwall_d, htop, Tcool)
                continue
            Acut = np.pi / 4 * (Dm ** 2 - (dm ** 2 if ztip <= s.get('prev_depth', H) + 1e-09 else 0.0))
            phi = float(np.sum(phi_mat[band] * V[band]) / np.sum(V[band]))
            P = p['eta'] * u * v * Acut * phi
            chi = p['chi_wall']
            wb_ = V * band
            l_eng = max(ztip - Dm / 2 * cotb, 0.0)
            wall = solid & (RR * 1000.0 >= Dm / 2) & (RR * 1000.0 < Dm / 2 + h_mm) & (ZZ * 1000.0 < ztip - Dm / 2 * cotb)
            ww = V * wall
            Q = (1 - chi) * P / (2 * np.pi) * wb_ / wb_.sum()
            Pw = chi * P * l_eng / (0.5 * H)
            if ww.sum() > 0 and Pw > 0:
                Q = Q + Pw / (2 * np.pi) * ww / ww.sum()
            P = (1 - chi) * P + (Pw if ww.sum() > 0 else 0.0)
            advance(dts, Q, hwall_d, htop, Tcool)
            E_step += P * dts
        energy_in += E_step
        step_log.append({'step': k, 'D': Dm, 'd_prev': dm, 'depth': H, 'rpm': s['rpm'], 'feed': v, 'f_rev_mm': f_rev, 'u_J_mm3': float(u), 't_drill_s': t_drill, 'E_bone_J': E_step, 'P_mean_W': E_step / t_drill, 't_start_s': t_start, 't_end_drill_s': t, 'removed_vol_mm3': float(np.sum(~solid * V) * 2 * np.pi * 1000000000.0)})
    advance(cool_after_s, Q0, p['h_irr'] if steps[-1].get('irrig', True) else p['h_air'], p['h_irr'] if steps[-1].get('irrig', True) else p['h_air'], float(steps[-1].get('Tcool', T0)))
    return {'t': np.array(th), 'probes': {k: np.array(v) for (k, v) in hist.items()}, 'Tmax': Tmax, 'cem43_min': cem, 't_above47_s': tab, 'solid': solid, 'r_mm': rc * 1000.0, 'z_mm': zc * 1000.0, 'dt_s': dt, 'h_mm': h_mm, 'E_bone_J': energy_in, 'steps': step_log, 'T0': T0}

def rosenthal_jaeger(steps, r_mm, z_mm, T0=20.0, p=None, k=None, rc=None, cool_after_s=60.0, n_t=4000):
    "BASLINJE : transient variable point source in infinitely homogeneous medium (Jaeger/Green function), same heating power P(t)\n    as the numerical model but without holes, layer, cooling or perfusion. The source of the shoulder at the depth of the tip.\n    Returns (t, T) at the point (r_mm, z_mm)."
    p = dict(PARAMS if p is None else p)
    k = p['k_cort'] if k is None else k
    rc = p['rho_cort'] * p['cp_cort'] if rc is None else rc
    alpha = k / rc
    (src_t, src_z, src_P) = ([], [], [])
    t0 = 0.0
    for s in steps:
        t0 += s.get('pause', 0.0)
        v = s['feed']
        f_rev = v / (s['rpm'] / 60.0)
        u = u_spec(f_rev, p)
        Acut = np.pi / 4 * (s['D'] ** 2 - s.get('d_prev', 0.0) ** 2)
        td = s['depth'] / v
        tt = np.linspace(0, td, 400, endpoint=False) + td / 800
        src_t += list(t0 + tt)
        src_z += list(v * tt)
        src_P += [p['eta'] * u * v * Acut] * len(tt)
        t0 += td
        src_dt = td / 400
        src_P[-len(tt):] = [x * src_dt for x in src_P[-len(tt):]]
    src_t = np.array(src_t)
    src_z = np.array(src_z) * 0.001
    src_E = np.array(src_P)
    T_end = t0 + cool_after_s
    tq = np.linspace(0, T_end, n_t)
    out = np.full(n_t, T0)
    rr = r_mm * 0.001
    zz = z_mm * 0.001
    for (i, t) in enumerate(tq):
        m = src_t < t
        if not m.any():
            continue
        tau = t - src_t[m]
        d2 = rr ** 2 + (zz - src_z[m]) ** 2
        out[i] += np.sum(src_E[m] / (rc * (4 * np.pi * alpha * tau) ** 1.5) * np.exp(-d2 / (4 * alpha * tau)))
    return (tq, out)

def selftest():
    "Verification of the articulated core: instantaneous ring source   holes .\n    Compare numerical solution (no removal, source in a cell on the axis, constant effect) with analytical\n    continuous point source dT = P/( 4 pi k R) erfc(R/( 2 sqrt(alpha t)))."
    from scipy.special import erfc
    p = dict(PARAMS)
    p['perf_cort_ml_min_100g'] = 0.0
    h = 0.1
    R = 10.0
    Z = 20.0
    (Nr, Nz) = (int(R / h), int(Z / h))
    dr = dz = h * 0.001
    rc = (np.arange(Nr) + 0.5) * dr
    rface = np.arange(Nr + 1) * dr
    k = p['k_cort']
    rcv = p['rho_cort'] * p['cp_cort']
    alpha = k / rcv
    kc = np.full((Nr, Nz), k)
    V = rc[:, None] * dr * dz * np.ones((1, Nz))
    rcV = rcv * V
    Gr = np.zeros((Nr, Nz))
    Gz = np.zeros((Nr, Nz))
    Gr[:-1] = k * (rface[1:-1, None] * dz) / dr
    Gz[:, :-1] = k * (rc[:, None] * dr) / dz
    solid = np.ones((Nr, Nz), bool)
    P = 1.0
    j0 = Nz // 2
    Q = np.zeros((Nr, Nz))
    Q[0, j0] = P / (2 * np.pi)
    T = np.zeros((Nr, Nz))
    Tn = T.copy()
    Gs = 2 * Gr + 2 * Gz + 4 * k * rc[:, None] * dr / dz
    dt = 0.3 * float(np.min(rcV / (Gs + 1e-30)))
    tend = 5.0
    n = int(tend / dt)
    for _ in range(n):
        _step(T, Tn, solid, kc, rcV, Gr, Gz, np.zeros((Nr, Nz)), 0.0, Q, 1e-09, 0.0, 0.0, dt, Nr, Nz, dr, dz, 1e-09, rface, rc)
        (T, Tn) = (Tn, T)
    t = n * dt
    res = []
    for dist_mm in (1.0, 2.0, 3.0):
        j = j0 + int(round(dist_mm / h))
        Rm = Nz and abs((j - j0) * dz)
        ana = P / (4 * np.pi * k * Rm) * erfc(Rm / (2 * np.sqrt(alpha * t)))
        res.append({'dist_mm': dist_mm, 'num': float(T[0, j]), 'ana': float(ana), 'rel_err': float(T[0, j] / ana - 1)})
    return {'t_s': t, 'h_mm': h, 'points': res}
if __name__ == '__main__':
    import json
    print(json.dumps(selftest(), indent=1))
