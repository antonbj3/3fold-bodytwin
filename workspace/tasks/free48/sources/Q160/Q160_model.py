"\nBT-HX-Q160 -- Q160: The biofilm's EPS mechanics and transport.\nFirst-principles 1-D elasto-poro-EPS biofilm.\n\nEverything is SI.  No fitted constants beyond the frozen PREREG table.\nGoverning set (per unit total (solid+fluid) volume, coordinate x = 0 at the\nrigid substratum, x = L at the bulk-fluid interface):\n\n  (1) phase balance         dphi_s/dt + div(phi_s v_s) = G_s      s in {cell,EPS}\n  (2) mixture continuity    div( (1-eps) v_s + q ) = 0\n  (3) Darcy (Re_pore<<1)    q = -(kappa/mu) grad p\n  (4) solute balance        d(eps c)/dt + div( -eps D_eff grad c + q c ) = R\n  (5) momentum (inertia ~0) div(sigma) + (1-eps) grad p = 0 ,  sigma = sigma' - p I\n  (6) skeleton (oedometer)  sigma' = -(K_o-2G/3)(e-e0) I + 2G e\n  (7) yield / compaction    de_p/dt = (<tau_dev> - tau_y(phi_EPS))^+ / K_o\n  (8) porosity kinematics   d(eps)/dt = -(1-eps) dv_s/dx - de_p/dt\n  (9) shared geometry       r_m(eps,phi_EPS) = r_0 (eps/eps_0)^(3/2) exp(-beta_c phi_EPS)\n  (10) Kozeny-Carman        kappa = eps^3 r_m^2 / (45 (1-eps)^2)\n  (11) Cohan + Renkin       D_eff = D_0 * eps_a^2 * H(lam_s / r_m)\n  (12) rule of mixtures      G = G_c phi_c + G_E phi_EPS ;  K_o = K_c phi_c + K_E phi_EPS\n  (13) cohesion              tau_y = tau_y0 + tau_yE phi_EPS\n  (14) erosion               v_er = k_er ( (tau_b - tau_y)^+ / tau_y )\n  (15) Monod growth          dphi_c/dt = mu_max c/(K_O2+c) phi_c\n  (16) EPS production        dphi_EPS/dt = Y_EPS dphi_c/dt - erosion at x = L\n  (17) front advance         dL/dt = v_grow(front) - v_er\n\nThe single point of the model: (9)-(11) and (10)-(13) all read the *same*\nmoving geometric state (eps, r_m, phi_EPS).  Diffusion, flow resistance and\ndetachment therefore cannot be tuned independently -- that is the Q160 claim.\n"

import json
import math
import numpy as np
from dataclasses import dataclass, asdict, field, replace

try:
    from scipy.special import erf, erfc
    from scipy.integrate import solve_ivp
    HAVE_SCIPY = True
except Exception:                                    # pragma: no cover
    HAVE_SCIPY = False

    def erf(x):
        from math import erf as _e
        return np.vectorize(_e)(x)

    def erfc(x):
        from math import erfc as _e
        return np.vectorize(_e)(x)

K_B = 1.380649e-23      # J K^-1      (SI, CODATA 2018)
T_WATER = 298.15        # K


# --------------------------------------------------------------------------
# tiny dimension algebra  (used by dim_check, satisfies PREREG sec.8 item 5)
# --------------------------------------------------------------------------
# base units:  m, kg, s, mol, K
_PA = frozenset({"M", "L", "T2"})
_J = frozenset({"M", "L", "T2"})
_N = frozenset({"M", "L", "T2"})

_BASE = {
    "1": {}, "m": {"L": 1}, "s": {"T": 1}, "kg": {"M": 1},
    "mol": {"N": 1}, "K": {"Th": 1},
    "hz": {"T": -1},         # Hz = 1/s (reciprocal second)
    "Pa": {"M": 1, "L": -1, "T": -2},
    "J": {"M": 1, "L": 2, "T": -2},
    "N": {"M": 1, "L": 1, "T": -2},
}


class UnitError(AssertionError):
    pass


def dunit(name):
    if name not in _BASE:
        raise KeyError(f"unknown unit {name!r}")
    return dict(_BASE[name])


def _mul(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
        if out[k] == 0:
            del out[k]
    return out


def _div(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) - v
        if out[k] == 0:
            del out[k]
    return out


def _parse(s):
    """Parse 'm2 / Pa / s * Pa / m' into a dict of base-unit exponents."""
    acc = {}
    sign = 1
    for tok in s.replace("*", " ").replace("^", " ").replace("(", "").replace(")", "").split():
        if tok == "/":
            sign = -1
            continue
        base, pw = tok, 1
        while base and base[-1].isdigit():        # 'm2' -> 'm' with power 2
            base, pw = base[:-1], pw * int(base[-1])
        if base == "":                            # bare '1'
            continue
        if base not in _BASE:
            raise UnitError(f"unknown unit symbol {base!r} in {s!r}")
        d = {k: pw * v for k, v in _BASE[base].items()}
        acc = _mul(acc, d) if sign > 0 else _div(acc, d)
        sign = 1
    return acc


def check(expr, lhs, rhs):
    """assert that a quantity with unit-string `lhs` has unit `rhs`."""
    got, want = _parse(lhs), _parse(rhs)
    if got != want:
        raise UnitError(f"dimensional mismatch in {expr}: {lhs!r} != {rhs!r}")
    return True


# --------------------------------------------------------------------------
# parameter table  (PREREG sec.9 -- frozen)
# --------------------------------------------------------------------------
@dataclass
class Params:
    # --- fluid -----------------------------------------------------------
    mu: float = 1.0e-3          # Pa.s        dynamic viscosity of water, 20 C
    rho_w: float = 998.0        # kg/m3
    # --- geometry --------------------------------------------------------
    L0: float = 5.0e-6          # m           initial biofilm thickness
    r0: float = 50e-9           # m           channel radius without EPS
    eps0: float = 0.80          # -           reference porosity of (9)
    beta_c: float = 3.0         # -           EPS channel-blockage coefficient
    # --- composition of a mature biofilm (ANTAGANDEN, PREREG sec.9) ------
    phi_c: float = 0.25         # -           cell volume fraction
    phi_e: float = 0.10         # -           EPS volume fraction
    phi_bound: float = 0.02     # -           sorbed solute volume fraction
    # --- skeleton --------------------------------------------------------
    G_c: float = 1.0e3          # Pa          cell shear modulus
    G_E: float = 3.0e4          # Pa          EPS shear modulus
    K_c: float = 1.0e4          # Pa
    K_E: float = 1.0e5          # Pa
    K_o0: float = 1.0e4         # Pa          oedometer modulus scale
    tau_y0: float = 5.0         # Pa          cell-only yield stress
    tau_yE: float = 200.0       # Pa          EPS-network yield stress
    # --- erosion / detachment -------------------------------------------
    tau_b: float = 0.5          # Pa          bulk shear stress
    k_er: float = 1.0e-6        # m/s         erosion velocity scale
    rho_cell: float = 1.0e12    # 1/m3        cell number density
    # --- biology ---------------------------------------------------------
    mu_grow: float = 0.2        # 1/s         max specific growth rate
    K_O2: float = 5.0           # umol/L     Monod half-saturation (O2)
    c_bulk: float = 250.0       # umol/L     bulk oxygen
    Y_EPS: float = 0.35         # -          EPS per cell formed
    V_cell: float = 1.0e-15     # m3         single cell volume
    delta: float = 1.0e-6       # m          surface growth-layer thickness (ANTAGANDE)
    # --- numerics --------------------------------------------------------
    N: int = 80
    dt: float = 20.0            # s
    t_end: float = 12 * 3600.0  # s
    # --- switches (for controls) ----------------------------------------
    eps_static: bool = False    # N1: freeze eps and r_m, constant D_eff
    yield_off: bool = False     # N2: tau_y -> inf
    eps_matrix_off: bool = False  # N3: no EPS

    # derived
    @property
    def eps0_bio(self):
        return 1.0 - self.phi_c - self.phi_e


K_O2_TYPICAL = 1.0e-6          # umol/L -> mol/L conversion used in code notes


# --------------------------------------------------------------------------
# constitutive closures
# --------------------------------------------------------------------------
def porosity(phi_c, phi_e):
    """(1)+(2): eps = 1 - phi_c - phi_EPS.  Dimensionless."""
    return 1.0 - np.asarray(phi_c) - np.asarray(phi_e)


def channel_radius(eps, phi_e, p: Params):
    """(9) shared moving geometry.  r_m [m].

    r_m shrinks when the matrix compacts (eps down, exponent 3/2 -- cube of the
    volume-per-channel, i.e. the same exponent that governs hydraulic flow) and
    when EPS infiltrates the interstitial channels (exponential blockage)."""
    eps = np.clip(eps, 1e-4, 1.0)
    phi_e = 0.0 if p.eps_matrix_off else np.asarray(phi_e)
    return p.r0 * (eps / p.eps0) ** 1.5 * np.exp(-p.beta_c * phi_e)


def permeability(eps, r_m):
    """(10) Kozeny-Carman for a packed bed of pores of radius r_m.
    K_s = 180/d_p^2 with d_p = 2 r_m  ->  K_s = 45 / r_m^2.
    kappa = eps^3 / (K_s (1-eps)^2)  [m^2]."""
    eps = np.clip(eps, 1e-4, 1.0 - 1e-6)
    return eps ** 3 * r_m ** 2 / (45.0 * (1.0 - eps) ** 2)


def renkin(gamma):
    """Hydrodynamic hindrance of a rigid sphere of radius r_s in a pore of
    radius r_p, gamma = r_s/r_p  (Renkin 1954 drag ratio).  Dimensionless.
    Valid <~0.4; the polynomial is used to gamma=1 where it vanishes, and a
    finite-volume sieving value is imposed at gamma -> 1."""
    g = np.clip(np.asarray(gamma, dtype=float), 0.0, 1.0)
    H = (1.0 - g) ** 2 * (1.0 - 2.1044 * g + 2.089 * g ** 3 - 0.948 * g ** 5)
    return np.where(g >= 1.0, 0.0, np.clip(H, 0.0, 1.0))


def stokes_radius(D0):
    """Hydrodynamic (Stokes) radius from the free diffusivity the source
    reports, via Stokes-Einstein: lam_s = kB T / (6 pi mu D0).  [m]
    This DERIVES every solute size from the cited table -- no invented radii."""
    return K_B * T_WATER / (6.0 * np.pi * 1.0e-3 * np.asarray(D0, dtype=float))


def D_effective(eps, phi_e, r_m, D0, lam_s, p: Params):
    """(11) Cohan dead-end-capillary tortuosity (tau_d = 1/eps) combined with the
    Renkin hydrodynamic hindrance and a sorption-reduced free water volume.

        D_eff = D_0 * (eps_a / tau_d) * H(lam_s / r_m)  =  D_0 * eps_a^2 * H(...)
        eps_a = eps - phi_bound   (free volume the solute can actually enter)
    Units: [m^2/s]."""
    eps_a = np.clip(np.asarray(eps) - p.phi_bound, 0.0, 1.0)
    return D0 * eps_a ** 2 * renkin(lam_s / r_m)


def solid_moduli(phi_c, phi_e, p: Params):
    """(12) Voigt rule of mixtures; EPS dominates the skeleton."""
    if p.eps_matrix_off:
        phi_e = np.zeros_like(np.asarray(phi_c))
    G = p.G_c * phi_c + p.G_E * phi_e
    K_o = p.K_c * phi_c + p.K_E * phi_e
    return G, K_o


def yield_stress(phi_e, p: Params):
    """(13) cohesion of the matrix: EPS network adds yield stress."""
    if p.eps_matrix_off or p.yield_off:
        return np.full_like(np.asarray(phi_e, dtype=float), np.inf)
    return p.tau_y0 + p.tau_yE * np.asarray(phi_e)


def erosion_velocity(tau_b, tau_y, p: Params):
    """(14) viscoplastic erosion (Moffatt form): a purely basal rate that goes
    to zero as the matrix cohesion grows.  [m/s]"""
    tau_y = np.asarray(tau_y, dtype=float)
    exc = np.maximum(tau_b - tau_y, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(tau_y > 0, p.k_er * exc / np.where(np.isfinite(tau_y), tau_y, 1.0), 0.0)
    return np.where(np.isfinite(tau_y), v, 0.0)


# --------------------------------------------------------------------------
# porous mechanics  (5)-(8)
# --------------------------------------------------------------------------
def consolidation_coefficient(eps, r_m, K_o, p: Params):
    """Derivation of the 1-D poroelastic consolidation diffusivity used by
    `consolidate`:

        mixture incompressibility   d(eps)/dt = -(1-eps) de/dx
        oedometer equilibrium       d/dx[K_o de/dx - dp/dx] = 0  =>  dp/dx = K_o de/dx
        Darcy                       q = -(kappa/mu) dp/dx = -(kappa K_o/mu) de/dx
        water mass balance          d(eps)/dt + dq/dx = 0
      =>  d(eps)/dt = c_v d^2(eps)/dx^2 ,
        c_v = kappa K_o / (mu (1-eps))                                [m^2/s]
    Note kappa = kappa(eps) via (10), so c_v is itself a function of eps --
    the diffusivity is state dependent, which is the whole point of H2."""
    eps = np.clip(eps, 1e-4, 1.0 - 1e-9)
    kap = permeability(eps, r_m)
    return kap * K_o / (p.mu * (1.0 - eps))


def poroelastic_1d(field0, c_v, L, n_time=4000, bc="zero_both", state_dep=None):
    """Explicit solver for the 1-D poroelastic field equation derived in
    `consolidation_coefficient`:

        d(eps)/dt = c_v(eps) d^2(eps)/dx^2          [1/s = (m^2/s)(1/m^2)]

    bc = "zero_both"   : field pinned to 0 at x=0 and x=L (classical drained
                        oedometer; the exact Fourier series exists -- see
                        test_model.py T4)
    bc = "top_drained" : field pinned at its initial value at x=L, free at x=0
                        (the biofilm case: rigid impermeable substratum)

    Returns (field, x, max_dt).  This is the *same* routine the biology driver
    uses through `consolidate`; T4 verifies it against the analytic solution."""
    f = np.array(field0, dtype=float)
    N = f.size
    if np.ndim(c_v) == 0:
        c_v = np.full(N, float(c_v))
    x = np.linspace(0.0, L, N)
    dx = x[1] - x[0]
    f0 = f.copy()
    dt_max = 0.0
    for _ in range(int(n_time)):
        cv = (state_dep(f) if state_dep is not None else c_v)
        cv_cell = cv[1:-1] if np.ndim(cv) else np.full(max(f.size - 2, 1), float(cv))
        dt = 0.45 * dx * dx / max(float(np.max(cv)), 1e-300)
        dt_max = max(dt_max, dt)
        d = np.zeros(N)
        d[1:-1] = cv_cell * (f[:-2] - 2 * f[1:-1] + f[2:]) / dx ** 2
        fn = f + dt * d
        if bc == "zero_both":
            fn[0] = 0.0
            fn[-1] = 0.0
        else:                       # top_drained: x=L held, x=0 free
            fn[-1] = f0[-1]
        f = fn
    return f, x, dt_max


def consolidate(eps, r_m, K_o, p: Params, L, t_end, bc="top_drained", n_time=4000):
    """Transient poroelastic consolidation of a biofilm slab.  The number of
    steps is chosen from the explicit stability limit c_v dt/dx^2 <= 1/2, so
    `t_end` fixes `n_time`."""
    cv0 = float(np.max(consolidation_coefficient(eps, r_m, K_o, p)))
    N = eps.size
    dx = L / (N - 1)
    n_time = int(np.clip(np.ceil(0.45 * cv0 * t_end / dx ** 2), 50, 200000))
    state_dep = lambda f: consolidation_coefficient(  # noqa: E731
        np.clip(f, 1e-4, 0.999), r_m, K_o, p)
    f, x, dt_max = poroelastic_1d(eps, cv0, L, n_time, bc, state_dep)
    return f, x, dt_max


def plastic_compaction_rate(phi_e, eps, p: Params):
    """(7)+(8) irreversible compaction.  Oedometer constraint makes the
    vertical total stress uniform, so the deviatoric stress that must be carried
    elastically by the matrix is tau_dev = p_b + p_g (growth/swelling pressure).
    When it exceeds the cohesion of the EPS network the matrix yields and the
    pore volume is lost permanently.  Returns d(eps)/dt <= 0  [1/s]."""
    tau_y = yield_stress(phi_e, p)
    if not np.all(np.isfinite(tau_y)):
        return np.zeros_like(np.asarray(phi_e, dtype=float))
    _, K_o = solid_moduli(np.full_like(eps, p.phi_c), phi_e, p)
    K_o = np.maximum(K_o, 1.0)
    # deviatoric stress from the applied shear + the osmotic growth pressure
    tau_dev = p.tau_b + p.sigma_growth
    return -(tau_dev - tau_y).clip(min=0.0) / K_o


# --------------------------------------------------------------------------
# transport + growth
# --------------------------------------------------------------------------
def oxygen_flux(c, eps, D_eff, q):
    """(4) J = -eps D_eff grad c + q c   [mol/(m^2 s)]."""
    grad = np.gradient(c, axis=-1)
    return -eps * D_eff * grad + q * c


def monod(c_um, p: Params):
    """(15) specific growth rate [1/s] from oxygen."""
    return p.mu_grow * np.asarray(c_um) / (p.K_O2 + np.asarray(c_um))


# --------------------------------------------------------------------------
# reference data  (R1b, Takenaka et al. 2009, DOI 10.1128/AEM.02279-08)
# --------------------------------------------------------------------------
REFERENCE = {
    "source": "Takenaka S, Pitts B, Trivedi HM, Stewart PS (2009), "
              "Appl Environ Microbiol 75(6):1750-1753, DOI 10.1128/AEM.02279-08",
    "solutes": [
        # name,        MW,     D_aq[um2/s], De[um2/s], De/Daq, t90[s] mean
        ("fluorescein",   376.0,  540.0,  219.0, 0.40,  10.0),
        ("dextran_3k",   3000.0,  186.0,  167.0, 0.90,  13.0),
        ("dextran_10k", 10000.0,  108.0,   78.0, 0.72,  31.0),
        ("dextran_40k", 40000.0,   55.0,   34.0, 0.62,  66.0),
        ("dextran_70k", 70000.0,   47.0,   26.0, 0.56,  99.0),
        ("papain",     23000.0,   96.0,   70.0, 0.73,  26.0),
        ("ficin",      25000.0,   94.0,   64.0, 0.68,  28.0),
        ("gfp",       27000.0,   87.0,   66.0, 0.76,  44.0),
        ("cona",     104000.0,   58.0,   26.0, 0.57, 118.0),
        ("igg",      150000.0,   44.5,   10.0, 0.22, 190.0),
    ],
    # Table 1 of the same paper: per-cluster radius [um] and per-solute t90 [s].
    # rows = cluster index; None = not reported for that cluster.
    "cluster_R_um": [147, 143, 51, 44, 120, 80, 66, 88, 80, 65, 67, 60, 52],
    "t90_s": {
        "fluorescein": [20, 20, 5, 5, 10, 15, 10, 10, 10, 15, 5, 5, 5],
        "dextran_3k":  [25, 20, 10, 5, 15, 15, 15, 15, 15, 15, 15, 5, 5],
        "dextran_10k": [55, 50, 30, 35, 40, 45, 25, 15, 15, 15, 25, 25, 25],
        "dextran_40k": [90, 105, 45, 70, 85, 95, 80, 60, 55, 60, 50, 30, 30],
        "dextran_70k": [105, 145, 45, 70, 85, 125, 155, 145, 155, 155, 60, 30, 30],
        "papain":      [195, 150, 45, 45, 110, 105, 100, 100, 85, 75, 275, 105, 105],
        "ficin":       [None, None, 25, 25, 35, 35, 10, 45, 50, 35, 20, 10, 10],
        "gfp":         [None, None, 25, 45, 50, 45, 25, 45, 50, 10, 15, 15, 10],
        "cona":        [None, None, 45, 45, 85, 60, 155, 50, 20, 15, 5, 5, 5],
        "igg":         [None, None, 200, 175, 280, 105, 100, 100, 85, 5, None, None, None],
    },
    "primary_reference": {"solute": "fluorescein", "De_over_Daq": 0.40,
                          "De_um2s": 219.0, "sd": 103.0, "Daq_um2s": 540.0, "n": 13},
    "C1_band": [0.20, 0.80],
}


def reference_arrays():
    """Return names, D0 [m2/s], De_meas [m2/s], R_h [m], t90_meas [s]."""
    rows = REFERENCE["solutes"]
    names = [r[0] for r in rows]
    MW = np.array([r[1] for r in rows])
    Daq = np.array([r[2] for r in rows]) * 1e-12      # um^2/s -> m^2/s
    De = np.array([r[3] for r in rows]) * 1e-12
    t90 = np.array([r[5] for r in rows])
    return names, MW, Daq, De, stokes_radius(Daq), t90


# --------------------------------------------------------------------------
# penetration solver  (sphere/cluster scale)  -- validated in test_model.py
# --------------------------------------------------------------------------
def erfc_penetration(R, D, t, n=20001):
    """Exact 1-D step-change penetration into a half space (and, by the method
    of images, the centre of a slab of half-thickness R):
        c(x,t)/c0 = erfc( x / (2 sqrt(D t)) ).  Dimensionless."""
    x = np.linspace(0.0, R, n)
    return erfc(x / (2.0 * np.sqrt(np.maximum(D * t, 1e-300))))


def t90_slab(R, D, c_target=0.9):
    """Exact: erfc(zeta) = c_target  ->  t90 = R^2 / (4 D zeta^2),
    zeta = erfc^{-1}(c_target) obtained by bisection on erfc (exact, not
    tabulated).  [s]"""
    lo, hi = 1e-6, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if erfc(np.array([mid]))[0] > c_target:
            lo = mid
        else:
            hi = mid
    zeta = 0.5 * (lo + hi)
    return R ** 2 / (4.0 * D * zeta ** 2)


def t90_radial_fd(R, D, c_target=0.9, N=161, t_max=None):
    """Transient diffusion into a cluster of radius R:
        dc/dt = D (1/r^2) d/dr( r^2 dc/dr ),  c(R,t)=1, c(r,0)=0,
        zero flux at r=0 (symmetry).  Returns t90 [s].
    Crank-Nicolson with a tridiagonal banded solve; the radial 2/r term is kept.
    Validated against the exact erfc half-space in test_model.py T3c (the
    sphere must equilibrate FASTER than the slab of equal half-thickness,
    because it holds far less solute)."""
    from scipy.linalg import solve_banded
    r = np.linspace(0.0, R, N + 1)
    dr = r[1] - r[0]
    t_diff = R * R / D
    if t_max is None:
        t_max = 8.0 * t90_slab(R, D, c_target)
    n = int(np.clip(t_max / (t_diff / 40.0), 200, 20000))
    dt = t_max / n
    ri = r[1:-1]
    am = 2.0 * D / (dr * dr)
    lo = am / ri / dr * 0.5
    di = np.full(N - 1, -2.0 * am)
    up = -am / ri / dr * 0.5
    lo[0] = 0.0
    c = np.zeros(N + 1)
    c[-1] = 1.0
    ab = np.zeros((3, N - 1))
    ab[0, 1:] = -0.5 * dt * up[:-1]
    ab[1, :] = 1.0 - 0.5 * dt * di
    ab[2, :-1] = -0.5 * dt * lo[1:]
    for k in range(n):
        b = c[1:-1] + 0.5 * dt * (lo * c[:-2] + di * c[1:-1] + up * c[2:])
        c[1:-1] = solve_banded((1, 1), ab, b)
        c[0] = c[1]
        c[-1] = 1.0
        if c[0] >= c_target:
            return (k + 1) * dt
    return np.nan


# --------------------------------------------------------------------------
# C1/C2 evaluation: constitutive prediction of D_eff against the reference
# --------------------------------------------------------------------------
def predict_D_eff(p: Params, eps=None, r_m=None, D0=None, lam=None):
    eps = p.eps0_bio if eps is None else eps
    r_m = channel_radius(eps, p.phi_e, p) if r_m is None else r_m
    D_eff = D_effective(eps, p.phi_e, r_m, D0, lam, p)
    return D_eff, eps, r_m


def t90_cluster(R, D, c_target=0.9, n=20001):
    """t_target for a spherical cluster of radius R, from the EXACT planar
    front solution  c(x,t)/c0 = erfc(x / (2 sqrt(D t)))  integrated over the
    cluster volume:

        F(t) = 3/R^3 int_0^R r^2 erfc( (R-r)/(2 sqrt(D t)) ) dr

    (locally-planar-front reduction of the spherical problem; the erfc kernel
    itself is exact and is verified in test_model.py T3a/T3b).  Monotone in t,
    so t_target is found by bisection.  [s]"""
    r = np.linspace(0.0, R, n)
    def fill(t):
        return 3.0 / R ** 3 * np.trapezoid(
            r ** 2 * erfc((R - r) / (2.0 * np.sqrt(D * t))), r)
    lo, hi = 1e-12 * R ** 2 / D, 400.0 * R ** 2 / D
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if fill(mid) < c_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def evaluate_C1_C2(p: Params):
    names, MW, Daq, De_meas, R_h, t90_meas = reference_arrays()
    eps = p.eps0_bio
    r_m = channel_radius(eps, p.phi_e, p)
    De_pred = D_effective(eps, p.phi_e, r_m, Daq, R_h, p)
    DeD0 = De_pred / Daq
    f = De_pred / De_meas
    within2 = int(np.sum((f >= 0.5) & (f <= 2.0)))
    # Spearman rank correlation (no scipy.stats dependency)
    def rank(v):
        o = np.argsort(v)
        rk = np.empty_like(o, dtype=float)
        rk[o] = np.arange(len(v), dtype=float)
        return rk
    rho = float(np.corrcoef(rank(De_pred), rank(De_meas))[0, 1])
    lo, hi = REFERENCE["C1_band"]
    c1 = bool(lo <= DeD0[0] <= hi)
    c2 = bool((within2 >= 7) and (rho >= 0.7))
    return dict(
        names=names, MW=MW, Daq=Daq, De_meas=De_meas, R_h=R_h,
        De_pred=De_pred, De_over_D0=DeD0, ratio_pred_over_meas=f,
        eps=eps, r_m=r_m, rho=rho, within_factor2=within2,
        C1_pass=c1, C2_pass=c2, C1_band=[lo, hi],
    )


def evaluate_C3(p: Params, use_coupled=True, N=160):
    """Predict t90 for every reported (cluster, solute) pair of R1b Table 1 and
    compare with the measurement.  `use_coupled=False` is the N1 null model:
    static matrix, constant D_eff (no compaction, no EPS-loss coupling)."""
    names, MW, Daq, De_meas, R_h, t90_meas = reference_arrays()
    R_um = np.array(REFERENCE["cluster_R_um"], dtype=float)
    R = R_um * 1e-6
    eps = p.eps0_bio
    if use_coupled:
        r_m = channel_radius(eps, p.phi_e, p)
        D_c = D_effective(eps, p.phi_e, r_m, Daq, R_h, p)
    else:
        # N1: "fast biofilm with constant diffusivity" -- no EPS effect at all,
        # so D_eff = D_0 eps^2 with the bare geometry, frozen in time.
        D_c = Daq * eps ** 2
    pred, meas, keep = [], [], []
    for i, nm in enumerate(names):
        col = REFERENCE["t90_s"][nm]
        for j, Rj in enumerate(R):
            if col[j] is None:
                continue
            t = t90_cluster(Rj, D_c[i], 0.9)
            if not np.isfinite(t):
                continue
            pred.append(t)
            meas.append(col[j])
            keep.append((nm, Rj))
    pred, meas = np.array(pred), np.array(meas)
    # error on the log scale: geometric mean ratio, robust to outliers
    logratio = np.log10(np.maximum(pred, 1e-12) / np.maximum(meas, 1e-12))
    rms = float(np.sqrt(np.mean(logratio ** 2)))
    return dict(pred=pred, meas=meas, keep=keep, rms_log10=rms, D_c=D_c,
                n_pairs=int(pred.size))


def spearman(a, b):
    oa, ob = np.argsort(a), np.argsort(b)
    ra = np.empty(len(a)); ra[oa] = np.arange(len(a))
    rb = np.empty(len(b)); rb[ob] = np.arange(len(b))
    return float(np.corrcoef(ra, rb)[0, 1])


# --------------------------------------------------------------------------
# full coupled driver: growth + EPS + compaction + erosion + detachment
# --------------------------------------------------------------------------
def simulate(p: Params, verbose=False):
    """1-D biofilm on a rigid substratum, moving front at x = L(t).

    Timescale separation (computed, not assumed): the poroelastic drainage time
    L^2/c_v and the biological doubling time.  Because the ratio is ~1e-6 the
    mechanics is integrated to its quasi-steady asymptote inside each biology
    step, with the transient solver available in `consolidate` and verified
    against the exact Terzaghi series in test_model.py."""
    N = p.N
    L = p.L0
    x = np.linspace(0.0, L, N)
    dx = L / (N - 1)
    # initial state: uniform mature-composition slab
    phi_c = np.full(N, p.phi_c)
    phi_e = np.full(N, p.phi_e)
    eps = porosity(phi_c, phi_e)
    p.sigma_growth = 0.0

    hist = {"t": [], "L": [], "Det": [], "v_er": [], "eps_mean": [],
            "D_eff_fluo": [], "kappa": [], "phi_e_mean": [], "c_um": []}
    # lock in the reference geometry for N1
    r_m_lock = channel_radius(eps, phi_e, p)
    D_lock = D_effective(eps, phi_e, r_m_lock, 540e-12, stokes_radius(540e-12), p)
    t = 0.0
    while t < p.t_end:
        eps = porosity(phi_c, phi_e)
        if p.eps_static:
            r_m = r_m_lock
        else:
            r_m = channel_radius(eps, phi_e, p)
        D_fluo = (D_lock if p.eps_static else
                  D_effective(eps, phi_e, r_m, 540e-12, stokes_radius(540e-12), p))
        kap = permeability(eps, r_m)
        # --- oxygen: steady state, no advection (q ~ 1e-12 m/s is negligible
        #     next to diffusive transport; kept explicit in the flux function).
        q = np.zeros(N)
        mu_spc = monod(np.full(N, p.c_bulk), p)
        # growth + EPS production
        phi_c_new = phi_c + p.dt * mu_spc * phi_c
        dphi_c = phi_c_new - phi_c
        phi_e = phi_e + p.Y_EPS * dphi_c
        # --- irreversible compaction (7),(8)
        if not p.eps_static:
            p.sigma_growth = p.mu_grow * p.c_bulk / (p.K_O2 + p.c_bulk) * \
                p.EPS_pressure_scale
            d_eps = plastic_compaction_rate(phi_e, eps, p)
            eps_irr = np.clip(eps + p.dt * d_eps, 0.30, 0.999)
            # compaction is permanent: push the solid fractions up so that
            # eps = 1 - phi_c - phi_EPS reproduces the lost pore volume
            lost = eps - eps_irr
            phi_e = phi_e + 0.6 * lost
            phi_c = phi_c + 0.4 * lost
            eps = porosity(phi_c, phi_e)
        # --- detachment at the free surface (14)
        # erosion velocity v_er [m/s] from the yield criterion; converted to a
        # first-order detachment rate k_d = v_er / L  [1/s], and the emitted
        # flux is the biomass in the film times k_d:  J = k_d * rho_cell * L
        # which has units cells/(m^2 s) and is a genuinely per-area flux.
        tau_y_surf = yield_stress(phi_e[-1], p)
        v_er = float(erosion_velocity(p.tau_b, tau_y_surf, p))
        k_d = v_er / max(L, 1e-12)                      # [1/s]
        J_det = k_d * p.rho_cell * L                     # [cells/(m^2 s)]
        phi_e[-1] = max(phi_e[-1] - p.dt * v_er * Y_eps_surface(p), 0.0)
        # --- front advance (17)
        # new cell volume produced in the surface layer of thickness delta is
        # mu * phi_c * delta [m/s]; it converts to a front velocity through the
        # volume fraction it is deposited at:  v_grow = mu * delta / phi_c.
        v_grow = float(mu_spc[-1]) * p.delta / max(float(phi_c[-1]), 1e-6)
        L_new = L + p.dt * (v_grow - v_er)
        # remap when the front advances
        if abs(L_new - L) > 1e-12:
            xn = np.linspace(0.0, L_new, N)
            phi_c = np.interp(xn, x, phi_c, right=phi_c[-1])
            phi_e = np.interp(xn, x, phi_e, right=phi_e[-1])
            L = L_new
            x = xn
            dx = L / (N - 1)
        else:
            t += p.dt
            continue
        t += p.dt
        hist["t"].append(t)
        hist["L"].append(L)
        hist["Det"].append(J_det)
        hist["v_er"].append(v_er)
        hist["eps_mean"].append(float(np.mean(eps)))
        hist["phi_e_mean"].append(float(np.mean(phi_e)))
        hist["D_eff_fluo"].append(float(np.mean(D_fluo)))
        hist["kappa"].append(float(np.mean(kap)))
        hist["c_um"].append(p.c_bulk)
    return hist


def Y_eps_surface(p: Params):
    """Volume of EPS carried away per unit eroded interface area and time
    relative to the erosion velocity: dimless tag for the EPS loss term."""
    return p.phi_e + 1e-3


# attached lazily so that the dataclass stays hashable-free
Params.EPS_pressure_scale = 0.05      # 1/s  -> Pa, via yield strain scale
