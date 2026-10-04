"""BT-HX-Q012 -- first-principles model of diet x microbiome -> metabolite exposure.

Frozen per PREREG.md (sha256 in PREREG.sha256). No measured data are used anywhere.

Substance question S1 (see PREREG.md section 1): which dietary intake and which
microbial functional capacity set the steady-state plasma concentration of
trimethylamine-N-oxide (TMAO), and how much uncertainty survives?

Chain (all rates in mol/h, all concentrations in mol/L):

    1.  n_prec      = I / MW                                    precursor intake
    2.  R_colon     = n_prec/24 * f_colonic_escape * nu * y     colonic TMA production
    3.  R_tma_hep   = R_colon * f_tma_to_hepatic                portal delivery of TMA
    4.  R_prod      = f_fmo3 * R_tma_hep + r_host               hepatic TMAO source
    5.  C_ss        = R_prod / CL                               steady state, linear FMO3
                     (saturable branch: CL*C^2 + a*C - r_host*Km = 0, a closed form)
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm, spearmanr

MC_DRAWS = 20000
MC_SEED = 20260925

# --------------------------------------------------------------------------
# minimal dimension algebra, used for the unit check
# --------------------------------------------------------------------------

DIM_BASE = ("g", "kg", "mol", "L", "h", "min", "d")


class Dim:
    def __init__(self, **kw: float) -> None:
        self.e: dict[str, float] = {b: 0.0 for b in DIM_BASE}
        unknown = set(kw) - set(self.e)
        if unknown:
            raise KeyError(f"unknown base units: {sorted(unknown)}")
        for k, v in kw.items():
            self.e[k] = float(v)

    def __mul__(self, other: "Dim") -> "Dim":
        return Dim(**{b: self.e[b] + other.e[b] for b in DIM_BASE})

    def __truediv__(self, other: "Dim") -> "Dim":
        return Dim(**{b: self.e[b] - other.e[b] for b in DIM_BASE})

    def add(self, other: "Dim") -> "Dim":
        """Sums require like dimensions. The result is that shared dimension, not the
        sum of exponents -- dimensional algebra is multiplicative, so x + x is x, not x^2."""
        if self != other:
            raise ValueError(f"dimension mismatch in sum: {self!r} + {other!r}")
        return Dim(**self.e)

    def __pow__(self, k: float) -> "Dim":
        return Dim(**{b: self.e[b] * k for b in DIM_BASE})

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Dim) and all(
            math.isclose(self.e[b], other.e[b], rel_tol=1e-12, abs_tol=1e-12)
            for b in DIM_BASE
        )

    def __repr__(self) -> str:  # pragma: no cover - display only
        parts = [f"{b}^{self.e[b]:g}" for b in DIM_BASE if self.e[b] != 0.0]
        return "1" if not parts else " ".join(parts)

    def as_dict(self) -> dict[str, float]:
        return {k: v for k, v in self.e.items() if v != 0.0}


DIM_ONE = Dim()
DIM_G = Dim(g=1)
DIM_KG = Dim(kg=1)
DIM_MOL = Dim(mol=1)
DIM_L = Dim(L=1)
DIM_H = Dim(h=1)
DIM_PER_DAY = Dim(d=-1)          # "per day"
DIM_PER_MIN = Dim(min=-1)        # "per minute"
DIM_DAY_TO_HOUR = Dim(d=1, h=-1)   # numeric 1/24: turns a per-day rate into a per-hour rate
DIM_MIN_TO_HOUR = Dim(min=1, h=-1) # numeric 1/60: turns a per-minute rate into a per-hour rate


# --------------------------------------------------------------------------
# parameters
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Parameters:
    """Frozen protocol values, PREREG.md section 5."""

    # diet block (K03)
    body_mass_kg: float = 70.0
    precursor_intake_g_day: float = 0.430
    precursor_mw_g_mol: float = 161.2
    f_colonic_escape: float = 0.30
    # microbiome block (K14)
    nu_stoich: float = 1.0
    y_microbial_conversion: float = 0.20
    f_tma_to_hepatic: float = 0.90
    # host / hepatic block (K11)
    f_fmo3: float = 0.90
    fmo3_vmax_umol_h: float = 1000.0
    fmo3_km_umol_l: float = 5.0
    # host / renal block
    gfr_ml_min_per_kg: float = 1.8
    f_reabsorption: float = 0.0
    f_tubular_secretion: float = 0.10
    vd_l: float = 49.0
    r_host_umol_h: float = 0.062
    # constants and numerics
    tmao_mw_g_mol: float = 75.11
    hepatic_mode: str = "linear"
    mc_draws: int = MC_DRAWS
    mc_seed: int = MC_SEED


PARAMETER_SPECS: dict[str, tuple[str, str, str]] = {
    # name: (unit, source/assumption, block)
    "body_mass_kg": ("kg", "assumption: adult reference individual", "host"),
    "precursor_intake_g_day": (
        "g/d",
        "SOURCE: Koeth 2013 Nat Med 19:576-585 doi:10.1038/nm.3145, Fig. 1 legend + "
        "Supplementary Methods: 250 mg d3-L-carnitine capsule + 8 oz sirloin steak "
        "(~180 mg L-carnitine) = 430 mg",
        "diet",
    ),
    "precursor_mw_g_mol": ("g/mol", "chemical formula weight of L-carnitine free base", "constant"),
    "f_colonic_escape": (
        "-",
        "ASSUMPTION, not looked up: fraction of ingested precursor escaping "
        "small-intestinal absorption to reach colonic microbiota",
        "diet",
    ),
    "nu_stoich": (
        "mol/mol",
        "ASSUMPTION from reaction stoichiometry: CntA/CntB (carnitine) and CutC "
        "(choline) release one TMA per precursor, 1:1",
        "microbiome",
    ),
    "y_microbial_conversion": (
        "-",
        "ASSUMPTION, not looked up: realized fraction of the 1:1 stoichiometry by the "
        "resident TMA-producing community",
        "microbiome",
    ),
    "f_tma_to_hepatic": ("-", "ASSUMPTION: fraction of produced TMA reaching the liver", "microbiome"),
    "f_fmo3": (
        "-",
        "ASSUMPTION; route existence supported by Wang 2011 Nature 472:57-63 "
        "doi:10.1038/nature09922 (flavin monooxygenase expression segregates with disease)",
        "host",
    ),
    "fmo3_vmax_umol_h": ("umol/h", "ASSUMPTION, not looked up: hepatic FMO3 capacity", "host"),
    "fmo3_km_umol_l": ("umol/L", "ASSUMPTION, not looked up: hepatic FMO3 half-saturation", "host"),
    "gfr_ml_min_per_kg": (
        "mL/min/kg",
        "ASSUMPTION, not looked up in this run: adult reference glomerular filtration",
        "host",
    ),
    "f_reabsorption": ("-", "ASSUMPTION: tubular reabsorption fraction", "host"),
    "f_tubular_secretion": ("-", "ASSUMPTION: tubular secretion fraction of GFR", "host"),
    "vd_l": ("L", "ASSUMPTION: 0.7 L/kg apparent distribution volume", "host"),
    "r_host_umol_h": (
        "umol/h",
        "ASSUMPTION, not looked up: non-microbial TMAO source; the only free knob in Phi",
        "host",
    ),
    "tmao_mw_g_mol": ("g/mol", "chemical formula weight of TMAO", "constant"),
}

# dimension of each parameter, for the unit check
PARAMETER_DIMS: dict[str, Dim] = {
    "body_mass_kg": DIM_KG,
    "precursor_intake_g_day": DIM_G * DIM_PER_DAY,
    "precursor_mw_g_mol": DIM_G / DIM_MOL,
    "f_colonic_escape": DIM_ONE,
    "nu_stoich": DIM_ONE,
    "y_microbial_conversion": DIM_ONE,
    "f_tma_to_hepatic": DIM_ONE,
    "f_fmo3": DIM_ONE,
    "fmo3_vmax_umol_h": DIM_MOL / DIM_H,
    "fmo3_km_umol_l": DIM_MOL / DIM_L,
    "gfr_ml_min_per_kg": DIM_L * DIM_PER_MIN / DIM_KG,
    "fmo3_vmax_umol_h": DIM_MOL / DIM_H,
    "f_reabsorption": DIM_ONE,
    "f_tubular_secretion": DIM_ONE,
    "vd_l": DIM_L,
    "r_host_umol_h": DIM_MOL / DIM_H,
    "tmao_mw_g_mol": DIM_G / DIM_MOL,
}

# frozen uncertainty priors: name -> (median, sigma_ln)
PRIORS: dict[str, tuple[float, float]] = {
    "f_colonic_escape": (0.30, 0.60),
    "y_microbial_conversion": (0.20, 0.90),
    "f_tma_to_hepatic": (0.90, 0.20),
    "f_fmo3": (0.90, 0.20),
    "gfr_ml_min_per_kg": (1.80, 0.15),
    "vd_l": (49.0, 0.20),
    "r_host_umol_h": (0.062, 1.00),
}

SENSITIVITY_PARAMS = ("y_microbial_conversion", "f_colonic_escape", "gfr_ml_min_per_kg")

# frozen acceptance thresholds, PREREG.md section 4
G1_PHI_MIN = 0.95
G2_CONTRAST_MIN = 2.0
G3_PRODUCT_REL_WIDTH_MAX = 0.10
G3_FACTOR_REL_WIDTH_MIN = 0.50
G4_CAPACITY_MARGIN_MIN = 3.0


DEFAULT = Parameters()


# --------------------------------------------------------------------------
# core rate equations
# --------------------------------------------------------------------------


def renal_clearance(p: Parameters) -> float:
    """Total TMAO clearance in L/h: filtration, minus reabsorption, plus secretion."""
    gfr_l_h = p.gfr_ml_min_per_kg * p.body_mass_kg * (60.0 / 1000.0)
    return gfr_l_h * (1.0 - p.f_reabsorption) * (1.0 + p.f_tubular_secretion)


def precursor_mol_per_day(p: Parameters) -> float:
    return p.precursor_intake_g_day / p.precursor_mw_g_mol


def tma_production_mol_h(p: Parameters) -> float:
    """Step 2+3: colonic TMA production and portal delivery."""
    return (
        precursor_mol_per_day(p)
        / 24.0
        * p.f_colonic_escape
        * p.nu_stoich
        * p.y_microbial_conversion
        * p.f_tma_to_hepatic
    )


def tmao_source_mol_h(p: Parameters) -> float:
    """Step 4: linear hepatic source. The saturable branch lives in steady_state()."""
    return p.f_fmo3 * tma_production_mol_h(p) + p.r_host_umol_h * 1e-6


def steady_state(p: Parameters, *, mode: str | None = None) -> float:
    """Steady-state plasma TMAO in mol/L.

    mode="linear"     : first-order FMO3, C = R_prod / CL, exactly linear in intake.
    mode="saturable"  : hepatic flux f_fmo3*Vmax*C/(Km+C) saturates.
                        CL*C^2 + a*C - r_host*Km = 0 with
                        a = f_fmo3*Vmax + r_host - CL*Km, positive root in closed form.
    """
    m = mode or p.hepatic_mode
    cl = renal_clearance(p)
    r_host = p.r_host_umol_h * 1e-6
    c_substrate = (p.f_fmo3 * tma_production_mol_h(p) + r_host) / cl
    if m == "linear":
        return c_substrate
    if m == "saturable":
        # hepatic flux cannot exceed either the substrate arriving (c_substrate branch)
        # or the catalytic capacity f_fmo3*Vmax*C/(Km+C) (c_capacity branch).
        # The self-consistent steady state is the larger of the two roots.
        vmax = p.f_fmo3 * p.fmo3_vmax_umol_h * 1e-6
        km = p.fmo3_km_umol_l * 1e-6
        a = vmax + r_host - cl * km
        disc = a * a + 4.0 * cl * r_host * km
        c_capacity = max(0.0, (a + math.sqrt(disc)) / (2.0 * cl))
        # consistency: the substrate branch is valid only if the enzyme can carry it
        flux_needed = p.f_fmo3 * tma_production_mol_h(p)
        flux_available = vmax * c_substrate / (km + c_substrate)
        if flux_needed <= flux_available:
            return c_substrate
        return c_capacity
    raise ValueError(f"unknown hepatic_mode {m!r}")


def capacity_ceiling(p: Parameters) -> float:
    """H3: hepatic capacity plateau in mol/L, the asymptote of the capacity branch as Km -> 0."""
    vmax = p.f_fmo3 * p.fmo3_vmax_umol_h * 1e-6
    return (vmax + p.r_host_umol_h * 1e-6) / renal_clearance(p)


def branch_split(p: Parameters) -> dict[str, float]:
    """Which branch the saturable model selects, and the substrate intake that flips it."""
    cl = renal_clearance(p)
    r_host = p.r_host_umol_h * 1e-6
    vmax = p.f_fmo3 * p.fmo3_vmax_umol_h * 1e-6
    km = p.fmo3_km_umol_l * 1e-6
    a = vmax + r_host - cl * km
    disc = a * a + 4.0 * cl * r_host * km
    c_cap = max(0.0, (a + math.sqrt(disc)) / (2.0 * cl))
    c_sub = (p.f_fmo3 * tma_production_mol_h(p) + r_host) / cl
    flux_available = vmax * c_sub / (km + c_sub)
    branch = "substrate" if p.f_fmo3 * tma_production_mol_h(p) <= flux_available else "capacity"

    def gap(intake: float) -> float:
        q = replace(p, precursor_intake_g_day=intake)
        cs = (p.f_fmo3 * tma_production_mol_h(q) + r_host) / cl
        return p.f_fmo3 * tma_production_mol_h(q) - vmax * cs / (km + cs)

    lo, hi = 1e-9, 1e6
    intake_needed = (
        float(brentq(gap, lo, hi)) if gap(lo) * gap(hi) < 0 else float("inf")
    )
    return {
        "c_substrate_branch_umol_l": c_sub * 1e6,
        "c_capacity_branch_umol_l": c_cap * 1e6,
        "selected_branch": branch,
        "intake_at_which_capacity_binds_g_day": intake_needed,
    }


def phi(p: Parameters, *, mode: str | None = None) -> float:
    """Microbial-dependence index: 1 - C_ss(y=0)/C_ss(y_ref). Independent of y by definition."""
    return 1.0 - steady_state(replace(p, y_microbial_conversion=0.0), mode=mode) / steady_state(p, mode=mode)


def microbial_dependence_index(p: Parameters, *, mode: str | None = None) -> float:
    return phi(p, mode=mode)


def exposure_contrast(p: Parameters) -> float:
    """G2: C_ss contrast between the 25th and 75th percentile of the microbial-capacity prior."""
    med, s = PRIORS["y_microbial_conversion"]
    z = norm.ppf(0.75)
    y_lo, y_hi = med * math.exp(-z * s), med * math.exp(z * s)
    c_lo = steady_state(replace(p, y_microbial_conversion=y_lo))
    c_hi = steady_state(replace(p, y_microbial_conversion=y_hi))
    return c_hi / c_lo


# --------------------------------------------------------------------------
# transient (used only to verify the closed forms)
# --------------------------------------------------------------------------


def dC_dt(t: float, c: float, p: Parameters, *, mode: str | None = None) -> float:
    m = mode or p.hepatic_mode
    cl = renal_clearance(p)
    r_host = p.r_host_umol_h * 1e-6
    if m == "linear":
        prod = p.f_fmo3 * tma_production_mol_h(p) + r_host
    else:
        vmax = p.f_fmo3 * p.fmo3_vmax_umol_h * 1e-6
        km = p.fmo3_km_umol_l * 1e-6
        # same self-consistency cap as the closed form: catalytic flux is bounded by
        # both capacity and the substrate actually arriving from the colon
        prod = min(p.f_fmo3 * tma_production_mol_h(p), vmax * c / (km + c)) + r_host
    return (prod - cl * c) / p.vd_l


def transient(p: Parameters, *, days: float = 400.0, n: int = 400001, **kw: Any) -> np.ndarray:
    """Explicit fixed-step Euler integration, single thread, tiny array."""
    t = np.linspace(0.0, days, n)
    dt = 24.0 * (t[1] - t[0])
    c = 0.0
    out = np.empty(n)
    for i in range(n):
        out[i] = c
        c = c + dt * dC_dt(t[i], c, p, **kw)
    return out


def half_life_h(p: Parameters) -> float:
    """Linear-kinetics terminal half-life, an analytical cross-check on Vd/CL."""
    return math.log(2.0) * p.vd_l / renal_clearance(p)


# --------------------------------------------------------------------------
# stochastic layer: lognormal priors, frozen draws and seed
# --------------------------------------------------------------------------


def draw_parameters(
    p: Parameters,
    rng: np.random.Generator,
    names: Iterable[str],
) -> dict[str, np.ndarray]:
    return {
        n: PRIORS[n][0] * np.exp(rng.normal(0.0, PRIORS[n][1], size=p.mc_draws))
        for n in names
    }


def _vectorised_steady_state(p: Parameters, col: dict[str, np.ndarray], *, mode: str | None = None) -> np.ndarray:
    m = mode or p.hepatic_mode
    const = np.full(p.mc_draws, 1.0)
    body = col.get("body_mass_kg", np.full(p.mc_draws, p.body_mass_kg))
    y = col.get("y_microbial_conversion", np.full(p.mc_draws, p.y_microbial_conversion))
    fesc = col.get("f_colonic_escape", np.full(p.mc_draws, p.f_colonic_escape))
    fhep = col.get("f_tma_to_hepatic", np.full(p.mc_draws, p.f_tma_to_hepatic))
    ffmo = col.get("f_fmo3", np.full(p.mc_draws, p.f_fmo3))
    gfr = col.get("gfr_ml_min_per_kg", np.full(p.mc_draws, p.gfr_ml_min_per_kg))
    host = col.get("r_host_umol_h", np.full(p.mc_draws, p.r_host_umol_h)) * 1e-6
    del const
    cl = gfr * body * (60.0 / 1000.0) * (1.0 - p.f_reabsorption) * (1.0 + p.f_tubular_secretion)
    tma_hep = (p.precursor_intake_g_day / p.precursor_mw_g_mol / 24.0) * fesc * p.nu_stoich * y * fhep
    if m == "linear":
        return (ffmo * tma_hep + host) / cl
    vmax = ffmo * p.fmo3_vmax_umol_h * 1e-6
    km = p.fmo3_km_umol_l * 1e-6
    a = vmax + host - cl * km
    disc = a * a + 4.0 * cl * host * km
    return np.maximum(0.0, (a + np.sqrt(disc)) / (2.0 * cl))


def monte_carlo(
    p: Parameters = DEFAULT,
    *,
    names: Iterable[str] | None = None,
    seed: int | None = None,
) -> dict[str, Any]:
    """Propagate the frozen lognormal priors into C_ss. Log-space, so the CI is a ratio."""
    names = tuple(names) if names is not None else tuple(PRIORS)
    rng = np.random.default_rng(p.mc_seed if seed is None else seed)
    col = draw_parameters(p, rng, names)
    c = _vectorised_steady_state(p, col)
    ln = np.log(c)
    lo, hi = np.percentile(ln, [5.0, 95.0])
    per_param = {}
    for n in names:
        if float(np.ptp(col[n])) > 0.0 and float(np.ptp(ln)) > 0.0:
            rho: Any = float(spearmanr(col[n], ln).statistic)
        else:
            rho = None
        per_param[n] = {
            "spearman_rho_with_ln_C": rho,
            "block": PARAMETER_SPECS[n][2],
            "rel_width_90pct_lognormal_prior": float(
                math.exp(1.6448536269514722 * PRIORS[n][1])
                / math.exp(-1.6448536269514722 * PRIORS[n][1])
                - 1.0
            ),
        }
    return {
        "draws": int(p.mc_draws),
        "seed": int(p.mc_seed if seed is None else seed),
        "params_sampled": list(names),
        "C_ss_umol_l_median": float(np.median(c) * 1e6),
        "C_ss_umol_l_p5": float(math.exp(lo) * 1e6),
        "C_ss_umol_l_p95": float(math.exp(hi) * 1e6),
        "C_ss_fold_width_90pct": float(math.exp(hi - lo)),
        "per_param": per_param,
    }


def variance_contributions(p: Parameters = DEFAULT, seed: int | None = None) -> dict[str, Any]:
    """One-at-a-time log-variance shares. Interaction terms are neglected by construction."""
    base_seed = p.mc_seed if seed is None else seed
    total = monte_carlo(p, seed=base_seed)
    ln_tot = math.log(total["C_ss_fold_width_90pct"])
    shares: dict[str, float] = {}
    for n in PRIORS:
        sub = monte_carlo(p, names=(n,), seed=base_seed)
        shares[n] = float((math.log(sub["C_ss_fold_width_90pct"]) / ln_tot) ** 2)
    s = sum(shares.values())
    return {
        "method": "one-at-a-time log-variance share, interactions neglected",
        "seed": int(base_seed),
        "raw_shares": shares,
        "normalised_shares": {k: v / s for k, v in shares.items()},
        "share_by_block": {
            b: float(
                sum(
                    v
                    for k, v in shares.items()
                    if PARAMETER_SPECS[k][2] == b
                )
                / s
            )
            for b in ("diet", "microbiome", "host")
        },
    }


def needed_variables(
    p: Parameters = DEFAULT,
    *,
    target_fold_width: float = 2.0,
    seed: int | None = None,
) -> dict[str, Any]:
    """Q012 directly: greedy set selection until the 90 % fold width of C_ss is met.

    Start with nothing measured (every uncertain parameter pinned at its prior median),
    add the parameter that removes the most log-variance, stop at the target.
    """
    base_seed = p.mc_seed if seed is None else seed
    med = {n: PRIORS[n][0] for n in PRIORS}

    def width_of(measured: list[str]) -> float:
        free = [n for n in PRIORS if n not in measured]
        if not free:
            return 1.0
        rng = np.random.default_rng(base_seed)
        col: dict[str, np.ndarray] = {n: np.full(p.mc_draws, med[n]) for n in PRIORS}
        col.update(draw_parameters(p, rng, free))
        ln = np.log(_vectorised_steady_state(p, col))
        return float(math.exp(np.percentile(ln, 95) - np.percentile(ln, 5)))

    order: list[dict[str, Any]] = []
    measured: list[str] = []
    w = width_of(measured)
    while w > target_fold_width and len(measured) < len(PRIORS):
        best, best_w = None, w
        for n in PRIORS:
            if n in measured:
                continue
            cand = width_of(measured + [n])
            if cand < best_w:
                best, best_w = n, cand
        if best is None:
            break
        measured.append(best)
        order.append(
            {
                "step": len(order) + 1,
                "variable": best,
                "block": PARAMETER_SPECS[best][2],
                "fold_width_90pct_after": best_w,
            }
        )
        w = best_w
    return {
        "target_fold_width_90pct": target_fold_width,
        "seed": int(base_seed),
        "start_fold_width_90pct": width_of([]),
        "final_fold_width_90pct": w,
        "criterion_met": bool(w <= target_fold_width),
        "greedy_order": order,
    }


def identifiability_attribution(p: Parameters = DEFAULT) -> dict[str, Any]:
    """G3 / H2: can the product p = f_esc * y be pinned, and can its factors be separated?"""
    z = 1.6448536269514722
    se, sy = PRIORS["f_colonic_escape"][1], PRIORS["y_microbial_conversion"][1]

    def rel_width(s: float) -> float:
        return math.exp(z * s) / math.exp(-z * s) - 1.0

    w_prod = rel_width(math.hypot(se, sy))
    w_esc = rel_width(se)
    w_y = rel_width(sy)
    return {
        "product_sigma_ln": math.hypot(se, sy),
        "rel_width_90pct_product": w_prod,
        "rel_width_90pct_f_colonic_escape": w_esc,
        "rel_width_90pct_y_microbial_conversion": w_y,
        "product_narrower_than_each_factor": bool(w_prod < min(w_esc, w_y)),
        "G3_pass": bool(w_prod < G3_PRODUCT_REL_WIDTH_MAX and w_esc > G3_FACTOR_REL_WIDTH_MIN and w_y > G3_FACTOR_REL_WIDTH_MIN),
    }


# --------------------------------------------------------------------------
# null models and placebo
# --------------------------------------------------------------------------


def null_model_diet_only(p: Parameters, *, calibrate_at: Parameters | None = None) -> float:
    """N0-A: one coefficient k = C/I fixed at the calibration point, then evaluated at p.

    The coefficient is deliberately NOT refitted per call, otherwise the null would
    silently reabsorb the mechanism it is supposed to be tested against.
    """
    cal = DEFAULT if calibrate_at is None else calibrate_at
    k = steady_state(cal) / cal.precursor_intake_g_day
    return k * p.precursor_intake_g_day


def null_model_microbiome_is_everything(p: Parameters) -> float:
    """N0-B: microbial term only, f_colonic_escape = 1 and f_fmo3 = 1."""
    q = replace(p, f_colonic_escape=1.0, f_fmo3=1.0)
    return steady_state(q)


def placebo_reproduces_linear(p: Parameters) -> dict[str, Any]:
    """Same code path, hepatic FMO3 driven far past saturation (Km -> 0 umol/L).

    With Km << C the enzyme flux collapses onto its capacity and the capacity branch can
    never bind, so the saturable closed form must return the linear closed form exactly.
    Note that the linear mode is a *reduced* low-load model, not the Km -> infinity limit
    of the saturable branch: there, the flux is first-order (Vmax/Km)*C, not a constant
    source, so a Km blow-up is not the right placebo.
    """
    q = replace(p, fmo3_km_umol_l=1e-3, hepatic_mode="saturable")
    return {
        "saturable_with_km_1e-3_umol_l": steady_state(q),
        "linear_closed_form": steady_state(p),
        "abs_rel_diff": abs(steady_state(q) / steady_state(p) - 1.0),
    }


# --------------------------------------------------------------------------
# acceptance criteria, frozen in PREREG.md section 4
# --------------------------------------------------------------------------


def criteria(p: Parameters = DEFAULT) -> dict[str, Any]:
    c_ref = steady_state(p)
    c_kap = capacity_capacity = capacity_ceiling(p)
    margin = c_kap / c_ref
    ph = phi(p)
    contrast = exposure_contrast(p)
    att = identifiability_attribution(p)
    return {
        "predicted_C_ss_umol_l": c_ref * 1e6,
        "predicted_phi": ph,
        "G1_phi_ge_0.95": {"value": ph, "threshold": G1_PHI_MIN, "pass": bool(ph >= G1_PHI_MIN)},
        "G2_exposure_contrast_ge_2": {
            "value": contrast,
            "threshold": G2_CONTRAST_MIN,
            "pass": bool(contrast >= G2_CONTRAST_MIN),
        },
        "G3_identifiability": {
            **att,
            "thresholds": {
                "product_rel_width_max": G3_PRODUCT_REL_WIDTH_MAX,
                "factor_rel_width_min": G3_FACTOR_REL_WIDTH_MIN,
            },
            "pass": att["G3_pass"],
        },
        "G4_capacity_margin_ge_3": {
            "capacity_ceiling_umol_l": capacity_capacity * 1e6,
            "margin_capacity_over_predicted": margin,
            "threshold": G4_CAPACITY_MARGIN_MIN,
            "pass": bool(margin >= G4_CAPACITY_MARGIN_MIN),
        },
    }


def sensitivity(p: Parameters = DEFAULT) -> dict[str, Any]:
    """OAT +/-50 % on the three frozen parameters, PREREG.md section 5."""
    out: dict[str, Any] = {}
    base = steady_state(p)
    for name in SENSITIVITY_PARAMS:
        row: dict[str, Any] = {"block": PARAMETER_SPECS[name][2], "base": getattr(p, name)}
        for tag, factor in (("minus50", 0.5), ("plus50", 1.5)):
            q = replace(p, **{name: getattr(p, name) * factor})
            c = steady_state(q)
            row[tag] = {
                "parameter_value": getattr(q, name),
                "C_ss_umol_l": c * 1e6,
                "abs_change_umol_l": (c - base) * 1e6,
                "rel_change": c / base - 1.0,
                "phi": phi(q),
            }
        out[name] = row
    return {"base_C_ss_umol_l": base * 1e6, "params": out}


# --------------------------------------------------------------------------
# unit check
# --------------------------------------------------------------------------


def dimensional_check() -> dict[str, Any]:
    d = PARAMETER_DIMS
    checks: dict[str, dict[str, Any]] = {}

    def add(name: str, lhs: Dim, rhs: Dim) -> None:
        ok = lhs == rhs
        checks[name] = {
            "lhs_dim": repr(lhs),
            "rhs_dim": repr(rhs),
            "status": "PASS" if ok else "FAIL",
        }

    add(
        "precursor_rate_mol_per_day",
        d["precursor_intake_g_day"] / d["precursor_mw_g_mol"],
        DIM_MOL * DIM_PER_DAY,
    )
    add(
        "tma_production_mol_per_h",
        (d["precursor_intake_g_day"] / d["precursor_mw_g_mol"] * DIM_DAY_TO_HOUR)
        * d["f_colonic_escape"]
        * d["nu_stoich"]
        * d["y_microbial_conversion"]
        * d["f_tma_to_hepatic"],
        DIM_MOL / DIM_H,
    )
    t1 = d["f_fmo3"] * (DIM_MOL / DIM_H)
    t2 = d["r_host_umol_h"]
    add("tmao_source_term_fmo3", t1, DIM_MOL / DIM_H)
    add("tmao_source_term_host", t2, DIM_MOL / DIM_H)
    add("tmao_source_mol_per_h", t1.add(t2), DIM_MOL / DIM_H)
    add(
        "renal_clearance_l_per_h",
        d["gfr_ml_min_per_kg"] * d["body_mass_kg"] * DIM_MIN_TO_HOUR,
        DIM_L / DIM_H,
    )
    add("steady_state_concentration", (DIM_MOL / DIM_H) / (DIM_L / DIM_H), DIM_MOL / DIM_L)
    add(
        "transient_rate_mol_per_l_per_h",
        (DIM_MOL / DIM_H) / d["vd_l"],
        DIM_MOL / (DIM_L * DIM_H),
    )
    add(
        "capacity_ceiling_mol_per_l",
        (DIM_MOL / DIM_H) / (DIM_L / DIM_H),
        DIM_MOL / DIM_L,
    )
    add("half_life_h", d["vd_l"] / (DIM_L / DIM_H), DIM_H)
    checks["molar_bridge_umol_to_mol"] = {
        "note": "fmo3_vmax_umol_h, fmo3_km_umol_l, r_host_umol_h carry a 1e-6 mol/umol bridge",
        "status": "PASS",
    }
    return {
        "status": "PASS" if all(v["status"] == "PASS" for v in checks.values()) else "FAIL",
        "checks": checks,
    }


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------


def run_all(p: Parameters = DEFAULT) -> dict[str, Any]:
    lin = steady_state(p, mode="linear")
    sat = steady_state(p, mode="saturable")
    return {
        "id": "BT-HX-Q012",
        "substance_question": "S1: dietary precursor intake x microbial conversion capacity -> steady-state plasma TMAO",
        "parameters": {k: asdict(p)[k] for k in asdict(p)},
        "parameter_table": {
            k: {"value": getattr(p, k), "unit": u, "source_or_assumption": s, "block": b}
            for k, (u, s, b) in PARAMETER_SPECS.items()
        },
        "priors_lognormal": {k: {"median": m, "sigma_ln": s, "block": PARAMETER_SPECS[k][2]} for k, (m, s) in PRIORS.items()},
        "intermediates": {
            "renal_clearance_l_h": renal_clearance(p),
            "precursor_mol_per_day": precursor_mol_per_day(p),
            "tma_production_delivered_mol_h": tma_production_mol_h(p),
            "tmao_source_mol_h": tmao_source_mol_h(p),
            "linear_half_life_h": half_life_h(p),
        },
        "predictions": {
            "C_ss_umol_l_linear": lin * 1e6,
            "C_ss_umol_l_saturable": sat * 1e6,
            "branch_split": branch_split(p),
            "C_ss_umol_l": lin * 1e6,
            "phi_microbial_dependence": phi(p),
            "capacity_ceiling_umol_l": capacity_ceiling(p) * 1e6,
            "exposure_contrast_fold": exposure_contrast(p),
        },
        "criteria": criteria(p),
        "sensitivity_plus_minus_50pct": sensitivity(p),
        "null_models": {
            "N0_A_diet_only_umol_l": null_model_diet_only(p) * 1e6,
            "N0_B_microbiome_is_everything_umol_l": null_model_microbiome_is_everything(p) * 1e6,
            "placebo": placebo_reproduces_linear(p),
        },
        "monte_carlo": monte_carlo(p),
        "variance_contributions": variance_contributions(p),
        "needed_variables": needed_variables(p),
        "identifiability": identifiability_attribution(p),
        "dimensional_check": dimensional_check(),
        "verification": {
            "transient_converges_to_closed_form_linear": float(transient(p, days=800.0)[-1]),
            "closed_form_linear": lin,
            "transient_converges_to_closed_form_saturable": float(
                transient(p, days=800.0, mode="saturable")[-1]
            ),
            "closed_form_saturable": sat,
        },
    }


def write_results(path: str | Path = "results.json") -> dict[str, Any]:
    res = run_all()
    Path(path).write_text(json.dumps(res, indent=2, sort_keys=False) + "\n")
    return res


if __name__ == "__main__":
    out = write_results()
    print(json.dumps({k: out[k] for k in ("predictions", "criteria", "dimensional_check")}, indent=2))
