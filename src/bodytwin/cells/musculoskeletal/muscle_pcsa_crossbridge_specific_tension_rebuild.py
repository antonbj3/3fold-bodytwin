#!/usr/bin/env python3
"""muscle_pcsa_crossbridge_specific_tension_rebuild.py., resolves
MODEL-MUSCLE-ARCHITECTURE-PCSA (cluster: model-mechanistic / MUSCULOSKELETAL-DYNAMICS -- converts a
molecular cross-bridge force budget into a whole-muscle specific-tension prediction).

STATUS BEFORE THIS SCRIPT: node's cert_design.verify said "SEED-DESIGN ... cited literature, not
executed". The node's text already names a fully-specifying parameter table (single myosin-head
force 3-4pN + 11nm step, Finer et al. 1994 Nature 368:113, PMID8139653; 294 myosin heads per thick
filament, filament length 1.6um, PMC5770512; d10 X-ray lattice spacing 35-40nm; isometric duty ratio
~0.20) and asserts a resulting sigma=22.3-29.8 N/cm^2 range (central 26.0) against two independent
literature anchors -- Persad et al. 2024 J Appl Physiol 137:945 (PMID39169839, PMC11486478, weighted
median 26.8 N/cm^2, IQR 20-43, n=96/30 studies) and Powell et al. 1984 J Appl Physiol (PMID6511546,
22.5 N/cm^2, r=0.99, n=26 guinea-pig muscles) -- but the ARITHMETIC was never executed as code (grep-
checked: no scripts/*/*.py computes a hexagonal-lattice-area-per-filament -> heads-per-CSA -> sigma
chain before this file; specific_tension_validation.py answers a DIFFERENT question, this
repo's gait-twin Fmax/PCSA decomposition, not this node's from-scratch molecular derivation).

GEOMETRY (derived here, not assumed): thick filaments pack on a 2D hexagonal lattice in cross-
section. X-ray diffraction measures d10, the spacing of the (1,0) Bragg planes of that hexagonal
lattice, which relates to the filament-to-filament nearest-neighbor spacing "a" by the standard
hexagonal-lattice identity d10 = a*sqrt(3)/2 (equivalently a = 2*d10/sqrt(3)). Each thick filament
occupies one hexagonal-lattice primitive cell, of area (sqrt(3)/2)*a^2 -- so:
    area_per_filament = (sqrt(3)/2) * (2*d10/sqrt(3))^2 = (2/sqrt(3)) * d10^2  ~= 1.1547 * d10^2
This is machine-verified below by an independent direct hexagonal-sum construction (place filaments
on an explicit hexagonal grid, Voronoi-cell area per filament by direct counting), not just trusted
algebra.

CHAIN: sigma [N/cm^2] = (heads_per_filament * duty_ratio * force_per_head) / area_per_filament
  heads_per_filament = 294 (PMC5770512), duty_ratio in [0.15,0.25] (isometric, node's band, from
  Finer1994/Molloy1995 duty-ratio literature reused per this node's cert_design), force_per_head
  in [3,4] pN (Finer1994), d10 in [35,40] nm (node's vertebrate bracket).

QUESTION (pre-registered before running): does this from-scratch geometric+force-budget derivation,
swept over the node's cited parameter RANGES (not a single fitted point), produce a sigma band
that contains BOTH the Persad-2024 (26.8 N/cm^2) and Powell-1984 (22.5 N/cm^2) anchors -- decorrelated
from each other and from every input to this derivation (X-ray/EM/optical-trap vs whole-muscle-force
measurement literature)?

GATE (pre-registered):
  G1: central estimate (d10=37.5nm mid, duty=0.20 mid, force=3.5pN mid) lands within [15,45] N/cm^2
      (this node's pre-registered falsifier band; a systematic-review median outside this would
      refute the derivation).
  G2: full parameter-sweep band [min,max] over the cited ranges CONTAINS both anchors (26.8, 22.5).
VOID FLOOR (pre-registered, node's disclosed adversary): using the SHORTENING duty ratio
  (0.01-0.05, from the SAME Finer1994/Molloy1995 source, but the wrong regime -- shortening, not
  isometric) must collapse sigma far below both anchors (node's claimed floor: 1.3-6.5 N/cm^2) --
  if it does NOT collapse, the isometric-vs-shortening duty-ratio distinction is not actually doing
  the work claimed and the whole derivation would be underdetermined by geometry alone.
"""
import numpy as np

PN_TO_N = 1e-12
NM_TO_CM = 1e-7  # 1 nm = 1e-7 cm
MYOSIN_MOLECULES_PER_FILAMENT = 294.0  # PMC5770512, thick filament ~1.6um length
HEADS_PER_FILAMENT = 2.0 * MYOSIN_MOLECULES_PER_FILAMENT  # 588: each myosin molecule is a DIMER
# with 2 heads (self-caught bug when this cell was written: 294 is molecules/dimers, not heads -- using 294
# directly as "heads" gave sigma_central=12.67 N/cm^2, ~1.9x below the falsifier band; the source
# material itself (search-corroborated, same PMC5770512 family) states "294 myosin dimers ...
# approximately 588 myosin heads total per thick filament" and separately "20-30% of the 294 myosin
# heads present in each HALF-thick-filament are attached" during isometric contraction -- confirming
# both the x2-heads correction AND this node's ~0.20 isometric duty-ratio estimate independently).

D10_RANGE_NM = (35.0, 40.0)          # vertebrate X-ray lattice spacing bracket (node's cite)
FORCE_PER_HEAD_RANGE_PN = (3.0, 4.0)  # Finer et al. 1994 Nature 368:113 PMID8139653
DUTY_ISOMETRIC_RANGE = (0.15, 0.25)   # node's "~0.20 isometric" band
DUTY_SHORTENING_RANGE = (0.01, 0.05)  # node's "shortening" void-floor band, same source family

ANCHOR_PERSAD_2024 = 26.8  # N/cm^2, PMC11486478, weighted median, n=96/30 studies
ANCHOR_POWELL_1984 = 22.5  # N/cm^2, PMID6511546, r=0.99 vs measured Po, n=26 guinea-pig muscles
FALSIFIER_BAND = (15.0, 45.0)

def area_per_filament_cm2_algebraic(d10_nm):
    """area_per_filament = (2/sqrt(3)) * d10^2, derived from d10 = a*sqrt(3)/2 (hexagonal lattice
    (1,0)-plane spacing vs nearest-neighbor spacing a), area of hex primitive cell = (sqrt(3)/2)*a^2."""
    d10_cm = d10_nm * NM_TO_CM
    return (2.0 / np.sqrt(3.0)) * d10_cm ** 2

def area_per_filament_cm2_direct_hexsum(d10_nm, n_side=25):
    """INDEPENDENT machine cross-check of the algebraic formula: build an explicit hexagonal lattice
    of filaments with nearest-neighbor spacing a = 2*d10/sqrt(3), and measure area-per-filament as
    (total lattice area covered) / (number of filaments), using a large patch to kill edge effects."""
    a_nm = 2.0 * d10_nm / np.sqrt(3.0)
    pts = []
    for i in range(-n_side, n_side + 1):
        for j in range(-n_side, n_side + 1):
            x = a_nm * (i + 0.5 * j)
            y = a_nm * (np.sqrt(3) / 2.0) * j
            pts.append((x, y))
    pts = np.array(pts)
    # bounding parallelogram area of the lattice patch / number of interior points approximates
    # area-per-lattice-point for a uniform lattice; use the exact primitive-cell determinant instead
    # for a noise-free cross-check (this IS the definition of lattice unit-cell area):
    v1 = np.array([a_nm, 0.0])
    v2 = np.array([a_nm * 0.5, a_nm * np.sqrt(3) / 2.0])
    cell_area_nm2 = abs(v1[0] * v2[1] - v1[1] * v2[0])
    # cross-check: count points strictly inside a large disk and compare to disk_area/cell_area
    R = a_nm * (n_side - 2)
    n_inside = int(np.sum(np.hypot(pts[:, 0], pts[:, 1]) < R))
    disk_area = np.pi * R ** 2
    empirical_cell_area = disk_area / n_inside
    return cell_area_nm2 * (NM_TO_CM ** 2), empirical_cell_area * (NM_TO_CM ** 2)

def sigma_n_per_cm2(d10_nm, force_per_head_pN, duty_ratio, heads=HEADS_PER_FILAMENT):
    area_cm2 = area_per_filament_cm2_algebraic(d10_nm)
    force_per_filament_N = heads * duty_ratio * force_per_head_pN * PN_TO_N
    pressure_N_per_cm2 = force_per_filament_N / area_cm2
    return pressure_N_per_cm2

def main():
    print("=== muscle_pcsa_crossbridge_specific_tension_rebuild: MODEL-MUSCLE-ARCHITECTURE-PCSA ===")

    # --- geometry self-check: algebraic vs direct hexagonal-lattice construction ---
    d10_mid = 37.5
    area_alg = area_per_filament_cm2_algebraic(d10_mid)
    area_cell_exact, area_cell_empirical = area_per_filament_cm2_direct_hexsum(d10_mid)
    rel_err_exact = abs(area_alg - area_cell_exact) / area_alg
    rel_err_emp = abs(area_alg - area_cell_empirical) / area_alg
    print(f"GEOMETRY CHECK d10={d10_mid}nm: algebraic area/filament={area_alg:.6e} cm^2, "
          f"exact-hex-cell={area_cell_exact:.6e} cm^2 (rel err {rel_err_exact:.2e}), "
          f"empirical-disk-count={area_cell_empirical:.6e} cm^2 (rel err {rel_err_emp:.2e}) "
          f"-> {'PASS (formula verified by independent construction)' if rel_err_exact < 1e-9 and rel_err_emp < 0.02 else 'FAIL'}")

    # --- G1: central estimate ---
    d10_c = 0.5 * sum(D10_RANGE_NM)
    force_c = 0.5 * sum(FORCE_PER_HEAD_RANGE_PN)
    duty_c = 0.5 * sum(DUTY_ISOMETRIC_RANGE)
    sigma_central = sigma_n_per_cm2(d10_c, force_c, duty_c)
    g1_pass = FALSIFIER_BAND[0] <= sigma_central <= FALSIFIER_BAND[1]
    print(f"\nG1 central estimate: d10={d10_c}nm force={force_c}pN duty={duty_c} "
          f"-> sigma={sigma_central:.2f} N/cm^2 (falsifier band {FALSIFIER_BAND}) "
          f"-> {'PASS' if g1_pass else 'FAIL'}")

    # --- G2: full sweep over cited ranges, isometric duty ratio ---
    d10s = np.linspace(*D10_RANGE_NM, 6)
    forces = np.linspace(*FORCE_PER_HEAD_RANGE_PN, 6)
    duties = np.linspace(*DUTY_ISOMETRIC_RANGE, 6)
    sigmas = np.array([sigma_n_per_cm2(d, f, u) for d in d10s for f in forces for u in duties])
    sigma_min, sigma_max = sigmas.min(), sigmas.max()
    contains_persad = sigma_min <= ANCHOR_PERSAD_2024 <= sigma_max
    contains_powell = sigma_min <= ANCHOR_POWELL_1984 <= sigma_max
    print(f"\nG2 full sweep (isometric duty {DUTY_ISOMETRIC_RANGE}, d10 {D10_RANGE_NM}nm, "
          f"force {FORCE_PER_HEAD_RANGE_PN}pN): sigma in [{sigma_min:.2f}, {sigma_max:.2f}] N/cm^2")
    print(f"  contains Persad-2024 anchor (26.8)? {'YES' if contains_persad else 'NO'}")
    print(f"  contains Powell-1984 anchor (22.5)? {'YES' if contains_powell else 'NO'}")
    g2_pass = contains_persad and contains_powell

    # --- void floor: shortening duty ratio instead of isometric ---
    duties_short = np.linspace(*DUTY_SHORTENING_RANGE, 6)
    sigmas_short = np.array([sigma_n_per_cm2(d, f, u) for d in d10s for f in forces for u in duties_short])
    s_short_min, s_short_max = sigmas_short.min(), sigmas_short.max()
    void_below_both = s_short_max < min(ANCHOR_PERSAD_2024, ANCHOR_POWELL_1984)
    print(f"\nVOID FLOOR (shortening duty ratio {DUTY_SHORTENING_RANGE}, same geometry/force): "
          f"sigma in [{s_short_min:.2f}, {s_short_max:.2f}] N/cm^2 "
          f"-> {'PASS (correctly collapses below both anchors)' if void_below_both else 'FAIL (does not discriminate)'}")

    verdict = "CONFIRMED" if (g1_pass and g2_pass and void_below_both and rel_err_exact < 1e-9) else "PARTIAL"
    print(f"\nVERDICT: {verdict}")
    print(f"DATAPOINTS: sigma_central={sigma_central:.2f}, sigma_sweep_min={sigma_min:.2f}, "
          f"sigma_sweep_max={sigma_max:.2f}, anchor_persad2024={ANCHOR_PERSAD_2024}, "
          f"anchor_powell1984={ANCHOR_POWELL_1984}, void_floor_max={s_short_max:.2f}, "
          f"geometry_formula_rel_err={rel_err_exact:.2e}")

if __name__ == "__main__":
    main()
