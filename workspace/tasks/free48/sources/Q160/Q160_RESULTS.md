BT-HX-Q160

# RESULTS — Q160 "Biofilm EPS mechanics and transport"

PREREF: `PREREG.md`, sha256 `4d705e6256accf61cc89083f9e7c9e5fa13218d51ba8cd4c29048000cf8b4aad`
(frozen 2026-09-25 13:01 UTC, **before the first run**). All numbers: `results.json`.
Model: `model.py` (1-D elasto-poro-EPS biofilm, 17 equations, see docstring).
Tests: `test_model.py` — **16/16 PASS**, including four analytical limiting cases.
Run: `run_experiment.py` → `results.json`, log `run.log`.

---

## 1. Prediction against reference

The model's mature biofilm (assumed composition φ_cell = 0,25, φ_EPS = 0,10):
**ε = 0,650**, **r_m = 27,1 nm**, **D_eff(fluorescein) = 2,01·10⁻¹⁰ m²/s**,
**κ = 3,67·10⁻¹⁷ m²**.

| ID | Quantity | Prediction | Reference (source) | Criterion | Outcome |
|----|---------|--------------|------------------|-----------|--------|
| **C1** | `D_eff/D_0`, fluorescein | **0,373** (−6,7 %) | **0,40** — Takenaka S, Pitts B, Trivedi HM, Stewart PS (2009), *Appl Environ Microbiol* **75**(6):1750–1753, DOI [10.1128/AEM.02279-08](https://doi.org/10.1128/AEM.02279-08), PMCID PMC2655469, Table 2. **Verified** (full text acquired 2026-09-25). | 0,20 ≤ x ≤ 0,80 | **MET** |
| **C2** | Size dependence, 10 solutes | 2/10 within factor 2; Spearman ρ = **0,976** | same DOI, Table 2 | ≥7/10 **and** ρ ≥ 0,7 | **NOT MET** (quantitatively) — the ρ part is met |
| **C3** | The trap: coupled mechanics vs null model N1 | RMS(log₁₀ t90): coupled 1,165 vs null 0,949 → **−22,7 %** | same DOI, Table 1 (119 pairs, 13 clusters) | ≥20 % better | **NOT MET** |

Robustness for C1: over the entire assumed composition interval
(φ_cell 0,15–0,30, φ_EPS 0,02–0,30), the model gives `D_eff/D_0` = **0,113–0,634**.
The band's lower edge is thus not reached over the entire interval ⇒ C1 holds **only for the
assumed composition**, not robustly. This is a dependence on an assumption, not a measured value.

## 2. What C2 and C3 actually say (mechanism, not score)

**C2 rejects a specific hypothesis, not the model.** Renkin hindrance is weak at
r_m = 27 nm (H = 0,94 for fluorescein, 0,85 for IgG — see `results.json:C2_diagnosis`).
Size dependence is thus carried almost **entirely by the size-independent tortuosity factor
ε_a²** (Cohan). A size-independent number cannot simultaneously give fluorescein 0,92 and
the dextrans 0,32–0,40 of the measured value. **Conclusion: wrong physics is localised — size-dependent
sorption/exclusion in the EPS network is missing, not tortuosity.** The fact that the model
underestimates the medium-sized substances by a factor 2,5–3,1 while matching fluorescein
(0,92) and IgG (0,75) is a diagnosis, not a curve fit.

**C3 is the question's own trap and it triggered.** The coupled mechanics makes
the `t90` curves **worse** (−22,7 %) than the null model (static matrix, constant D_eff).
The model thus does **not** explain the held-out penetration and detachment curves.
It is a negative result and is reported as one.

**Controls N2/N3:** N2 (yield disabled) gives identical `D_eff/D_0` = 0,3731 ⇒
the yield logic does **not** affect D_eff (it affects only detachment, see §3).
N3 (EPS removed, ε recalculated to 0,75) gives 0,5136 ⇒ **EPS lowers D_eff by 27 %**
and lowers κ simultaneously. It is the only quantified EPS effect that holds.

**H1's "interior optimum" is falsified in the D_eff part.** `eps_tradeoff_scan` shows that
`D_eff/D_0` is **monotonically decreasing** in φ_EPS: 0,514 (φ_EPS = 0) → 0,347 (0,10) →
0,042 (0,42). No optimum. The opposing direction exists only between D_eff/κ (minus) and
cohesion τ_y (plus), and since τ_y = 5 + 200·φ_EPS = 25 Pa ≫ τ_b = 0,5 Pa, it gives
**no observable detachment at all**.

## 3. Outputs according to the question

| Output | Value | Status |
|---|---|---|
| Biofilm thickness `L` | 5,86·10⁻³ m after 12 h | **FALSIFIED (upper bound).** Front speed becomes 0,136 m/s — ~5 orders of magnitude above observed biofilm growth (µm/h). The growth closure is not substrate-limited. |
| Effective diffusion `D_eff` | 2,01·10⁻¹⁰ m²/s (fluorescein) | Verified against literature (C1). |
| Detachment flux `J_det` | 0 cells/(m²·s) | **NOT MEASURABLE with frozen τ_y.** τ_y = 25 Pa exceeds τ_b = 0,5 Pa; the model gives 0. No published detachment reference was found within the time budget ⇒ **UNKNOWN**. |

Time scales (calculated, not assumed): `c_v` = κ·K_o/(µ(1−ε)) = **1,05·10⁻⁹ m²/s**,
drainage time L²/c_v = **0,024 s**, doubling time = 3,5 s ⇒ ratio **6,9·10⁻³**.
The mechanics is ~150× faster than the biology, which justifies using only
the quasi-static asymptote in the driver run; the transient solver exists in
`consolidate` and is verified against an exact Terzaghi series (`test_model.py` T4, error 0,001 %).

## 4. Sensitivity ±50 % (which 2–3 govern)

| Parameter | Δ`D_eff/D_0` | Δκ | Roll |
|---|---|---|---|
| **φ_cell** (0,25) | **0,83** | **6,80** | Largest. Porosity governs both, through ε² and ε³ respectively. |
| **φ_EPS** (0,10) | **0,35** | **2,52** | Second largest; the EPS recipe. |
| **r_0** (50 nm) | 0,08 | **2,00** | Strongly governs κ, weakly D_eff (γ is small). |

`Y_EPS` and `τ_yE` govern only `L` (0,60 and 0,50 respectively) — they are growth-rate and
cohesion parameters, not transport parameters. `J_det` is undefined (zero) at
baseline, so its sensitivity is **UNKNOWN** until τ_b is scaled so that detachment starts.

## 5. Next resolution step

1. **Measure ε and r_m, don't assume them.** PFG-NMR or FRAP of a reference dextran
   with known size gives r_m directly from Renkin inversion; CLSM with lectin/SYTO staining
   gives ε. Staudt et al. (2004, DOI 10.1002/bit.20241) has the method but **the numerical value could
   not be acquired** from the abstract ⇒ the porosity interval here is an ASSUMPTION, marked
   UNKNOWN in PREREG §5/R3.
2. **Size-dependent sorption.** Partition coefficient K_d(c_s) for fluorescein vs
   dextran in EPS is exactly the quantity C2 identifies. Without it, no tortuosity model
   is right for both ends.
3. **Scale detachment in.** τ_b must exceed τ_y; measure τ_y(φ_EPS) (AFM shear or
   macroscopic low-stress testing). Without this, `J_det` is a dead branch.
4. **Growth closure.** Replace `v_f = µ·δ/φ_c` with a substrate-limited field so that
   the front becomes O(µm/h); thickness in its current form is an upper bound, not a
   prediction.
5. **3D.** The velocity field near the substrate is 3D; the 1D slab is an average representation
   across EPS and cell channels.
6. **Cross-couplings** (the question): Q159 (separate biofilm from luminal bacteria —
   `Y_eps_surface` and `ρ_cell` need that boundary), Q161/Q149 (the time scales
   c_v-vs-growth and mechanics–transport).

## 6. What was NOT done

No calibration against internal data (no bodytwin repo exists in the sandbox).
No 3D-FEM. No fabricated measured data. No verdict words: C2 and C3 stand as **not
met** and every open question above is a node.
