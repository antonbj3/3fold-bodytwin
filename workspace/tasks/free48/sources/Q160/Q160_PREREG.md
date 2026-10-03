# PREREG — BT-HX-Q160 "Biofilm EPS mechanics and transport"

Frozen **before the first run**. Hash: see `PREREG.sha256`. All numbers below are inputs, not results.

## 0. Nod

`BT-HX-Q160` expands to: `BT-HX-Q160/P1` (EPS→geometry→D_eff), `BT-HX-Q160/P2`
(EPS→geometry→kappa), `BT-HX-Q160/P3` (geometry→yield stress→detachment),
`BT-HX-Q160/P4` (cross-coupling to Q149/Q159/Q161).

## 1. Builds on

- `inputs/QUESTION.md` (Q160), `inputs/NIGHT_PREAMBLE.md`.
- **No internal code anchor.** `~/projects/bodytwin` does not exist in this sandbox (`ls` → missing),
  `MECHANISM_ANCHOR_GRAPH.json` is therefore unreadable. The task itself says "Ingen direkt kodankre
  kartlagd i denna bounded review". No `import *` from others' code.
- Published literature, verified online this session (see §3).

## 2. Not redone

- No 3D-FEM. 1D slab (x through the thickness). Reason: ≤45 min, 1 thread.
- No stochastic/agent-based detachment model; detachment = continuous erosion rate.
- No data calibration against internal data. No fabricated measured data.
- No Explicit-the reference model comparison (the question is not musculoskeletal; BT-B24/N12 does not apply to Q160).

## 3. Hypothesis (mechanism, not correlation)

> **H1 (shared moving geometry).** EPS does not act through an arbitrary scalar "D-reducing factor".
> EPS and deformation change **one single geometric quantity — channel radius r_m [m]** — and from the same
> geometry *both* hydraulic permeability κ (∝ ε³r_m²) and effective diffusivity are derived:
> D_eff (through pore hindrance λ_s/r_m). The consequence is a **coupled, opposing** effect:
> more EPS ⇒ smaller r_m ⇒ **lower κ AND lower D_eff AND higher cohesion (lower yield stress)**
> ⇒ an interior optimum in the EPS fraction.
>
> **H2 (irreversible compaction).** Bulk shear stress τ_b exceeds the matrix yield stress
> τ_y(φ_EPS) ⇒ irreversible compaction (ε falls, cannot be restored by diffusion on
> the hour scale) ⇒ κ collapses and D_eff falls. It is the same yield stress that governs detachment, so
> **the same parameter appears in the transport and detachment parts** — this is the core of the coupling.
>
> **H3 (the trap).** If coupled mechanics does not give better penetration and detachment curves than
> a static matrix with constant D_eff, H1 is false. This is tested, not assumed.

## 4. Predicted quantity (predictand)

**Primary:** `D_eff / D_0` for **fluorescein** (MW 376 g/mol) in a mature biofilm.
Dimensionless ratio. Calculated by the model from porosity, sorption, channel radius and Renkin hindrance.

**Secondary (model outputs, not primary criterion):** biofilm thickness `L [m]`,
hydraulic permeability `κ [m²]`, detachment flux `J_det [celler/(m²·s)]`.

## 5. Reference value (looked up, cited)

**R1 — primary.**
> Takenaka S, Pitts B, Trivedi HM, Stewart PS (2009). *Diffusion of macromolecules in model
> oral biofilms.* Applied and Environmental Microbiology **75**(6):1750–1753.
> DOI: **10.1128/AEM.02279-08** — PMID 19168660, PMCID PMC2655469.
> **Fluorescein (376 g/mol): D_e = 219 ± 103 µm²/s, D_aq = 540 µm²/s, D_e/D_aq = 0,40** (n = 13 clusters).
> Acquired: `https://pmc.ncbi.nlm.nih.gov/articles/PMC2655469/` (full text, Table 2), 2026-09-25.
> **Verified (not OVERIFIERAD).**

**R1b — same source, size dependence (secondary predictand).**
D_e/D_aq: dextran 3 kDa 0,90 | 10 kDa 0,72 | 40 kDa 0,62 | 70 kDa 0,56 | papain 0,73 |
ficin 0,68 | GFP 0,76 | ConA 0,57 | IgG 0,22. D_aq [µm²/s]: 540, 186, 108, 55, 47, 96, 94,
87, 58, 44,5. Table 1 in the same paper gives `t90 [s]` per cluster and cluster radius R [µm]
(13 clusters, R = 44…147 µm). Source: same DOI.

**R2 — context, NOT used as a criterion value.**
Stewart PS (2003). *Diffusion in biofilms.* J Bacteriol **185**(5):1485–1491.
DOI 10.1128/JB.185.5.1485-1491.2003. Used only to support the tortuosity form.

**R3 — composition method, NOT numerical value.**
Staudt C, Horn H, Hempel DC, Neu TR (2004). *Volumetric measurements of bacterial cells and
extracellular polymeric substance glycoconjugates in biofilms.* Biotechnol Bioeng
**88**(5):585–592. DOI 10.1002/bit.20241. **The numerical value could not be acquired from the abstract →
the porosity interval below is AN ASSUMPTION, not a cited measurement. Marked UNKNOWN.**

## 6. Frozen acceptance criteria (binding, frozen before the run)

| ID | Storhet | Kriterium | Granskas mot |
|----|---------|-----------|--------------|
| **C1** | `D_eff/D_0`, fluorescein, mature biofilm | **0,20 ≤ pred ≤ 0,80** (factor 2 around R1 = 0,40) | R1 |
| **C2** | Size dependence: `D_eff` for all 10 solutes | **≥ 7 of 10** within factor 2 of measured `D_e`; **and** Spearman correlation(pred `D_e`, measured `D_e`) ≥ 0,7 | R1b |
| **C3** | The trap (the question's own): coupled mechanics vs null model | RMS error in predicted `t90` must **decrease ≥ 20 %** compared with the null model (static matrix, constant `D_eff`, no compaction) | R1b Table 1 |

**The C1 band is freely chosen (factor 2), NOT fitted to model outcomes.** If the model falls outside
the band, it is reported as **NOT MET** and the criterion is not changed afterwards.

## 7. Null models / placebo

- **N1 (the question's "Comparison"):** fixed biofilm volume, constant diffusivity. Porosity and r_m
  are frozen at t=0; no compaction; no EPS loss. Compared against C3.
- **N2 (placebo):** same model but `τ_y → 10¹² Pa` (compaction disabled). The N1/N2 difference
  isolates the yield logic; if they are equal, the yield logic is dead code.
- **N3 (pure EPS-off):** `φ_EPS = 0` ⇒ r_m = r_0, cohesion = the cells'. Isolates the EPS effect.

## 8. What counts as error

1. C1 outside [0,20; 0,80] ⇒ C1 NOT MET (reported as an error, not adjusted).
2. C2 < 7/10 ⇒ size dependence is not explained ⇒ H1 partly false.
3. C3 not met ⇒ **H1/H2 false for detachment and penetration curves** (the question's trap).
4. N1 ≈ N2 ⇒ the yield parameter has no effect ⇒ reported as dead code.
5. Any dimensional analysis that does not match ⇒ stop, no run.
6. Analytical limiting case not reproduced ⇒ stop, no result reporting.

## 9. Preselected parameters with allowed uncertainty

Entry from literature or assumption; **no values adjusted after the run**.

| Symbol | Value | Unit | Source / assumption |
|---|---|---|---|
| `D_0` per substance | 44,5…540 | µm²/s | **R1b Table 2** (the source's own values) |
| `λ_s` | from Stokes–Einstein from `D_0` | m | Derivation from R1b, not assumption |
| `τ_b` (bulk shear stress) | 0,5 | Pa | ASSUMPTION (channel flow typically 10⁻³ m/s, µ = 10⁻³ Pa·s) |
| `φ_EPS` mature biofilm | 0,10 | – | ASSUMPTION (interval 0,02–0,30, see R3 → UNKNOWN) |
| `φ_cell` mature biofilm | 0,25 | – | ASSUMPTION, same interval logic |
| `ε` | 0,65 (=1−0,25−0,10) | – | derived from the above |
| `r_0` (channel radius without EPS) | 50 | nm | ASSUMPTION |
| `β_c` (EPS blocking) | 3,0 | – | ASSUMPTION |
| `φ_bound` (sorption, volume fr.) | 0,02 | – | ASSUMPTION (Takenaka 2009 suggests sorption as an explanation for the fluorescein deviation) |
| `μ` water | 1,0·10⁻³ | Pa·s | ASSUMPTION (20 °C) |
| `K_o` (oedometer modulus) | 1,0·10⁴ | Pa | ASSUMPTION (biofilm E, order of magnitude) |
| `τ_y0` | 5,0 | Pa | ASSUMPTION |
| `τ_yE` | 200 | Pa | ASSUMPTION (the EPS network's contribution) |
| `k_er` | 1,0·10⁻⁶ | m/(Pa·s) | ASSUMPTION |
| `m` (shear exponent) | 1,0 | – | ASSUMPTION |
| `Y_EPS` | 0,35 | g/g cells | ASSUMPTION (GPS-like) |
| `ρ_cell` | 10¹² | cells/m³ | ASSUMPTION |
| `q_max` | 0,2 | 1/s | ASSUMPTION |
| `K_O2` (half saturation) | 5 | µmol/L | ASSUMPTION |
| `Y_O2` | 8,6·10⁻⁴ | mol O₂/mol cell | ASSUMPTION (cells ≈ 3·10¹¹/mL ≈ 3·10¹⁷/m³) |

## 10. What the next step requires

See `RESULTS.md` §"Next resolution step". Briefly: (a) measured ε and r_m directly (CLSM + PFG-NMR),
(b) measured yield stress τ_y(φ_EPS), (c) scale points for detachment flux, (d) 3D geometry.
