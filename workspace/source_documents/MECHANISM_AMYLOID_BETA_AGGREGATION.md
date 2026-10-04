# MECHANISM AMYLOID-BETA AGGREGATION — Knowles/Cohen master-equation kinetics, Aβ42-vs-Aβ40, CSF/PET decorrelated anchor (2026-07-22)

New hard-disease domain: a **quantitative, geometry-derived model of Aβ42 aggregation kinetics**
(primary nucleation / elongation / fibril-surface-catalyzed secondary nucleation — the
Knowles/Cohen master-equation framework), attacked as an overshoot-target adversary. Script:
`scripts/msk/amyloid_beta_aggregation.py`. Raw machine-checked results:
`reports/probes/amyloid_beta_aggregation_results.json`. Consolidated citation ledger:
`docs/MECHANISM_AMYLOID_BETA_AGGREGATION_evidence.json`.

**Scope note before anything else** (confirmed by direct read of `data/MECHANISM_ANCHOR_GRAPH.json`,
not assumed): this repo already carries `NEU-ALZHEIMERS-MULTIMODAL` (OPEN, cluster
`cancer-immune-wargame`), which triangulates the amyloid-hypothesis-vs-tau-vs-neuroinflammation
**causality** debate (anti-amyloid RCT effect sizes, Braak/plaque-vs-cognition, APOE4/TREM2
genetics) — that node explicitly holds the Lesné-fraud and RCT-effect-size adjudication and is
**not re-litigated here**. It does *not* contain a quantitative aggregation-kinetics forward
model. This doc is a **companion**: it builds and falsifies the mechanistic master equation
itself — the piece that node's own `couples_to`/`decorrelated_anchor` prose gestures at but never
constructs — plus a fresh CSF/amyloid-PET decorrelated leg. `NEU-TAU-PROPAGATION` (OPEN) and
`NEU-GLIA-NEUROINFLAMMATION` (OPEN, glymphatic-clearance sub-cell already gates
"amyloid/metabolite clearance rate by sleep-stage architecture") are the existing downstream/
clearance couplings this doc plugs into (§ couples_to), not duplicates.

## Part 1 — The master equation, derived from the geometry (not curve-fit)

### Method, in one paragraph

Unseeded amyloid self-assembly has two structurally different feedback topologies. **Primary
nucleation alone** (monomer → new fibril, rate `kn·mⁿᶜ`) has no dependence on existing fibril
mass `M` — the number of growth-competent ends `P` grows independently of `M`, so mass grows only
**polynomially** in time; there is no autocatalytic loop. **Secondary nucleation** (existing
fibril *surface* catalyzes new-aggregate formation from monomer, rate `k2·mⁿ²·M`) closes a
feedback loop: more `M` → more new `P` → more `M`. Linearizing the `(P,M)` system at early time
gives a 2×2 Jacobian whose eigenvalue is

```
κ = sqrt( 2 · (k+·k2) · m0^(n2+1) )         =>  M(t) ~ exp(κt)  at small t
```

so the aggregation half-time inherits κ's own concentration dependence: **t½ ~ κ⁻¹ ~
m0^-(n2+1)/2**. This exponent — not a fitted curve shape — is the geometric signature this script
tests. A rescaling `P' ≡ k+·P` turns the ODE into one that depends **only on the combined rate
constants** `(k+·kn)`, `(k+·k2)` — exactly the quantities the source papers report, and exactly
why unseeded kinetics alone cannot resolve `k+` individually (verified live in Meisl 2014 PNAS'
own stated methodological limitation, see evidence JSON) — so no elongation-rate assumption is
needed anywhere in this model.

**Geometric self-consistency check** (two independent computational paths, machine-compared, not
eyeballed): the numerically-fit early-time exponential growth rate of `P'(t)` (pure-secondary
limit) vs the closed-form eigenvalue above, at `m0=4µM`: fit = `0.00226274 s⁻¹`, theory =
`0.00226274 s⁻¹`, **relative error 1.2×10⁻⁷ — PASS** (threshold: <0.02). This formula is then
checked against the source papers' **own** independently-published formula (not this script's
invention): Cohen 2013 PNAS states `γ = -(n2+1)/2` for the pure-secondary limit and Meisl 2014
PNAS derives, from its own Eq. 11 (`γ = -½·(n2/(1+[m]0ⁿ²/KM) + 1)`), the identical two limits
(`-(n2+1)/2` at low `[m]0`, `-1/2` once secondary nucleation saturates) — **two independently
authored papers, six years apart, using the same functional form this script derived from
scratch.**

## Part 2 — Falsifier: does the model reproduce the MEASURED Aβ42 scaling exponent?

**Pre-registered claim (C):** a mechanistic ODE integrated with Cohen et al. 2013 PNAS' own
globally-fitted combined rate constants (their Fig. 1 panel "D", the full model:
`√(k+kn)=30 M⁻¹s⁻¹`, `√(k+k2)=2×10⁵ M⁻³ᐟ²s⁻¹`, `nc=n2=2` — **not tuned by this script**) reproduces
a log-log half-time-vs-concentration slope consistent with the **measured** exponent (Cohen 2013
PNAS Fig. 1A real ThT data: **γ = -1.33 ± 0.03**), and this is measurably distinguishable from a
mechanism lacking the fibril-surface feedback loop.

**Pre-registered threshold:** simulated γ ∈ [-1.7, -1.1] (bracketing both the measured value and
the pure-theory -1.5 with margin).

**The adversary, forced to its strongest fair form:** Cohen et al. 2013 PNAS itself fitted **two
no-secondary-nucleation alternative models** to the same data — reported in their own Fig. 1
legend, not invented here — (B) primary-nucleation-only with an *elevated* reaction order `nc=3`
(`√(k+kn)=8×10³ M⁻³ᐟ²s⁻¹`, `k2=0`) — the paper's own best attempt to let primary nucleation alone
explain a steep slope — and (C) primary nucleation (`nc=2`) + fragmentation (no fibril-surface
catalysis). This script runs adversary (B) — the strongest no-secondary-nucleation fit the source
paper produced — through the **identical ODE code**, no separate adversary-specific machinery.

**Results (machine-computed, `reports/probes/amyloid_beta_aggregation_results.json`):**

| quantity | value | anchor / threshold | verdict |
|---|---|---|---|
| True-model (secondary-nucleation) simulated γ, `m0`=1-8µM | **-1.320** (R²=0.99994, 8 pts) | measured -1.33±0.03; theory -1.5; band [-1.7,-1.1] | **PASS** |
| Adversary-B simulated γ (primary-only, nc=3, k2=0) | **-1.500** (R²=1.00000 exactly*) | — | see below |
| Slope-alone gap \|γ_true − γ_adversary\| | 0.180 | — | slope alone is NOT tied here (0.18 ≫ measured σ=0.03) |
| Shape metric: normalized 10→90% rise / t½ @ 4µM, true model | **0.890** | — | — |
| Shape metric: same, adversary B | **2.100** | — | — |
| Shape-metric relative gap | **57.6%** | pre-registered ≥20% | **PASS — adversary's transition is 2.4x more gradual** |

`*` R²=1.00000 for adversary B is not a fitting artifact: a single-power-law ODE with no second
term is **exactly scale-invariant** (the system `dP'/dt=kn_comb·mⁿᶜ, dM/dt=2mP'` has no timescale
other than the one set by `m0` itself, so *any* conversion-fraction threshold, not just the
half-time, obeys a perfect power law) — confirmed algebraically, not just numerically, as part of
this session's own derivation (§ evidence JSON).

**Honest reframing of the adversary (symmetric QC on my own pre-registration):** I pre-registered
an expectation that the slope alone might be *degenerate* (an adversary could always pick some
`nc` to mimic any target slope) — that did **not** concretely materialize in this specific,
non-reverse-fitted comparison (the two models' own independently-reported orders give slopes 0.18
apart, resolvably different at the measured ±0.03 precision). The **conceptual** vulnerability
remains real, though (nothing stops an adversary from choosing `nc` to exactly match -1.33), which
is exactly why Cohen et al. 2013 PNAS themselves state a no-secondary-pathway fit "is not able to
describe even qualitatively **the full time courses**" — not just the slope. The shape metric
above is the machine-computed form of that same claim, and it is the more decisive, non-gameable
discriminator: **PASS**, forced-adversary falls.

### Honest gaps (Part 2)
- The elongation rate constant `k+` is never individually assumed anywhere in this model (the
  rescaling removes the need) — but this also means this script cannot independently verify `k+`
  itself; only the *combined* constants and the reaction orders are tested.
- Cohen 2013 PNAS' exact ThT concentration grid for Aβ42 was not extracted this session (only
  "10 concentrations," "4µM representative," and "no saturation up to 6µM" were machine-verified
  text/numbers); this script's own 1-8µM grid is representative, not a pixel-exact replica of
  their Fig. 1A points.
- The absolute simulated half-times (0.3–4.6 hours across 8→1µM) are in the physically plausible
  range for a ThT Aβ42 assay — a sanity check, not an independently cited absolute-timescale
  anchor (only the *combined* rate constants, not a standalone timescale reference, were verified
  live).
- Numerical resolution: interpolated `t½` values carry a few-percent grid-discretization
  sensitivity (checked directly: re-running with a different `stop_frac` shifted one reported
  `t½` by ~3%); every pre-registered threshold band in this doc has margin ≥10x that noise floor,
  so no PASS/FAIL verdict here is sensitive to it — disclosed, not hidden.

## Part 3 — Aβ42-vs-Aβ40 propensity difference

**Claim (C):** Aβ42 aggregates substantially faster than Aβ40 at matched concentration, and this
difference is disproportionately carried by a **shift in mechanism** (relatively weaker primary
nucleation in Aβ40, not just uniformly slower rates) — Meisl et al. 2014 PNAS, verified live
(PMID 24938782): "the significant differences...originate...from a shift of **more than one order
of magnitude** in the relative importance of primary...versus...secondary nucleation," on top of
combined rate constants that are themselves "**both an order of magnitude smaller**" for Aβ40.
Aβ40's secondary nucleation is also **concentration-dependently saturating** (Michaelis-Menten-like,
`KM = 31 µM²`, half-saturation at 5-6µM — a real, live-verified number, not assumed), unlike
Aβ42 ("no saturation is observed at concentrations up to 6µM").

**Reconstruction (disclosed as order-of-magnitude, not digit-exact):** Cohen 2013's Fig.1-panel-D
combined constants, divided by 10 (`k+k2`) and by 100 (`k+kn` — the extra order of magnitude for
the "relative importance" shift), plus the verified `KM=31µM²`, `nc=n2=2` (Meisl 2014 explicitly
fixed these equal to Aβ42's values "to allow direct comparison"). Fig. 3B's exact fitted table is
a rendered figure/table graphic in PMC, not machine-extractable text this session — an honest,
disclosed gap; the reconstruction tests **direction and rough magnitude**, not exact digits.

**Pre-registered threshold:** simulated `t½(Aβ40)/t½(Aβ42)` at matched concentration (3.5µM, the
lowest concentration both peptides were tested near) ∈ [3, 300] (deliberately wide, given the gap
above).

| quantity | value | verdict |
|---|---|---|
| `t½` Aβ42 @ 3.5µM (simulated) | 3228 s (0.90 h) | — |
| `t½` Aβ40 @ 3.5µM (simulated) | 15824 s (4.40 h) | — |
| **Fold difference** | **4.90x** | **PASS** (inside [3,300], real and non-trivial) |
| Aβ40 slope, 3.5-10µM (straddles KM crossover √31≈5.6µM) | -1.010 (R²=0.997) | between the -1.5/-0.5 asymptotes, as expected |
| Aβ40 slope, 20-70µM (saturated regime) | -0.680 (R²=0.99997) | approaching the -0.5 saturated limit |
| Aβ40 single power-law fit, full 3.5-70µM range | -0.798, **R²=0.9909** | — |
| Aβ42 single power-law fit, its own 1-8µM range | -1.320, **R²=0.99994** | — |
| **R² gap (curvature signature)** | **0.0091** | **PASS** (≥1e-4 threshold): Aβ40 shows real, machine-measurable curvature (concentration-dependent effective order) that Aβ42 does not — the quantitative form of Meisl 2014's own qualitative claim. |

**Decorrelated cross-check (different technique entirely, not ThT kinetics):** Bitan et al. 2003
PNAS (PMID 12506200, verified live), using PICUP photo-crosslinking + size-exclusion
chromatography + dynamic light scattering + circular dichroism + electron microscopy: Aβ40 exists
as a monomer/dimer/trimer/tetramer rapid equilibrium, while Aβ42 preferentially forms distinct
pentamer/hexamer "paranuclei" that assemble into protofibrils — a **structurally** distinct
early-assembly pathway, from a technique fully decorrelated from fluorescence kinetics, agreeing
in the same Aβ42-more-aggregation-prone direction. Jarrett, Berger & Lansbury 1993 Biochemistry
(PMID 8490014, verified live) is the historical anchor establishing nucleation-dependence itself
and that Aβ42/43 (not Aβ40) seeds cross-aggregation — qualitative, no fold-number extracted, but
the earliest evidence for the same direction.

### Honest gaps (Part 3)
- The Aβ40 rate-constant reconstruction is an order-of-magnitude, disclosed placeholder, not the
  paper's exact fitted digits (Fig. 3B table not machine-extractable this session).
- The 4.90x matched-concentration fold-difference is this script's own **prediction**, not an
  independently-cited literature fold-number (the source text describes the direction and
  "order of magnitude smaller rate constants" qualitatively; no single citable "Aβ40 is Nx slower"
  headline number was found/verified this session to check 4.90 against directly).
- pH differs between the two peptides' original datasets (Aβ42 at pH 8.0, Aβ40 at pH 7.4, per
  Meisl 2014's own disclosure) — the paper itself flags that a matched-pH comparison would likely
  show an even larger gap; not modeled here (an honest, disclosed asymmetry in the underlying data
  this reconstruction inherits).

## Part 4 — Decorrelated in-vivo check: CSF Aβ42 and the Aβ42/Aβ40 ratio vs amyloid-PET

This is a **genuinely decorrelated** anchor: a different organism (living human brain, not a test
tube), different technique (immunoassay + PET imaging, not ThT fluorescence), different timescale
(years of deposition, not hours of in-vitro polymerization). It is **not** used to validate Part
1-3's kinetics quantitatively — only the qualitative direction (Aβ42 pathologically altered) and
the biomarker-ratio logic (normalizing for inter-individual production variance) are shared
concepts.

**Claim (C1) — the paradoxical CSF Aβ42 decrease:** Olsson et al. 2016 Lancet Neurology
(PMID 27068280, verified live) — systematic review + meta-analysis, **231 articles, 15,699 AD
patients, 13,018 controls** — CSF Aβ42 average AD/control ratio = **0.56** (95% CI 0.55–0.58,
p<0.0001), i.e. a ~44% concentration decrease, paradoxically *as* Aβ deposits into insoluble
plaque (soluble CSF pool depletes as it sequesters into brain parenchyma). Contrast: CSF T-tau
ratio 2.54, P-tau ratio 1.88 (**both increased**, opposite direction) — same meta-analysis, same
patients, direct evidence the two markers move oppositely. **Pre-registered threshold:** ratio <1
with CI excluding 1 in a large-N meta-analysis → **PASS**, trivially and robustly met.

**Claim (C2) — the Aβ42/Aβ40 ratio outperforms Aβ42 alone against amyloid-PET (the anchor):**
Janelidze et al. 2016 Ann Clin Transl Neurol (PMID 27042676, verified live, Table 1 extracted
directly from full text): against amyloid PET status (n=108 PET-positive + comparison cohort),
Aβ42 alone AUC=**0.894** (cutoff <507.5 pg/mL, sens 0.832, spec 0.833); Aβ40 alone AUC=**0.556**
(near chance); **Aβ42/Aβ40 ratio AUC=0.912** (cutoff <0.16, sens 0.900, spec 0.900) —
**significantly better than Aβ42 alone (p=0.002)**. Independently, Hansson et al. 2018
Alzheimers & Dementia (PMID 29499171, verified live, Table 2): Aβ42 alone vs PET AUC=86.5%; a
tau/Aβ42 ratio (a different, but conceptually analogous, normalizing ratio) reaches AUC=94.4-96.3%
— same qualitative "a ratio beats the raw concentration" finding, independently replicated with a
different ratio and a different cohort (BioFINDER/ADNI). **Pre-registered threshold:** ratio-AUC
significantly (p<0.05) exceeds Aβ42-alone AUC in a real, PET-anchored cohort → **PASS** (p=0.002).

### Symmetric QC — what Part 4 does NOT prove
- **In-vitro kinetics ≠ in-vivo deposition timescale.** Parts 1-3 model hours; Part 4 reflects
  years-to-decades of accumulation. No claim is made that the master equation's rate constants
  govern the in-vivo timescale — held explicitly OPEN.
- **In-vitro aggregation is highly prep/seed/agitation-sensitive** (Cohen 2013's own Fig. 2 shows
  the scaling exponent itself drifts from -1.33 (quiescent) toward the fragmentation limit -0.5
  under shear) — reproducibility of any single lab's absolute rate constants across labs/batches
  is **not** independently tested here; held OPEN per the task's own brief.
- Amyloid-PET itself is an imperfect proxy for true plaque burden (tracer-specific, threshold-
  dependent) — treated here as a real, independent, but not infallible, in-vivo anchor.
- The CSF Aβ42 decrease and the ratio's AUC advantage are **correlational clinical biomarker
  facts**, not evidence that the master-equation *mechanism* (secondary nucleation) is what
  drives in-vivo plaque growth — that inferential leap is not made.

## Symmetric QC — overall claimed vs. not claimed

**Claimed (PASS, machine-verified, externally anchored):**
1. A master equation derived from first principles (linearized-Jacobian eigenvalue) reproduces —
   using the source paper's own fitted combined rate constants, untuned — the **measured** Aβ42
   ThT half-time-vs-concentration scaling exponent (simulated -1.320 vs measured -1.33±0.03).
2. The forced adversary (Cohen 2013's own no-secondary-nucleation alternative fit) is
   distinguishable via a machine-computed curve-shape metric (57.6% relative gap in normalized
   rise time), matching the source paper's own qualitative "cannot describe the full time course"
   verdict in machine-checkable form.
3. Aβ42 aggregates measurably faster than Aβ40 (4.90x at matched concentration, inside a wide
   pre-registered band), Aβ40 shows a real machine-measured curvature signature (R² gap 0.0091)
   Aβ42 lacks, and an entirely decorrelated structural technique (Bitan 2003, PICUP/SEC/DLS/CD/EM)
   independently confirms the same directional propensity difference via distinct oligomeric
   pathways.
4. CSF Aβ42 is robustly decreased in AD (0.56x, N=28,717, meta-analytic) while tau markers
   increase (opposite direction, same meta-analysis); the Aβ42/Aβ40 ratio significantly
   outperforms Aβ42 alone against amyloid-PET as an external, decorrelated anchor (AUC 0.912 vs
   0.894, p=0.002), independently corroborated by a second cohort/ratio pairing (Hansson 2018).

**NOT claimed (explicitly OPEN):**
1. No claim that in-vitro kinetics quantitatively predicts in-vivo deposition rate/timescale.
2. No claim of cross-lab reproducibility of absolute rate constants (prep/seed/agitation
   sensitivity is a known, disclosed confound, not tested here).
3. Aβ40's exact fitted rate-constant table (Fig. 3B) was not machine-extracted (rendered graphic);
   the reconstruction is order-of-magnitude, disclosed as such throughout.
4. The slope-alone adversary check does not conceptually rule out a differently-tuned adversary
   matching the exact -1.33 value — only the shape metric is treated as the non-gameable
   discriminator, consistent with the source paper's own reasoning.
5. Amyloid-PET is treated as a real but imperfect anchor, not ground truth.
6. No claim about disease-CAUSALITY (amyloid-hypothesis-vs-tau-vs-inflammation priority) — that
   adjudication lives in `NEU-ALZHEIMERS-MULTIMODAL` and is not re-opened here.

## Confidence tier

**In-vitro-anchored (ThT kinetics, Parts 1-3) + in-vivo biomarker (CSF/PET, Part 4) — two
tiers, explicitly not fused into one number.** Part 1-3's mechanism is anchored to real ThT
kinetic data (Cohen 2013, Meisl 2014) with a machine-reproduced exponent within the measured
error bar; Part 4 is anchored to real clinical cohorts against amyloid-PET. The two tiers are
**decorrelated by design** (different technique, organism, timescale) and are reported
side-by-side, never combined into a false single confidence score.

## `couples_to`
- `NEU-ALZHEIMERS-MULTIMODAL` (OPEN) — the disease-causality triangulation cell (RCT effect
  sizes, Braak/plaque-cognition, APOE4/TREM2). This doc supplies the quantitative aggregation-
  kinetics forward model that node's own prose gestures at but does not build; no overlap
  confirmed by direct read (that node holds zero rate-constant/scaling-exponent numbers).
- `NEU-TAU-PROPAGATION` (OPEN) — the explicit downstream partner: that node's own La Joie 2020
  leg (PMID 31894103) found tau-PET, not amyloid-PET, predicted 15-month structural atrophy in
  its cohort — consistent with this doc's Part 4/Jack 2013 framing of Aβ as the earlier-changing,
  upstream biomarker and tau as the more proximate driver of neurodegeneration (Jack et al. 2013
  Lancet Neurol, PMID 23332364, verified live: "the sequence of events...is Aβ pathophysiology
  first, then tau related neurodegeneration" — with that same paper's own disclosed caveat that
  autopsy evidence in young individuals shows tauopathy can precede Aβ deposition; a population-
  level tendency, not a universal per-individual law).
- `NEU-GLIA-NEUROINFLAMMATION` (OPEN) — its glymphatic-clearance sub-cell already gates
  "amyloid/metabolite clearance rate by sleep-stage architecture"; this doc's Part 4 CSF-Aβ42-
  decrease finding is the clearance/deposition-balance observable that sub-cell's clearance-rate
  parameter would feed into (not built or edited here).
- **Aging** (domain-level coupling, not a single tight node this session): age is the dominant
  epidemiological risk factor for the whole cascade this doc models; `MOL-CELLULAR-SENESCENCE-
  SASP` (OPEN) is this repo's nearest existing aging/senescence axis, but the link here is
  domain-level (both are age-associated proteinopathy/degeneration processes), not a mechanistic
  coupling verified this session — disclosed as such, not overstated.

## Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python scripts/msk/amyloid_beta_aggregation.py   # writes reports/probes/amyloid_beta_aggregation_results.json
```
Pure numpy/scipy (`solve_ivp`) ODE integration — no OpenSim, no external data download. Runtime:
under 30 seconds. No git operations; no writes outside
`reports/probes/amyloid_beta_aggregation_results.json` and this doc pair (this `.md` +
`docs/MECHANISM_AMYLOID_BETA_AGGREGATION_evidence.json`).
