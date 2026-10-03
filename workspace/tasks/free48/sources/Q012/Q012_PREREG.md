# BT-HX-Q012 — preregistration (frozen before first run)

## 0. Scope

The question from BodyTwin's functionality space (`inputs/QUESTION.md`):

> **Q012** Which diet and microbiome variables are needed for a specific substance question?
> Must be able to represent: *metabolite and exposure outcomes with remaining uncertainty.*
> Entry point in the material: **K03, K11, K14**.

**K03/K11/K14 lack labels in this workspace.** Searching `../../notes/`, `../../data/`,
`../../docs/`, `../../scripts/`, `../../START.md`, `../../STARTUP_MESSAGE.md` and `../../tasks/`
for `K03|K11|K14` gives zero hits; no `results/K03|K11|K14` exist (3105 directories in `results/`
searched without a hit). The labels are therefore mapped **positionally** according to `inputs/QUESTION.md`'s own
column "What it must be able to represent" — "**Diet**, **substance/metabolites** and comparable
**microbiome function**" ↔ K03, K11, K14. The content of the K nodes themselves is `UNKNOWN`; see §9.

This is a **first executable mechanistic model**, not an overview. No measured data, no internal
data, no personal information. Only published literature.

## 1. The specific substance question

For Q012 to be decidable, the substance question must be locked. Here:

> **S1.** Which dietary dose (precursor intake) and which microbial functional capacity
> (`y_microbial_conversion`) determine the steady-state plasma concentration of
> **trimethylamine-N-oxide (TMAO)**, and how large is the remaining uncertainty in this concentration?

The choice is motivated by TMAO specifically having (a) a stoichiometric, stepwise chain that can be written
from first principles, (b) quantitative human interventions with isotope labeling, and (c) a
quantified clinical anchor. See §3.

## 2. Hypotes (mekanistisk, falsifierbar)

> **H1 (microbial dependence).** Steady-state plasma concentration of TMAO is linearly proportional
> to the product `f_colonic_escape · y_microbial_conversion` at hepatic load far below capacity.
> Consequently, **Φ = 1 − C_ss(y→0)/C_ss(y_ref)**, the fraction of steady-state TMAO that requires
> bacterial conversion, is **≥ 0.95**.

> **H2 (identifiability).** From a single steady-state plasma concentration without isotope labeling,
> **only the product `p = f_colonic_escape · y_microbial_conversion` is identifiable**; the dietary factor
> and microbiome factor are separable only if precursor passage is also observed
> (labeled test) or TMA flow separately. The model must show this numerically, not assert it.

> **H3 (capacity ceiling).** Above hepatic FMO3 capacity (`f_fmo3·Vmax`), C_ss becomes independent of
> dietary intake. It is the only mechanism that can make the dietary variable irrelevant to this particular
> substance question, and its threshold must be calculated.

H1–H3 are preregistered here, before running. No parameter may be chosen after a result has been
seen.

## 3. Reference values (looked up, primary sources)

### R1 — Koeth RA, Wang Z, Levison BS, Buffa JA, Org E, Sheehy BT, Britt EB, Fu X, Wu Y, Li L,
Smith JD, DiDonato JA, Chen J, Li H, Wu GD, Lewis JD, Warrier M, Brown JM, Krauss RM, Tang WH,
Bushman FD, Lusis AJ, Hazen SL. "Intestinal microbiota metabolism of L-carnitine, a nutrient in
red meat, promotes atherosclerosis." *Nature Medicine* 2013;19(5):576–585.
**DOI 10.1038/nm.3145**, PMID 23563705, PMC3650111. Fulltext verifierad via Europe PMC
(2026-09-25).

- **Fig. 1 / Supplementary Methods:** L-carnitine test = 250 mg d₃-L-carnitine in a capsule **+**
  227 g (8 oz) sirloin steak, corresponding to an estimated **180 mg L-carnitine**. Unit: mg, per test.
- **Fig. 1b–e:** after **1 week** of oral **metronidazole 500 mg 2×/day + ciprofloxacin 500 mg
  2×/day**: *"near complete suppression of endogenous TMAO in both plasma and urine"* and
  *"virtually no detectable formation of either native or d3-labeled TMAO"*. Verbal, no number.
- **Fig. 2c:** fasting baseline plasma TMAO, **n = 26** vegans/vegetarians versus **n = 51** omnivores,
  significantly lower in vegans/vegetarians.
- **Fig. 3a:** after 15 weeks of L-carnitine supplementation, d₃-TMA/d₃-TMAO formation is induced
  **≈10-fold** in mice.
- **Fig. 5d:** reverse cholesterol transport (RCT) dropped **≈30 %** (choline and carnitine diets respectively)
  and **35 %** (TMAO diet) versus normal chow.
- **Fig. 5a,b:** dietary L-carnitine approximately doubled aortic root plaque burden (≈2×).

### R2 — Tang WH, Wang Z, Levison BS, Koeth RA, Britt EB, Fu X, Wu Y, Hazen SL. "Intestinal
microbial metabolism of phosphatidylcholine and cardiovascular risk." *New England Journal of
Medicine* 2013;368(17):1575–1584. **DOI 10.1056/NEJMoa1109400**, PMID 23614584. Abstrakt
verifierat via Europe PMC (2026-09-25).

- **Methods:** phosphatidylcholine challenge = **two hard-boiled eggs** + d₉-labeled phosphatidylcholine;
  plasma **and urine TMAO** were quantified before and after oral broad-spectrum antibiotics.
- **Results:** *"Plasma levels of TMAO were markedly suppressed after the administration of
  antibiotics and then reappeared after withdrawal of antibiotics."* Verbal, no number.
- **Results (PRIMARY NUMERICAL VALUE):** hazard ratio for MACE, highest versus lowest TMAO quartile over
  3 years of follow-up, n = 4007 (selective coronary angiography):
  **HR = 2.54, 95 % CI 1.96–3.28, P < 0.001.** Unit: dimensionless hazard ratio.

### R3 — Wang Z, Klipfell E, Bennett BJ, Koeth R, Levison BS, Dugar B, Feldstein AE, Britt EB,
Fu X, Chung YM, Wu Y, Schauer P, Smith JD, Allayee H, Tang WH, DiDonato JA, Lusis AJ, Hazen SL.
"Gut flora metabolism of phosphatidylcholine promotes cardiovascular disease." *Nature*
2011;472(7341):57–63. **DOI 10.1038/nature09922**, PMID 21475195, PMC3086762. Abstrakt
verifierat via Europe PMC (2026-09-25).

- Choline, TMAO and betaine were identified as metabolites of dietary phosphatidylcholine that predict
  CVD risk; germ-free mice confirm that **both dietary choline and gut flora** are required for TMAO formation;
  genetic variation in **flavin monooxygenase** expression segregated with atherosclerosis → the liver is
  the TMA→TMAO step.

### What could NOT be looked up

Absolute levels for plasma TMAO (µM), renal excretion ratios for TMAO and a **numerical**
suppression ratio for the antibiotic existence are **not in the text I could fetch**
(NCBI abstracts and PMC3650111 full text; PMC3086762 and PMC3701945 returned HTTP 500/404 from
Europe PMC REST, supplement PDFs not fetched). This is marked **UNVERIFIED** and is
not used as a number — only as a direction. No measurements are invented to fill the gap.

## 4. Frozen acceptance review metric

Calculated on `model.py` with the frozen values in §5. No retuning afterward.

- **G1 (microbial dependence index, primary).** Φ ≥ **0.95**.
  The reference (R1 Fig. 1b–e) is verbal ("near complete suppression", "virtually no detectable
  formation"); the threshold 0.95 is a **declared interpretation** of the verbal expression, not a
  measured number. The condition is thus conditional on `r_host_umol_h` (§5) and this must be stated explicitly.
- **G2 (exposure contrast, secondary).** The model must produce a
  concentration contrast between low- and high-exposure modes of **≥ 2×**. Rationale: R2's
  HR = 2.54 (1.96–3.28) for highest versus lowest quartile requires an exposure contrast across
  the quartile range to exist and have an unambiguous direction. **HR and concentration ratio are different quantities;
  the model does not identify HR from concentration, and G2 only tests existence and direction.**
- **G3 (identifiability, H2).** The attribution run with `f_colonic_escape` and
  `y_microbial_conversion` must give the product `p` with **relative 90 % interval width < 0.10** while
  each factor alone has **> 0.5**. If any part becomes < 0.10, H2 is falsified in that run.
- **G4 (capacity ceiling, H3).** The hepatic capacity rate `C_kap = (f_fmo3·Vmax + r_host)/CL`
  must exceed the predicted steady-state concentration under the frozen protocol by
  **≥ 3×**; then the diet is linear and H3 is triggered only outside the protocol. If the ratio < 1,
  G4 FAIL and H3 already applies under the protocol.

## 5. Frozen protocol and parameter values

Reference individual: 70 kg. Precursor: L-carnitine (free base, MW 161.2 g/mol).
Diet block (K03), microbiome block (K14), host block (K11).

| Symbol | Value | Unit | Source / assumption |
|---|---:|---|---|
| `body_mass_kg` | 70 | kg | assumption: adult reference individual |
| `precursor_intake_g_day` | 0.430 | g/d | R1 Fig. 1 + Supp. Methods: 250 mg + 180 mg (source) |
| `precursor_mw_g_mol` | 161.2 | g/mol | chemical formula weight of L-carnitine, free base (source) |
| `f_colonic_escape` | 0.30 | – | **assumption**, not looked up |
| `nu_stoich` | 1.0 | mol/mol | stoichiometry of CntA/CntB and CutC respectively, 1:1 (assumption from reaction stoichiometry) |
| `y_microbial_conversion` | 0.20 | – | **assumption**, not looked up |
| `f_tma_to_hepatic` | 0.90 | – | assumption: fraction of TMA reaching the liver |
| `f_fmo3` | 0.90 | – | assumption: fraction of TMA taking the FMO3 pathway (R3 supports FMO3 as a TMAO source) |
| `fmo3_vmax_umol_h` | 1000 | µmol/h | assumption: hepatic capacity, not looked up |
| `fmo3_km_umol_l` | 5.0 | µmol/L | assumption, not looked up |
| `gfr_ml_min_per_kg` | 1.8 | mL/min/kg | assumption: adult reference value (not looked up in this run) |
| `f_reabsorption` | 0.0 | – | assumption |
| `f_tubular_secretion` | 0.10 | – | assumption |
| `vd_l` | 49.0 | L | assumption: 0.7 L/kg |
| `r_host_umol_h` | 0.062 | µmol/h | **assumption**: non-microbial TMAO source, calibrated so that Φ≈0.98; **not looked up** |
| `tmao_mw_g_mol` | 75.11 | g/mol | chemical formula weight of TMAO (source) |
| `mc_draws` | 20000 | – | frozen; `mc_seed = 20260925` |
| `hepatic_mode` | `linear` | – | frozen for the main run |

Priors for uncertainty propagation (lognormal, median/σ_ln, **all declared assumptions**):

| Parameter | median | σ_ln | Block |
|---|---:|---:|---|
| `f_colonic_escape` | 0.30 | 0.60 | diet |
| `y_microbial_conversion` | 0.20 | 0.90 | microbiome |
| `f_tma_to_hepatic` | 0.90 | 0.20 | microbiome |
| `f_fmo3` | 0.90 | 0.20 | host |
| `gfr_ml_min_per_kg` | 1.80 | 0.15 | host |
| `vd_l` | 49.0 | 0.20 | host |
| `r_host_umol_h` | 0.062 | 1.00 | host |

Sensitivity: OAT ±50 % on exactly three **frozen preselected** parameters —
`y_microbial_conversion`, `f_colonic_escape`, `cl_renal_l_h`. No parameter is chosen after results.

## 6. Null model and placebo

- **Null model A (diet-proportional, "diet is enough"):** `C_ss = k·I` with a single calibrated
  coefficient, no microbial term, no capacity mode. It may pass G2 but **must miss G1
  if H1 applies**; otherwise H1 is falsified.
- **Null model B (microbial fraction, "the microbiome handles everything"):** microbial term, but
  `f_colonic_escape` set to 1.0 and `f_fmo3` to 1.0. Must overpredict C_ss.
- **Placebo:** the same code with `hepatic_mode="linear"` and `fmo3_km_umol_l = 1e12` (non-saturating
  FMO3) must reproduce the closed expression for `linear` to numerical tolerance. A placebo that
  is not reproduced is FAIL.
- **Countertest of placebo identifiability:** vary only `y_microbial_conversion` and check
  that Φ does not change (Φ is independent of y by definition). If Φ moves, the code is wrong.

## 7. What counts as an error

1. Fitting any parameter to R1/R2's qualitative directions after running.
2. Reporting Φ ≥ 0.95 without also stating that `r_host_umol_h` is a free assumption.
3. Using R2's HR = 2.54 as though it were a concentration ratio.
4. Claiming a numerical antibiotic suppression ratio (it is UNVERIFIED).
5. Filling a missing absolute level with an invented plasma TMAO value.
6. Reporting a Monte Carlo run with a different number of draws or seed than §5.

## 8. What counts as support (non-requirement)

That G1 and G4 both pass, that G3 shows weak identifiability of individual factors, and that
variance decomposition ranks the microbiome parameter above the diet parameter. These are reported
as they turn out; they are not conditions for `PASS`.

## 9. Builds on / not redone

- **Builds on:** no previous result files in `results/BT-HX-Q012/` (at
  start the directory contained only `BRIEF.md`, empty `ALLOW_WEB`, `inputs/NIGHT_PREAMBLE.md`, `inputs/QUESTION.md`).
  No previous model, no calibration from `results/` (3105 directories list only
  ID services, no `K03|K11|K14`). The form follows the house `model.py`/`PREREG.md`/`RESULTS.md`
  (e.g. `results/BT-HX-Q100`), **structure only** — no code or parameter imported.
- **Not redone:** BT-B24/N12 null models and the reference model baseline are **irrelevant here** and
  deliberately not acquired; they concern joint force, not metabolite kinetics. No cloud jobs, no
  LOSO sweeps, no `*` import, no `pgrep -f`, no `nvidia-smi -q`.
- **K node labels:** `UNKNOWN`, see §0.
