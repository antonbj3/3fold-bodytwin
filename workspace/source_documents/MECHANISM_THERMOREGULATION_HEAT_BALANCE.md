# MECHANISM THERMOREGULATION — HEAT-BALANCE (M−W=R+C+K+E) EXTENSION (2026-07-22)

**This is a companion/delta doc, not a replacement.** `docs/MECHANISM_THERMOREGULATION.md` +
`scripts/msk/thermoregulation.py` (2026-07-21) already exist, are real, and were **re-run fresh
this session — every headline number reproduced byte-for-byte** (735.7 / 593.5 / 352.5 g/h required
sweat rate, all 8 gates PASS). That script answers a related but distinct question (adiabatic
dT/dt upper bound + an independent Saltin-Hermansen empirical core-temperature equilibrium). This
doc answers the specific, literal brief given for this task: the named steady-state equation
**M − W = R + C + K + E**, DuBois body-surface-area, the 1-MET/58.2-W/m² anchor, explicit core
(37°C)/skin (33°C) reference temperatures, and two decorrelated falsifiers stated in Watts. Per
this repo's isolation rule ("touch only files you create," another instance writes concurrently),
the pre-existing doc/script/JSON are **read-only inputs here, never edited**. Script:
`scripts/msk/thermoregulation_heat_balance.py`. Evidence:
`data/msk_smoketest/subject2_walking1/thermoregulation_heat_balance/thermoregulation_heat_balance_results.json`.

## The governing equation, and how it unifies the sibling script's own two calculations

Classical human heat balance (ASHRAE Fundamentals / ISO 7933 / Gagge two-node lineage — the same
lineage Malchaire 2006, already cited by the sibling script, belongs to):

```
S = M − W − R − C − K − E        (S = rate of body heat storage, W)
```

This is a **conservation law** (first law of thermodynamics on the body as a control volume), not
an empirical fit — nothing to "verify" about the equation itself. What requires external anchoring
is whether the computed **M** is a correct metabolic rate, and whether the required **E** is
achievable against a real measured max sweat-rate capacity. Two limits of this ONE equation,
previously computed by the sibling script as if they were separate calculations, are unified here
for the first time:

| limit | condition | what it gives | source |
|---|---|---|---|
| (a) zero dissipation | R=C=K=E=0 | all of M−W stored → dT/dt bound | sibling's own "whole-body dT/dt adiabatic upper bound," **reused unmodified** |
| (b) steady state | S=0, R+C+K held flat at its resting value (disclosed upper-bound-on-E simplification) | required E = exercise-increment heat | sibling's own "required_sweat_rate_g_h," **independently reimplemented + cross-checked exactly (rel_diff=0.0, all 3 configs)** in Step 5 below |

R/C/K's **individual** split is deliberately, explicitly held **OPEN** here — per this task's own
pre-registration ("R/C/K coefficients depend on environment/clothing/posture — hold open") and
because two live-fetch attempts on Malchaire (2006)'s own full-text convective/radiative equations
both failed to yield readable equations this session (disclosed in full below, not hidden). Only
the **lumped** (R+C+K) dry-heat term, which conservation alone determines once M, W, E are fixed,
is reported — as a rigorous bound, not an invented 3-way split.

## Citations — verified LIVE this session (WebSearch budget was exhausted mid-session — a shared,
session-wide cap; fell back to direct WebFetch against pubmed.ncbi.nlm.nih.gov / NCBI eutils /
Wikipedia, the same fallback `docs/MECHANISM_METABOLIC_CROSS_ACTIVITY.md` already used for the same reason)

| source | verified as | what it anchors | this session's OWN catches |
|---|---|---|---|
| Du Bois D, Du Bois EF (1916, reprinted *Nutrition* 1989;5(5):303-11) "A formula to estimate the approximate surface area if height and weight be known." **PMID 2520314** | Bibliographic (title/authors/journal/year) verified live; PubMed record itself carries no abstract (pre-abstracting-era reprint) | Existence/identity of the classic BSA reference | — |
| DuBois BSA coefficient: `BSA(cm²) = 71.84 × mass(kg)^0.425 × height(cm)^0.725` | **TEXTBOOK-GRADE** (Wikipedia "Body surface area," live-fetched) — not independently re-derived from Du Bois's own 1916 primary text this session, flagged | The exact formula used below | Sanity-checked: a reference 70kg/170cm adult gives BSA=1.81 m², matching the standard "~1.8 m²" figure independently |
| Jetté M, Sidney K, Blümchen G (1990) "Metabolic equivalents (METS) in exercise testing..." *Clin Cardiol* 13(8):555-65. **PMID 2204507** | **Full abstract fetched & read live**: "One metabolic equivalent (MET) is defined as the amount of oxygen consumed while sitting at rest and is equal to 3.5 ml O2 per kg body weight × min." | 1 MET = 3.5 mL O2/kg/min | An earlier guessed PMID (2225418) for this same citation resolved to an **unrelated 1990 cardiology paper** (myocardial actuation during ventricular fibrillation, dogs) — caught via live fetch, corrected via a fresh NCBI esearch author-query to the true PMID, not silently kept |
| 1 MET = 1.162 W/kg = **58.2 W/m²** | **TEXTBOOK-GRADE** (Wikipedia "Metabolic equivalent," live-fetched) — same conversion `metabolic_cost.py`'s own docstring already uses/flags | The task brief's own "≈58 W/m², ≈100 W for 1.8 m²" figure — **confirmed exactly**: 58.2×1.8=104.76 W | — |
| **Mifflin MD, St Jeor ST, Hill LA, Scott BJ, Daugherty SA, Koh YO (1990)** "A new predictive equation for resting energy expenditure in healthy individuals." *Am J Clin Nutr* 51(2):241-7. **PMID 2305711** | **Full abstract fetched & read live**: n=498 healthy adults (247F/251M, age 19-78y), REE **measured by indirect calorimetry**, regression R²=0.71: `REE(males)=10×weight(kg)+6.25×height(cm)−5×age(y)+5` | **PRIMARY external anchor for Falsifier 1** — a real, decorrelated (different population, different measurement modality — whole-body indirect calorimetry, not a model-internal constant), directly-measured-RMR regression | — |
| Wikipedia "Human body temperature" (live) | "generally given as 36.5–37.5 °C"; rectal/core "37.0 °C" | CORE = 37 °C reference (task brief's own value, confirmed) | — |
| Wikipedia "Skin temperature" (live) | "Normal human skin temperature on the trunk of the body varies between 33.5 and 36.9 °C" | SKIN = 33 °C reference (task brief's own value sits at ~the measured lower bound of a real range) | — |
| Latent heat of sweat, 2.43 kJ/g (task brief's own value) | **REUSED**, not re-fetched: essentially identical to the sibling script's own already-verified 2426 J/g (ISO 7933/Gagge lineage, textbook-grade, flagged there) | 0.16% difference is rounding, not conflict; both used explicitly, never conflated | — |
| Malchaire JBM (2006) PMID 16922181 — max non-acclimatized sweat rate 650-1,000 g/h; observed lab sweat rate 424±172 g/h (n=672) | **REUSED bibliographically** from the sibling script's own JSON (not re-fetched) | Falsifier-2 evaporative-capacity anchor | Two fresh attempts this session to independently pull Malchaire's own convective/radiative (R/C) linearized-coefficient *equations* both **failed** (PDF returned as unreadable compressed binary; the HTML abstract page carries no equations) — disclosed, not silently patched with a recalled number; this is exactly why R/C/K's individual split stays OPEN below, not asserted from memory |
| Wikipedia "Perspiration" (already live-fetched by the sibling script, reused) | "maximum sweat rates of an adult... up to 2-4 litres per hour"; "average intensity exercise... up to 2 litres... per hour" | Falsifier-2 absolute-ceiling anchor | — |

Recall discipline confirmed again this session: **two of my own guessed PMIDs drifted** (Jetté et
al. and a discarded "IUPS thermal glossary" guess that resolved to an unrelated veterinary
anaesthesia paper, PMID 11599678 — dropped entirely rather than kept) — consistent with this
project's own previously-measured ~62% recall-citation-drift rate (`docs/MECHANISM_METABOLIC_COST.md`).
Every number above was independently re-verified live, not recalled.

## Headline: subject2 (mass 78.2 kg, height 1.96 m, male — from `sessionMetadata.yaml`, external
read-only mount; mass cross-checked exactly against `metabolic_cost_results.json`'s own value)

**DuBois BSA = 71.84 × 78.2^0.425 × 196^0.725 = 2.1031 m²** — 116.8% of the task brief's own
"average person" 1.8 m² reference, because this subject's height (1.96 m) is well above an average
adult male: the geometrically expected consequence of the formula's own height^0.725 scaling, not
a bug (sanity-gated: PASS, 1.3–2.8 m² adult band).

### M via four routes (all W, total body)

| route | value | vs Mifflin-St Jeor range [85.3, 92.6] W | Falsifier-1 gate (±20%) |
|---|---:|---:|:---:|
| twin's own OpenSim-internal basal (`M_rest`, reused from `metabolic_cost_results.json`) | **93.84 W** | 105% of mid | **PASS** |
| 1 MET × task-brief generic BSA (1.8 m²) | 104.76 W | 118% of mid | PASS |
| 1 MET × subject2's **own** DuBois BSA (2.103 m²) | **122.40 W** | **138% of mid** | **FAIL** |
| 1 MET × subject2's own mass × 1.162 W/kg | 90.87 W | 102% of mid | PASS |

## FALSIFIER 1 (resting): does modeled resting heat match the metabolic thread's measured RMR? — PASS, with a forced, disclosed adversary

**PRIMARY, pre-registered comparator**: twin's own `M_rest` (93.84 W, from `metabolic_cost.py`'s
own already-published basal constant, reused not recomputed) vs. **Mifflin-St Jeor's real,
indirect-calorimetry-measured regression** (PMID 2305711, n=498) evaluated at subject2's own actual
mass/height across the full plausible healthy-adult age range (age not recorded in
`sessionMetadata.yaml` — disclosed gap): **85.33–92.59 W**, an **age-invariant** range (only 8.5%
wide across 20–50y, so the missing age does not decide this comparison). 93.84 W sits **within
±20% of this real, decorrelated, indirect-calorimetry-measured anchor (105% of its midpoint) →
PASS.** This directly answers the task's own falsifier ("does modeled resting heat production
match the metabolic thread's MEASURED resting metabolic rate, indirect calorimetry") with a
stronger anchor than the task's own generic "~1 kcal/min ≈ 70-80 W" rule of thumb — that generic
figure is a **population-average** approximation; Mifflin-St Jeor's subject-specific prediction
(which properly accounts for this subject's unusually large body size) is the correct,
apples-to-apples comparator, and it is **higher** than 70-80 W for exactly the expected reason
(a materially larger-than-average person has a materially larger-than-average absolute RMR even at
a similar per-kg rate).

**Forced adversary (not swept under the rug)**: the more "rigorous-looking" route — correcting the
generic 58.2 W/m² MET shortcut to this subject's **own** DuBois BSA (122.40 W) — **FAILS** the same
gate, overshooting Mifflin-St Jeor by 38%. Properly personalizing the BSA makes the estimate
**worse**, not better. The generic-BSA route's coincidental proximity to Mifflin-St Jeor (104.76 W,
118% of mid, a PASS) is not evidence the BSA-scaling method is right for this subject — it is
evidence that DuBois BSA scaling (mass^0.425 × height^0.725) does not track fat-free/metabolically-
active mass 1:1 for a height outlier (1.96 m). **Reported as an open, unresolved discrepancy**, not
hidden behind the flattering aggregate.

## Reimplementation cross-check (Step 5): an independent codebase agrees exactly

Before trusting a NEW formula on a NEW (higher-intensity) input, this script's own
`required_evaporative_w` was run on the **sibling script's own three already-published configs**
(Umberger-primary 733.6 W, Bhargava-primary 609.9 W, combined-corrected 400.4 W) and compared
against the sibling's own already-published `required_sweat_rate_g_h` numbers:

| config | this script | sibling (2026-07-21) | rel. diff |
|---|---:|---:|---:|
| umberger_primary | 735.7190 g/h | 735.7190 g/h | **0.0** |
| bhargava_primary | 593.4575 g/h | 593.4575 g/h | **0.0** |
| combined_corrected | 352.5089 g/h | 352.5089 g/h | **0.0** |

**Exact agreement, two independently-written codebases, same physics.** This is the machine
cross-check that licenses reusing the same formula on a genuinely new, higher-intensity input below
— confidence in the *code*, decorrelated from confidence in the input.

## FALSIFIER 2 (sustained exercise): does required E fall within measured max sweat capacity? — CLOSES, with a precise, non-trivial margin

The task's own target band is **500-1000 W** sustained exercise heat; the sibling script's own
walking-based numbers (331–590 W across its three configs) sit at/below this band's low end.
Rather than run a new simulation, this doc reuses the twin's **own already-published, different-
activity** number: squat (Umberger2010, `docs/MECHANISM_METABOLIC_CROSS_ACTIVITY.md`,
`metabolic_cross_activity_combined.json`, gross rate 13.842 W/kg × 78.2 kg = **1082.5 W total**) —
genuinely decorrelated from walking (different SO solve, different movement, same subject/mass).

| η convention | required E | = sweat rate | verdict vs task's own 675-1350 W (1-2 L/h) band | vs Malchaire "typical" non-acclim. max (650-1000 g/h) |
|---|---:|---:|:---:|:---:|
| η=0.225 (conventional) | 766.2 W | 1135.1 g/h (1.135 L/h) | **inside** | **exceeds** |
| η=0 (rigorous, no net external work — squat starts/ends standing, same logic as level walking) | 988.6 W | 1464.6 g/h (1.465 L/h) | **inside** | **exceeds** |

**Both bookends fall entirely inside the task's own stated 675-1350 W / 1-2 L/h capacity band** —
the decorrelated energy-balance closure the falsifier asked for **holds**. At the same time, both
bookends **exceed** Malchaire's own measured "typical non-acclimatized" sweat-rate ceiling
(650-1,000 g/h) while staying comfortably inside Wikipedia's broader absolute physiological ceiling
(2-4 L/h). **Read precisely, not rounded up**: holding core temperature steady at squat-intensity
sustained exercise is **physiologically possible but non-trivial** — it requires heat-acclimatized
or above-typical sweating capacity, not "any person's" typical output. A genuine, borderline,
falsifiable, machine-checked finding — not a rubber-stamped PASS.

**Disclosed limitation on this falsifier**: squat's own Static-Optimization solution is flagged
**INHERITED-IMPLAUSIBLE** in `docs/MECHANISM_METABOLIC_CROSS_ACTIVITY.md` (muscle-pinning in 4
muscles) — that flag concerns whether the specific 1082.5 W number is the *exact right* squat
metabolic rate, not whether squatting produces vigorous-exercise-range heat at all (a well-
established qualitative fact independent of this specific SO solve's fidelity). A cleaner version
of this falsifier would use a fully externally-sourced (not twin-simulated at all) high-intensity
MET-table value as a second, even-more-decorrelated cross-check — **not done this session**
(the specific numeric MET-table codes were not independently live-verified, and asserting one from
memory here would reintroduce exactly the citation-drift risk this doc is built to avoid) — a
disclosed, not silently patched, gap.

## R+C+K — held OPEN, exactly as pre-registered

At rest (W=0, steady state): **(R+C+K)_rest + E_rest = M_rest = 93.84 W** (exact, conservation).
E_rest (insensible loss — evaporation that occurs without active sweating — is physiologically
distinct from the exercise-driven sweating Falsifiers 1-2 concern) was **not** independently
live-verified this session (the two failed Malchaire full-text fetches above). **(R+C+K)_rest is
therefore reported only as a rigorous upper bound, ≤93.84 W** (since E_rest≥0) — with **no**
individual R-vs-C-vs-K numeric split asserted. This is the task's own explicitly pre-registered
open scope boundary ("R/C/K coefficients depend on environment/clothing/posture — hold OPEN"), kept
open by design, not by omission.

**Core (37.0 °C) − skin (33.0 °C) = 4.0 °C** internal gradient drives the conductive+convective heat
transport from the core to the skin surface *before* any of R/C/K/E can act at the body-environment
boundary — a real transport step this 0-D (no spatial resolution, no perfusion) model does not
resolve, structurally the same gap the sibling script already discloses ("no vascular/perfusion
heat transport... the single biggest structural gap"), now named via the explicit core-skin ΔT
rather than only via the muscle-adiabatic-bound illustration.

## Void-floor / non-degeneracy check (forced, machine-gated, not eyeballed)

Sweeping M∈[200,1200] W (n=25, this script's own new `required_evaporative_w` formula, independent
of the sibling's own already-published sweep): required-E range **[82.3, 857.3] W**, **strictly
increasing** across the whole sweep (PASS — a real function of the input, not a pinned constant),
and the numerical derivative matches the closed-form analytical derivative **dE/dM=(1−η)=0.7750**
exactly (PASS, `rtol=1e-3`).

## Confidence tier — symmetric-QC correction against this repo's own established convention

The task brief's own suggested label is "in-vivo-anchored." Per `docs/MECHANISM_TRUST_LEDGER.md`'s
own tier legend — and the **identical** precedent already set by
`docs/MECHANISM_METABOLIC_CALORIMETRY.md` for this exact same question (metabolic-vs-calorimetry) —
this repo reserves **in-vivo-anchored** specifically for direct in-vivo instrumented/telemetry
force measurement (OrthoLoad implant telemetry is the source for every row in that tier). Mifflin-
St Jeor and Malchaire are real, external, decorrelated, indirect-calorimetry-measured human datasets
— genuinely anchored, but a different modality (whole-body respiratory/regression measurement, not
implant telemetry). Per this repo's own established convention, **not** the task brief's generic
phrasing: **cadaveric-or-published-plausibility** — muscle-energetics-model-limited (inherits
`metabolic_cost.py`'s own disclosed ~1.5-2× spread across its primary configs) and
BSA-scaling-limited (this doc's own forced-adversary finding: DuBois BSA over-predicts for this
subject's specific anthropometry).

## Honest gaps (first-step extension, not final)

1. **R/C/K individual split is not resolved** — by task design, and because two live-fetch attempts
   on Malchaire (2006)'s own explicit convective/radiative equations both failed technically this
   session (disclosed above, not silently patched with a recalled coefficient).
2. **E_rest (insensible loss) is unbounded below** — only the upper bound (R+C+K)_rest≤M_rest is
   reported; no live-verified insensible-water-loss figure was found this session (two Wikipedia
   fetch attempts, both empty of the needed number).
3. **Subject2's age is not recorded** in `sessionMetadata.yaml` — Falsifier 1 sweeps the full
   20-50y healthy-adult range instead (shown to be near-age-invariant, 8.5% wide, so this does not
   change the verdict, but it is a real, disclosed data gap, not a silently-assumed age).
4. **Falsifier 2 reuses squat's own INHERITED-IMPLAUSIBLE-flagged metabolic rate** (muscle-pinning
   in the underlying Static Optimization solve) rather than a fully independent, externally-sourced
   high-intensity MET-table figure — disclosed, not hidden; a decisive follow-up would substitute or
   corroborate with a live-verified MET-table value once one is fetched without drift risk.
5. **0-D lumped, no environment, no perfusion** — inherits every one of the sibling script's own
   already-disclosed structural gaps (no spatial resolution, no clothing/wind/humidity, no
   vascular/perfusion heat transport) unchanged; this doc adds the named governing equation and two
   new anchors, it does not remove any prior limitation.
6. **Single subject (subject2), one body composition, one sex** — same scope boundary as every
   other cert in this family; DuBois BSA and Mifflin-St Jeor are population regressions, applied
   here to n=1.
7. **The ±20% Falsifier-1 tolerance and the "675-1350 W"/"650-1000 g/h" Falsifier-2 bands are the
   task's/literature's own pre-registered numbers** — genuine, generic engineering/physiology
   tolerances, not reverse-tuned from this subject's own answer (confirmed non-vacuous: the
   subject-specific-BSA route independently FAILS the same Falsifier-1 gate at the same tolerance,
   demonstrating real discriminating power, not a rubber stamp).

## Repro

```
cd ~/projects/bodytwin
.venv-msk/bin/python scripts/msk/thermoregulation_heat_balance.py
```
Reads (all read-only, none modified): `/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/sessionMetadata.yaml`
(external mount), `data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json`,
`data/msk_smoketest/metabolic_cross_activity/metabolic_cross_activity_combined.json`,
`data/msk_smoketest/subject2_walking1/thermoregulation/thermoregulation_results.json`. Writes only
`data/msk_smoketest/subject2_walking1/thermoregulation_heat_balance/thermoregulation_heat_balance_results.json`.
Pure Python/numpy, no OpenSim call, runs in under a second. Does **not** edit
`scripts/msk/thermoregulation.py`, `docs/MECHANISM_THERMOREGULATION.md`, or their JSON. No git
operations (isolation respected).

## Files

- `scripts/msk/thermoregulation_heat_balance.py` — new, re-runnable, self-contained script (created
  this session; does not edit any pre-existing file).
- `data/msk_smoketest/subject2_walking1/thermoregulation_heat_balance/thermoregulation_heat_balance_results.json`
  — every number in this doc, machine-written: subject anthropometry, DuBois BSA, all 4 M routes,
  the Mifflin-St Jeor age sweep, both falsifiers' full detail, the reimplementation cross-check, the
  void-floor sweep, and all gate verdicts.
- Read-only, never modified: `docs/MECHANISM_THERMOREGULATION.md`, `scripts/msk/thermoregulation.py`,
  `data/msk_smoketest/subject2_walking1/thermoregulation/thermoregulation_results.json` (the sibling
  this doc extends), `data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json`,
  `data/msk_smoketest/metabolic_cross_activity/metabolic_cross_activity_combined.json`,
  `docs/MECHANISM_TRUST_LEDGER.md` (tier vocabulary), `docs/MECHANISM_METABOLIC_CALORIMETRY.md` (the
  tier-correction precedent this doc follows).
