# MECHANISM RESPIRATORY — ventilation / O2-uptake layer (2026-07-21)

Couples the musculoskeletal twin's metabolic-cost capability (`scripts/msk/metabolic_cost.py`,
`docs/MECHANISM_METABOLIC_COST.md`) to whole-body **gas exchange** — a sibling of
`scripts/msk/thermoregulation.py` (which took the same metabolic-rate input to the *thermal* side;
this takes it to the *respiratory* side). From the twin's own already-computed metabolic rate:
oxygen uptake (VO2) via the caloric equivalent of O2, and required minute ventilation (VE) via the
ventilatory equivalent VE/VO2. Script: `scripts/msk/respiratory.py`. Evidence:
`data/msk_smoketest/subject2_walking1/respiratory/respiratory_results.json`.

**NO re-solve, no OpenSim.** Opens `metabolic_cost_results.json` (plain JSON, read-only) and does
pure arithmetic on numbers already there. Never touches the `.osim` model or any `.sto` file.

## Method, in one paragraph

VO2 [mL/(kg·min)] = P [W/kg] × 60 / E_O2, where E_O2 = 20.9 kJ/L O2 is this task's own specified
caloric equivalent of oxygen (moderate-exercise, carbohydrate-leaning fuel mix) — cross-checked,
not just asserted, by inverting Wikipedia's "Indirect calorimetry" RQ-interpolation formula between
its own quoted carbohydrate (21.13 kJ/L, RQ=1.0) and fat (19.62 kJ/L, RQ=0.7) bounds: 20.9 kJ/L
implies **RQ ≈ 0.954**, physiologically appropriate for moderate, carb-leaning exercise (vs. a
lower, more fat-oxidizing ~0.80–0.85 typical at rest — applying one shared E_O2 to both rest and
exercise is a disclosed, <2%-sized simplification, see Honest gaps). VE [L/min] = VO2 [L/min] ×
(VE/VO2), swept across this task's own given range (25 / 27.5 / 30) rather than collapsed to one
number. Both are applied to **four** metabolic-rate inputs read directly from
`metabolic_cost_results.json` (no new metabolic number invented): **rest** (the OpenSim probes' own
basal coefficient, 1.2 W/kg, cross-model self-consistent to machine precision), and **three walking
configs** — Umberger-primary, Bhargava-primary, and **combined-corrected** (tendon-compliance +
literature-anchored muscle-mass correction — `metabolic_cost.py`'s own disclosed sensitivity
variant, and this task brief's own named **"~400 W gross rate"**: 5.1197 W/kg × 78.2 kg = 400.4 W,
read from the JSON's own `sensitivity.combined_tendon_and_mass_correction_w_per_kg` + basal, not
retyped) — carried forward end-to-end, the same three-input discipline `thermoregulation.py`
already established, rather than silently picking one.

## Citations — verified LIVE this session, not recalled

| source | verified as | what it anchors |
|---|---|---|
| Ainsworth BE et al. (2011) "2011 Compendium of Physical Activities: a second update of codes and MET values," *Med Sci Sports Exerc* 43(8):1575–81. **PMID 21681120**, DOI 10.1249/MSS.0b013e31821ece12 | Title/authors/journal/year/volume/pages/DOI all matched live on pubmed.ncbi.nlm.nih.gov | Source of the MET framework. The Compendium's own official public companion database (pacompendium.com, fetched live) gives: code 17255 "Walking, self-selected speed, indoor track or outdoors, firm surface" = **4.0 METs**; code 17190 "Walking, 2.8–3.4 mph, level, moderate pace" = **3.8 METs**; code 17170 "Walking, 2.5 mph, firm, level" = **3.0 METs**. Via 1 MET = 3.5 mL O2/kg/min (a physiological *definition*, not a separately-cited number — the same constant `metabolic_cost.py` already used for its own basal cross-check), these convert to VO2 = **14.0, 13.3, 10.5 mL/kg/min** — the "self-selected speed" code lands at dead center of the task's 12–15 mL/kg/min anchor. |
| Waters RL, Mulroy S (1999) "The energy expenditure of normal and pathologic gait," *Gait Posture* 9(3):207–31. **PMID 10575082**, DOI 10.1016/s0966-6362(99)00009-0 | Title/authors/journal/year/volume/pages/DOI matched live via TWO independent NCBI passes (PubMed-page fetch + esummary) | The classic, most-cited gait-energetics reference. No abstract/numeric table accessible via this session's fetch route (paywalled Elsevier) — cited for topical grounding, **not** as the source of an independently re-extracted number, flagged as such (same discipline `metabolic_cost.py` applied to Janssen et al. 2000). |
| Browning RC, Kram R (2005) "Energetic cost and preferred speed of walking in obese vs. normal weight women," *Obes Res* 13(5):891–9, **PMID 15919843**, DOI 10.1038/oby.2005.103; Browning RC, Baker EA, Herron JA, Kram R (2006) "Effects of obesity and sex on the energetic cost and preferred speed of walking," *J Appl Physiol* 100(2):390–8, **PMID 16210434**, DOI 10.1152/japplphysiol.00767.2005 | Both verified live via two independent passes, abstracts fetched and read | Source of **preferred/self-selected walking speed ≈ 1.40–1.47 m/s** in healthy adults — used to note the twin's own trial speed (1.065 m/s) is *slower* than typical "preferred" speed, making a match to preferred-speed VO2 ranges a genuine, non-guaranteed finding. |
| Wasserman K, Whipp BJ, Koyl SN, Beaver WL (1973) "Anaerobic threshold and respiratory gas exchange during exercise," *J Appl Physiol* 35(2):236–43. **PMID 4723033**, DOI 10.1152/jappl.1973.35.2.236 | Title/authors/journal/year/volume/pages/DOI matched live | The classical primary reference establishing the ventilatory-equivalent (VE/VO2, VE/VCO2) method. No abstract accessible (pre-abstracting era) — the specific "25–30 at moderate exercise" figure is this task brief's own stated, textbook-consensus value associated with this reference, **not** independently re-extracted from the primary text, flagged as such. |
| ATS/ACCP Statement on Cardiopulmonary Exercise Testing (2003), *Am J Respir Crit Care Med* 167(2):211–77. **PMID 12524257**, DOI 10.1164/rccm.167.2.211 | Title/journal/year/volume/pages/DOI matched live | The standard clinical-guideline reference for CPET normal values incl. ventilatory equivalents; no abstract accessible (guideline document) — cited as the standard-reference pointer, not a specific extracted number. |
| Péronnet F, Massicotte D (1991) "Table of nonprotein respiratory quotient: an update," *Can J Sport Sci* 16(1):23–9. **PMID 1645211** (no DOI on record — pre-dates DOIs, confirmed via NCBI esummary's own empty articleids field, not just an omission) | Title/authors/journal/year/volume/pages matched live twice independently, partial abstract read | Provenance/context for the caloric-equivalent-of-O2 concept (updates the classical 1924 Zuntz/Lusk table). The specific numeric cross-check used here comes from Wikipedia (below), not re-extracted from this paper's own table this session. |
| Wikipedia "Indirect calorimetry" (fetched live) | **TEXTBOOK-GRADE, flagged** (same discipline as `metabolic_cost.py`'s 1-MET↔Watts conversion, `thermoregulation.py`'s specific-heat/latent-heat constants) | "21.13 kilojoules (5.05 kcal)... per litre of oxygen by the oxidation of carbohydrate" (RQ=1.0); "19.62 kilojoules (4.69 kcal), the value for fat" (RQ=0.7) + the linear RQ-interpolation formula — used to machine-check the task's 20.9 kJ/L constant (→ implied RQ=0.954). |
| Wikipedia "Minute ventilation" (fetched live) | **TEXTBOOK-GRADE, flagged** (unsourced within the article itself) | "A normal minute volume while resting is about 5–8 liters per minute"; "light activities... around 12 litres"; "moderate exercise may be between 40 and 60 litres per minute" — used as sanity **bands**, not a precise walking-specific value. |

Recall discipline: this project's own prior finding is a **measured ~62% citation-drift rate** from
memory (`MECHANISM_METABOLIC_COST.md`). Every PMID/DOI above was fetched live this session (NCBI
eutils esearch/esummary, direct PubMed-page fetches, the Compendium's own companion site, Wikipedia)
— not recalled. `WebSearch` itself was unavailable this session (shared per-session budget already
exhausted, same constraint `MECHANISM_METABOLIC_COST_CROSS_ACTIVITY.md` hit) — all verification below
used direct `WebFetch` calls instead, the same NCBI-eutils-first-party discipline `metabolic_cost.py`
already established.

## Headline results (subject2/walking1, mass 78.2 kg; four metabolic-rate inputs)

| quantity | rest (basal) | Umberger-primary | Bhargava-primary | **combined-corrected** |
|---|---:|---:|---:|---:|
| Input metabolic rate | 1.200 W/kg (93.8 W) | 9.381 W/kg (733.6 W) | 7.799 W/kg (609.9 W) | **5.120 W/kg (400.4 W)** |
| **VO2** | **3.445 mL/kg/min** (269 mL/min) | 26.930 mL/kg/min (2106 mL/min) | 22.389 mL/kg/min (1751 mL/min) | **14.698 mL/kg/min** (1149 mL/min) |
| **VE** (VE/VO2=25 / 27.5 / 30) | 6.73 / **7.41** / 8.08 L/min | 52.65 / 57.91 / 63.18 L/min | 43.77 / 48.15 / 52.52 L/min | **28.73 / 31.61 / 34.48 L/min** |

All numbers, gates, and the full sensitivity sweep are on disk in `respiratory_results.json`;
nothing above is hand-computed prose — the script prints and JSON-serializes every one of these
from its own arithmetic (`vo2_ml_per_kg_min()`, `ve_l_per_min()`).

## External anchors — three, decorrelated, and a forced discrimination check

1. **Rest vs the 1-MET definition** (a self-consistency check, not an independent discovery — 1 MET
   is *defined* as 3.5 mL O2/kg/min): twin's own basal → **3.445 mL/kg/min**, ratio to definition
   **0.984** (1.6% low) — inside the ±10% pre-registered tolerance. Mostly validates the E_O2
   choice's reasonableness and the pipeline's own arithmetic, reported plainly as such.
2. **Walking VO2 vs the MET/Compendium anchor** (Ainsworth et al. 2011, PMID 21681120; task's own
   stated **12–15 mL/kg/min** range, independently corroborated here by the Compendium's own
   "self-selected speed" code = 14.0 mL/kg/min): **combined-corrected (14.698 mL/kg/min) PASSES**
   — inside the range, ratio to range-midpoint 1.09. **Both primary configs FAIL** — Umberger
   26.930 mL/kg/min (1.8× the range's own upper bound), Bhargava 22.389 (1.5×).
3. **Ventilation vs Wikipedia's textbook bands**: rest VE (mid, 7.41 L/min) sits inside the quoted
   5–8 L/min resting band. Combined-corrected walking VE (31.61 L/min, mid) sits in the
   light-to-moderate transition zone; **both primary configs' VE (57.91, 48.15 L/min) land inside
   the 40–60 L/min "moderate exercise" band instead** — i.e. they read like a jog, not a walk, by
   ventilation alone, independently of the VO2 comparison above.

**The forced adversary**: is this just a generic linear rescaling that makes *any* plausible
metabolic-rate input "look reasonable," discriminating nothing about the twin's own specific
number? Answer: **no** — the published-VO2-range gate passes for exactly one of the three configs
and fails for the other two, in the theoretically-expected direction (over-prediction, not
under-prediction). This is the **same discrimination pattern** `metabolic_cost.py`'s own
Koelewijn-cost-of-transport anchor and `thermoregulation.py`'s own Saltin-Hermansen-validity-domain
anchor already found — now via a **third**, differently-sourced external literature base (MET
activity classification + textbook ventilation ranges), not a re-run of the same comparison.

**SHARED-INPUT CAVEAT, stated once, applies throughout — this is not three independent
replications.** All three layers (COT, core-temperature, VO2/VE) still judge the **same single**
upstream `metabolic_cost.py` OpenSim-probe output. Per this project's own prior finding
(`bt_memory/a-conservation-law-closure-is-a-tautology-if-all-terms-flow-from-one-instrument-chain...md`,
2026-07-18): apparent multi-signal agreement is only as strong as the number of *physically
independent channels* behind it. Here, the three **external anchors** genuinely are decorrelated
(different studies, different measurement techniques: indirect-calorimetry cost-of-transport,
directly-measured core temperature, MET activity classification) — so their **mutual agreement on
which config is more externally consistent** is real, convergent evidence. But the **twin-internal
number itself is one shared instrument chain** across all three layers; if it has an undiagnosed
bug (beyond the already-disclosed Fmax-muscle-mass/tendon-compliance issue), all three inherit it
identically, not independently. Reported as "3 external anchors agree," never as "3 independent
measurements of the twin."

## Void-floor / non-degeneracy check (forced, not eyeballed)

25-point sweep of input metabolic power across [90, 1200] W total (near-rest to near-maximal
sustained exercise, matching `thermoregulation.py`'s own upper bound): VO2 spans
[3.30, 44.05] mL/kg/min, VE(mid) spans [7.11, 94.74] L/min, both **strictly monotonic** (real
function of the twin's own number, not pinned), and the numerically-estimated derivative matches
the closed-form analytical derivative (`d(VO2)/dP = 60/E_O2 = 2.8708` mL/kg/min per W/kg;
`d(VE)/dP = 6.1737` L/min per W/kg) to machine precision (`rtol=1e-6`) — this conversion chain is
exactly linear by construction, so this check mainly guards against a units/coding bug (e.g. a
stray factor of 60 or 1000), which is exactly the kind of error this class of layer is prone to.

## Scope, honest gaps (first step, not final)

- **0-D, STEADY-STATE ONLY.** No breath-by-breath dynamics: no VO2-kinetics time constant (the real
  ~20–45 s exponential rise at exercise onset, "O2 deficit"), no respiratory-rate/tidal-volume
  decomposition of VE — this script predicts only their *product*.
- **NO dead-space modeling.** VE here is *total* (minute) ventilation, not alveolar ventilation
  (VA = VE − VD·f). Real dead space is a genuine fraction of each breath (~30% of resting tidal
  volume, proportionally less at exercise's larger tidal volumes) — not separated here at all, a
  real structural gap if this layer is ever used to reason about alveolar gas exchange specifically.
- **NO diffusion-limitation modeling.** Assumes ventilation and O2 uptake are always adequately
  coupled (stoichiometry alone suffices) — cannot represent a diffusion-limited state (lung disease,
  extreme altitude, elite-athlete exercise-induced arterial hypoxemia at very high cardiac output).
- **NO altitude / inspired-O2 coupling.** Sea-level, normoxic constants only.
- **GENERIC constants** (E_O2, VE/VO2 range) — not this subject's own measured RER, ventilatory
  efficiency, age, sex, fitness, or lung function. Applying the **same** E_O2=20.9 kJ/L to both rest
  and exercise is itself a disclosed simplification: true resting RQ is typically lower (~0.80–0.85,
  more fat oxidation) than the ~0.95 this exercise-tuned constant implies — a <2%-sized effect on
  the rest-VO2 number at this resolution, not separately corrected for.
- **SHARES THE METABOLIC-COST COMMON-MODE** — explicitly, not in passing (see External anchors
  section above): every number here is a linear re-expression of the same single
  `metabolic_cost.py` OpenSim-probe output `thermoregulation.py` already consumes. Not a fresh,
  independent measurement of this subject's respiration.
- **FIRST STEP**: no coupling yet to environment/altitude, lung pathology, or — despite
  `thermoregulation.py` already existing — the respiratory contribution to evaporative/convective
  heat loss (real exercise hyperpnea measurably increases respiratory heat/water loss; a genuine,
  not-yet-closed link, since `thermoregulation.py`'s sweat-rate calculation implicitly assumes *all*
  evaporative loss is cutaneous). VCO2/RER is used only as a conceptual cross-check constant here
  (Step 2 of the script), never separately computed as its own twin output.
- **Single subject (subject2), single trial (walking1), one speed (1.065 m/s)** — inherits, does not
  fix or expand, `metabolic_cost.py`'s own scope limits. This speed is itself *slower* than the
  "preferred/self-selected" walking speeds reported in the literature used above (Browning & Kram
  2005/2006: ~1.40–1.47 m/s) — so a match to "preferred-speed" VO2 ranges is a genuine, non-forced
  finding, not a foregone one given how the anchor was chosen.
- **`WebSearch` unavailable this session** (shared budget exhausted) — every citation above was
  verified via direct `WebFetch` to NCBI eutils / PubMed / the Compendium's own site / Wikipedia,
  not a broad web search; a few candidate primary-source numeric tables (e.g. Waters & Mulroy 1999's
  own reported VO2 value, the ATS/ACCP statement's own reference-range table) remained paywalled/
  inaccessible via this route and are flagged as bibliographically-verified-only, per citation table.

## Repro

```
cd ~/projects/bodytwin
.venv-msk/bin/python scripts/msk/respiratory.py
```
Requires `data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json` to
already exist (run `scripts/msk/metabolic_cost.py` first if not). Writes
`data/msk_smoketest/subject2_walking1/respiratory/respiratory_results.json`. Pure Python/numpy, no
OpenSim call, runs in under a second. No git operations; reads the metabolic-cost JSON read-only;
writes only under `data/msk_smoketest/subject2_walking1/respiratory/`.
