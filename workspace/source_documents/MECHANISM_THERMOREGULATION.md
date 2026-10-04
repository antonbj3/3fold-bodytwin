# MECHANISM THERMOREGULATION — metabolic-heat / 0-D heat-balance capability (2026-07-21)

Couples the musculoskeletal twin's metabolic-cost capability (`scripts/msk/metabolic_cost.py`,
`docs/MECHANISM_METABOLIC_COST.md`) to whole-body **heat balance** — the mission's "battery/thermal
verticals = metabolism/calorimetry cert-twins" direction (COORDINATOR.md Sec.5), taken one physical step
further: from metabolic *rate* to metabolic *heat* to a predicted temperature-rise rate and required
sweat rate. Script: `scripts/msk/thermoregulation.py`. Evidence:
`data/msk_smoketest/subject2_walking1/thermoregulation/thermoregulation_results.json`.

**NO re-solve, no OpenSim.** This layer opens `metabolic_cost_results.json` (plain JSON, read-only)
and does pure arithmetic on numbers already computed there. It never touches the `.osim` model or
any `.sto` file.

## Method, in one paragraph

Muscle contraction converts only ~20-25% of metabolic energy into actual mechanical work; the rest
becomes heat (this task's own framing, anchored to the classical measured-efficiency literature,
see Citations). Heat production is split physiologically into two terms: `H_prod = M_rest` (basal/
organ metabolism — ~100% heat, no mechanical work) `+ (1-eta)*(M_gross - M_rest)` (the *exercise
increment* only gets the muscular-efficiency credit, since that's where the actual mechanical work
happens). Both inputs (`M_gross`, `M_rest`) are read directly from the twin's own already-computed
`metabolic_cost_results.json` — no new number is invented for the metabolic rate itself. Two 0-D
(lumped, no spatial resolution) compartments are then evaluated as **zero-dissipation adiabatic
upper bounds** on rate of temperature rise (`dT/dt = H / (mass * specific_heat) * 60`, °C/min): the
**whole body** (mass = this subject's own 78.2 kg, heat includes basal) and the **local working
muscle only** (mass = 17.2 kg, metabolic_cost.py's own literature-anchored lower-limb+hip
contractile-mass estimate, re-derived here from that JSON's own stored correction-scale factor
rather than re-typed as a fresh number; heat = the exercise increment only, since basal heat is
generated throughout the whole body, not concentrated in leg muscle). A **required steady-state
sweat rate** is then computed from the exercise-increment heat alone (dry heat loss assumed flat at
resting level — a disclosed, conservative simplification that over-, not under-, estimates the
requirement) via the latent heat of vaporization of sweat. Finally, an **independent, empirical**
cross-check — Saltin & Hermansen's (1966) directly-measured equilibrium core-temperature-vs-
metabolic-rate relationship, `tcor = 36.6 + M/500` (M in Watts, ~10 min time constant), as verbatim-
quoted and live-fetched from Malchaire (2006) — gives a genuinely different-mechanism (statistical
fit, not this script's own energy-balance arithmetic) prediction of equilibrium core temperature and
a transient time-course, plus an out-of-domain flag against that equation's own stated validity
range (100-450 W).

## Citations — verified LIVE this session, not recalled

| source | verified as | what it anchors |
|---|---|---|
| Cramer MN, Jay O (2016) "Biophysical aspects of human thermoregulation during heat stress," *Auton Neurosci* 196:3-13. **PMID 26971392**, DOI 10.1016/j.autneu.2016.03.001 | Abstract fetched & read live (pubmed.ncbi.nlm.nih.gov) | The CONCEPT: core temp = dynamic balance of production vs. dry+evaporative dissipation. No PMC full text (paywalled, confirmed via NCBI ID-converter) — no numeric constants taken from it. |
| Malchaire JBM (2006) "Occupational Heat Stress Assessment by the Predicted Heat Strain Model," *Industrial Health* 44(3):380-7. **PMID 16922181**, DOI 10.2486/indhealth.44.380 | **FULL TEXT fetched live** (J-STAGE open access) and read directly, not just abstract | The equation `tcor=36.6+M/500` + ~10 min time constant (quoting Saltin & Hermansen 1966 below); max non-acclimatized sweat rate 650-1,000 g/h (citing Araki et al. 1979, Gosselin 1947); their own **measured** validation-dataset sweat rates (lab 424±172 g/h n=672, field 317±187 g/h n=237); the PHS model's own stated metabolic-rate validity domain (100-450 W); the WHO (1969) 38°C occupational core-temp ceiling. |
| Saltin B, Hermansen L (1966) "Esophageal, rectal, and muscle temperature during exercise," *J Appl Physiol* 21(6):1757-62. **PMID 5929300** | Verified live via NCBI esummary (title/journal/year/pages match Malchaire's own citation exactly) | ORIGINAL source of the `tcor=36.6+M/500` equation. Full text NOT independently fetched this session — used via Malchaire (2006)'s verified direct quotation (same secondary-citation-chain discipline `muscle_fatigue.py` already applies to Xia & Frey-Law 2008 via Frey-Law et al. 2012). |
| Whipp BJ, Wasserman K (1969) "Efficiency of muscular work," *J Appl Physiol* 26(5):644-8. **PMID 5781619** | Verified live (bibliographic: title/journal/year/pages match) | The classical reference establishing muscular efficiency as a directly MEASURED quantity. No abstract/full text accessible live (pre-1969-abstracting era, no PMC copy) — the ~20-25% figure used is this task's own stated textbook-consensus value, **not independently re-extracted from this paper's primary text**, flagged exactly as such (same discipline metabolic_cost.py applied to Janssen et al. 2000). |
| Saltin B, Gagge AP, Stolwijk JA (1968) "Muscle temperature during submaximal exercise in man," *J Appl Physiol* 25(6):679-88. **PMID 5727193** | Verified live (bibliographic) | Cited for scope: real muscle temperature during dynamic exercise reaches a bounded, modest steady state (unlike this script's adiabatic bound) because of perfusion. No specific number asserted from it (full text not fetched). |
| Gonzalez-Alonso J (2012) "Human thermoregulation and the cardiovascular system," *Exp Physiol* 97(3):340-6. **PMID 22227198** | Verified live (bibliographic) | Scope pointer to the not-yet-built muscle-perfusion coupling. |
| Sawka MN et al. (2007) ACSM position stand "Exercise and fluid replacement," *Med Sci Sports Exerc* 39(2):377-90. **PMID 17277604** | Verified live (bibliographic) | Context citation for the sweat-rate/fluid-balance consensus domain. |
| Nadel ER (1985) "Recent advances in temperature regulation during exercise in humans," *Fed Proc* 44(7). **PMID 3884384** | Verified live (bibliographic; corrects a common "Nadel ES" mis-recall — author is Nadel ER) | This task brief's own suggested classic-reference family. |
| Specific heat of body tissue (3.49 kJ/(kg·K)), latent heat of vaporization of sweat (2,426 J/g) | **TEXTBOOK-GRADE, flagged, NOT independently re-verified live** | Standard constants from the Gagge two-node / ISO 7933 lineage Malchaire (2006) itself belongs to. Live fetch attempts (Wikipedia Perspiration/Thermoregulation pages, same discipline as metabolic_cost.py's MET conversion) did not surface either exact number — disclosed as a genuine open gap rather than silently asserted from recall. Wikipedia's Perspiration page DID yield an independent sweat-rate cross-check used below (fetched live): "maximum sweat rates... up to 2-4 L/h," "average intensity exercise... up to 2 L/h." |

Recall discipline: this project's own prior finding is a measured ~62% citation-drift rate from
memory (see `MECHANISM_METABOLIC_COST.md`). Every PMID/DOI/equation above was fetched live this
session (NCBI eutils, NCBI PMC ID-converter, a full-text PDF fetch+read for Malchaire 2006, and
Wikipedia for the two flagged textbook constants) — not recalled.

## Headline results (subject2/walking1, mass 78.2 kg; three metabolic-rate inputs)

Three inputs are carried forward **end-to-end**, not just one, because `metabolic_cost.py`'s own
prior finding was an honest SURPRISE (primary/uncorrected config landed ABOVE, not below, its
external anchor): **Umberger-primary** (733.6 W), **Bhargava-primary** (609.9 W), and
**combined-corrected** (400.4 W, = that script's own disclosed tendon-compliance + muscle-mass
sensitivity correction). Basal rate (93.8 W) is cross-model self-consistent (derived from the JSON
itself, not re-typed) to machine precision.

| quantity (η=0.225 mid) | Umberger-primary | Bhargava-primary | combined-corrected |
|---|---:|---:|---:|
| Heat production `H_prod` | 589.6 W (80.4% of M_gross) | 493.8 W (81.0%) | 331.4 W (82.8%) |
| Whole-body dT/dt (zero-dissipation bound) | 7.78 °C/hour | 6.51 °C/hour | 4.37 °C/hour |
| Local working-muscle dT/dt (adiabatic, no perfusion) | 29.7 °C/hour | 24.0 °C/hour | 14.2 °C/hour |
| Required steady-state sweat rate | 735.7 g/h (0.74 L/h) | 593.5 g/h (0.59 L/h) | 352.5 g/h (0.35 L/h) |
| Saltin-Hermansen equilibrium core temp | 38.07 °C (+1.28 °C) | 37.82 °C (+1.03 °C) | 37.40 °C (+0.61 °C) |
| Within Malchaire's own M-validity domain (100-450 W)? | **NO** (733.6 W) | **NO** (609.9 W) | **YES** (400.4 W) |

All numbers, gates, and the full sensitivity sweep are on disk in
`thermoregulation_results.json`; nothing above is hand-computed prose — the script prints and JSON-
serializes every one of these from its own arithmetic.

## The decisive external anchors (not a tautology gate)

1. **Malchaire's own MEASURED validation dataset**, an anchor this script's arithmetic was never fit
   to: at M=450 W (the top of the PHS model's own stated validity range), this script's required-
   sweat-rate calculation gives **409.6 g/h**, against Malchaire's own independently measured
   laboratory dataset (n=672): **424 ± 172 g/h** — inside 1 SD. A second, textbook-grade cross-check
   (Wikipedia, flagged) — "average intensity exercise up to 2 L/h, adult max 2-4 L/h" — the twin's
   full 0.35-0.74 L/h span across all three metabolic-rate inputs sits comfortably inside.
2. **Saltin-Hermansen's empirical equilibrium-core-temperature equation**, an entirely different
   physiological mechanism (a directly-measured statistical fit vs. this script's own energy-balance
   arithmetic): at pure rest (M=93.8 W) it predicts tcor=36.79 °C — matching normal resting core
   temperature as a sanity floor. At the twin's three exercise inputs it predicts +0.61 to +1.28 °C
   rises, landing in the range well-documented for moderate-to-vigorous sustained exercise.

**SURPRISE, symmetric with metabolic_cost.py's own prior surprise, corroborated via a completely
independent mechanism.** `metabolic_cost.py` already disclosed that its PRIMARY configuration
over-estimates metabolic rate relative to the (Koelewijn, cost-of-transport) external anchor, with
a "combined correction" sensitivity variant landing much closer to it. This script's OWN,
independently-sourced anchor (Saltin-Hermansen/Malchaire, a core-temperature equation, not a
cost-of-transport comparison) reaches the **same qualitative conclusion by a different route**: the
Umberger- and Bhargava-primary inputs (733.6 W, 609.9 W) both fall **outside** the PHS model's own
stated M-validity domain (100-450 W) — consistent with those inputs representing more like a jog
than a walk at 1.065 m/s — while the corrected input (400.4 W) stays inside it. Two decorrelated
external anchors, from two different physiological literatures, now agree that the corrected number
is the more externally-consistent one to carry forward.

## Void-floor / non-degeneracy check (forced, not eyeballed)

Concern (the adversary this kind of thermal-conversion layer is tempted to skip): is this just
reproducing generic textbook thermophysiology ranges regardless of the twin's own specific number —
i.e. would ANY plausible input produce a "plausible-looking" answer, making the whole cross-check
non-discriminating? Forced via a 25-point sweep of the input metabolic rate across
[150 W, 1200 W] (light activity to near-maximal sustained human exercise, machine-checked in
`thermoregulation.py`'s Step 7): both the required-sweat-rate and dT/dt outputs are **strictly
increasing** across the entire sweep (required sweat rate spans 64.6-1272.1 g/h; dT/dt spans
0.030-0.209 °C/min), and the numerically-estimated derivative matches the closed-form analytical
derivative (`d(sweat)/dM = (1-η)/L_vap*3600 = 1.1500 g/h per W`) to machine precision. **A real,
non-pinned function of the twin's own metabolic-rate number, not a constant dressed up as a
response.**

One genuine bug was caught and fixed by this same check, not swept under the rug: the first attempt
swept from M=90 W, which is *below* this subject's own ~93.8 W resting rate — the required-sweat-
rate function has a physically-correct `max(excess_heat, 0)` floor there (you cannot have negative
required sweat rate when there is a heat *deficit*, not excess), creating a real, expected kink at
M=M_rest that a blanket derivative-match check straddled and failed. Diagnosed via direct
`numpy.gradient` inspection (confirmed: mismatch was confined to exactly the kink point and its
central-difference neighbor), not patched over — the fix was to sweep the physiologically meaningful
exercise-intensity domain (M > M_rest), which is both the more correct test design and the more
correct physiological question for this layer.

## Scope, honest gaps (first-step, not final)

- **0-D lumped, not spatial.** No temperature field, no FEM/PDE conduction, no skin-layer resolution.
  "Muscle temperature" is one lumped number for 17.2 kg of working muscle, not a spatial map.
- **No vascular/perfusion heat transport modeled — the single biggest structural gap, and the one
  this script's own numbers make legible rather than hiding in prose.** The local-muscle adiabatic
  bound (14-30 °C/hour depending on input) is *physiologically unsustainable if taken literally* —
  real muscle temperature during dynamic exercise reaches a bounded, modest new steady state within
  minutes (Saltin, Gagge, Stolwijk 1968, PMID 5727193) precisely because continuously increasing
  local blood flow carries heat away, which this model does not represent at all. This is the
  natural next coupling: the muscle-perfusion agent's (not yet built) convective output would replace
  this script's "zero heat removal" assumption with a real term, turning this diagnostic upper bound
  into an actual transient prediction.
- **Generic constants** (efficiency, specific heat, latent heat) — not this subject's own measured
  body composition, sex, age, VO2max, or heat-acclimatization state (Malchaire's own review notes
  actual acclimatized sweat rate can run ~2× non-acclimatized in a given environment, though max
  *capacity* only rises ~25%; none of that individual variation is modeled here).
- **No environment.** No ambient temperature, humidity, wind, or clothing insulation — "required
  sweat rate" is evaporative *demand* only; whether that demand is physically achievable depends on
  ambient wet-bulb conditions this script does not represent.
- **Required-sweat-rate is itself a disclosed upper bound**, not a best estimate: it assumes dry
  heat loss (radiation+convection+conduction) stays flat at the resting value, when in reality dry
  heat loss rises somewhat with exercise-elevated skin temperature — so true required sweat rate is
  probably somewhat *below* the numbers reported here.
- **The (1-η) efficiency convention itself is a simplification specific to how this task frames the
  problem, disclosed as such in the code**: in level, constant-speed walking there is no net external
  mechanical work (nothing permanently leaves the body-environment system), so a fully rigorous
  accounting would eventually count closer to 100% of the exercise increment as heat, not 75-80%. The
  script computes this η=0 alternative explicitly (+21-24% higher heat production than the η=0.225
  primary) as a disclosed sensitivity bound, not a hidden caveat.
- **Single subject (subject2), single trial (walking1), one speed (1.065 m/s)** — inherits, does not
  fix or expand, `metabolic_cost.py`'s own scope limits.
- **Directly inherits and propagates `metabolic_cost.py`'s own disclosed input uncertainty** (Fmax-
  derived muscle-mass over-estimate, tendon-compliance choice, SO's co-contraction under-estimate)
  rather than silently picking one number — shown explicitly via the three-input comparison above.
- **Transient time-course borrows a time constant from a different-population empirical fit**
  (Saltin & Hermansen's ~10 min), applied to this twin's own metabolic rate as an explicit projection,
  not a measurement over an actual sustained trial (the underlying measured stride is 1.57 s; every
  *rate* reported is well-defined at that instant and needs no extrapolation, but sustaining that rate
  for 10-40 minutes to approach equilibrium is a disclosed assumption).

## Repro

```
cd ~/projects/bodytwin
.venv-msk/bin/python scripts/msk/thermoregulation.py
```
Requires `data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json` to
already exist (run `scripts/msk/metabolic_cost.py` first if not). Writes
`data/msk_smoketest/subject2_walking1/thermoregulation/thermoregulation_results.json`. Pure
Python/numpy, no OpenSim call, runs in under a second. No git operations; reads the metabolic-cost
JSON read-only; writes only under `data/msk_smoketest/subject2_walking1/thermoregulation/`.
