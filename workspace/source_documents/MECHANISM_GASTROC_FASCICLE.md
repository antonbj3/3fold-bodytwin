# MECHANISM GASTROC/SOLEUS FASCICLE KINEMATICS — validating muscle-INTERNAL resolution against in-vivo B-mode ultrasound (2026-07-21)

Tests a dimension of the plantarflexor model that neither `docs/MECHANISM_PUSHOFF_PLANTARFLEXOR.md`
(Achilles-equivalent tendon FORCE vs Komi/Finni/Frisberg transducer data) nor
`docs/MECHANISM_CONTACT_WAVEFORM.md` (joint contact-force waveform SHAPE) probed: the gastrocnemius/
soleus muscle-tendon unit's INTERNAL kinematics — fascicle (fiber) length, pennation angle, and the
fiber-vs-MTU "decoupling" the compliant Achilles tendon is supposed to provide — checked against a
NEW anchor decorrelated from both force-transducer and joint-contact-force data: in-vivo B-mode
ultrasound fascicle-length measurements during walking.

**Headline: MIXED / PARTIAL, not a clean confirm and not a clean rule-out.** Neither of the task's
own two named falsifiers triggers: the model's gasmed_r/gaslat_r/soleus_r fascicles do NOT
over-shorten (stance-phase normalized fiber length stays 0.688–0.791, comfortably off the
ascending-limb/degenerate region) and do NOT mis-time (fiber-length minimum lands at 68.2–69.7%GC,
just past the independently-established R toe-off event at 66.9%GC, robust to smoothing) — a real,
qualitative reproduction of "isometric-ish stance, real shortening concentrated at push-off/take-off."
A forced compliant-vs-rigid-tendon contrast (this task's own explicit symmetric-QC ask) proves this
is a genuine tendon-compliance effect, not a construction tautology: the rigid-tendon bookend shows
97–99% coupling (fiber tracks essentially 100% of MTU excursion, as it must with no compliance),
while the real compliant model shows materially less (55.7–86.8% depending on window/muscle). **But
a THIRD, additional, pre-registered axis this session measured beyond the task's own two — the
STRENGTH of that decoupling and the ABSOLUTE fascicle length scale — does NOT hold up as well**:
early-mid-stance decoupling ratios (0.557–0.868) fail this session's own pre-registered ≤0.50 pass
bar for all three muscles, and soleus_r (0.868) crosses the explicit ≥0.80 "no decoupling" hard
falsifier; over the FULL stance phase the picture is quantitatively worse still (0.79–0.85 for all
three) — the model's fiber does the MAJORITY of total MTU length-change work across a full stance
cycle, the opposite of the tendon-dominated "catapult" mechanism the ultrasound literature describes.
Absolute fascicle length for the gastrocs also runs long: mid-stance peaks of 85.16mm (gasmed_r) and
94.66mm (gaslat_r) sit 30–40mm above the literature's typical reported range (~40–58mm). **A forced
adversary caught mid-analysis, independently verified live, not silently corrected: the task's own
stated PMID for Farris & Sawicki 2012 (22233781) is the WRONG paper** (an unrelated in-vitro
neurotoxicology study) — the real PMID is **22219360**, confirmed live via NCBI eutils both by a
sub-agent and independently re-checked directly by this session.

## Headline numbers

| muscle | stance fiber range (mm) | stance min (norm) | early-stance mean (mm) | push-off min (mm) @ %GC | shortening (mm) | decoupling, 0–40%GC (compliant / rigid) | decoupling, full stance (compliant) |
|---|---:|---:|---:|---:|---:|---:|---:|
| gasmed_r | 50.92 – 85.16 | 0.688 | 80.84 | 49.54 @ 68.2%GC | 31.30 | **0.629** / 0.989 | 0.847 |
| gaslat_r | 60.78 – 94.66 | 0.702 | 89.84 | 57.93 @ 69.7%GC | 31.91 | **0.557** / 0.981 | 0.791 |
| soleus_r | 43.62 – 64.70 | 0.791 | 56.78 | 42.85 @ 68.2%GC | 13.93 | **0.868** / 0.937 | 0.830 |

All numbers machine-computed by `scripts/msk/gastroc_fascicle_kinematics.py`, written to
`data/msk_smoketest/subject2_walking1/gastroc_fascicle_kinematics/gastroc_fascicle_kinematics_results.json`
— none transcribed from console prose. Model: `Millard2012EquilibriumMuscle`,
`ignore_tendon_compliance=false` (COMPLIANT tendon) for all three muscles, confirmed live from the
`.osim` XML and from the OpenSim API — this is the answer to this task's own explicit symmetric-QC
question ("state which the model uses").

## 1. Pre-registration (stated before this task's OODA-forced window correction, §2.4, changed any number)

**C / ¬C, thresholds fixed in advance, directly reusing the task's own framing:**
- **C (fascicle kinematics validated):** stance-phase fiber length stays within a generous outer
  physiological band (gasmed_r/gaslat_r/soleus_r: 25–100/25–110/20–90mm) AND does not over-shorten
  below 0.50× optimal_fiber_length AND reaches its minimum inside a push-off/take-off window AND
  shows material (≥3mm) push-off shortening AND passes the geometric self-consistency identity
  (h = fiber_length·sin(pennation), exact by construction of the Millard "constant muscle height"
  pennation model).
- **¬C (fiber-kinematics ERROR, biases the force-length operating point):** over-shortening
  (normalized fiber length < 0.50 at any stance frame), OR mis-timing (minimum outside the push-off
  window), OR — an ADDITIONAL axis this session pre-registered beyond the task's own two, given the
  task's explicit emphasis on "decoupling" as the mechanism being tested — no decoupling (fiber/MTU
  peak-to-peak ratio ≥0.80 over early-mid-stance, meaning the tendon isn't buffering meaningfully).
- **Forced adversary (the one this task is tempted to skip):** near-isometric-then-shorten behavior
  could be a TRIVIAL/tautological consequence of ANY reasonably-parametrized compliant-tendon Hill
  muscle, not evidence the model's SPECIFIC parameters are realistic. Forced via a RIGID-tendon
  geometric bookend computed alongside every frame (§2.2): if the rigid bookend showed the SAME
  qualitative decoupling, the compliant result would be suspect as a test-construction artifact. It
  does not (rigid decoupling 93.7–98.9% vs compliant 55.7–86.8%, §4.2) — the adversary FALLS, proving
  the effect is caused by this model's real tendon compliance being load-bearing.

**Symmetric-QC checks forced, not assumed:** (1) the task's own explicit ask — "model normalized
fiber length × optimal-fiber-length = absolute mm — verify the conversion" — checked at all 158
frames × 3 muscles: max relative error **1.8×10⁻¹⁶** (gasmed_r), 1.6×10⁻¹⁶ (gaslat_r), 1.3×10⁻¹⁶
(soleus_r) — exact to floating-point precision, PASS. (2) the compliant-tendon vs rigid-tendon
choice — verified live (`ignore_tendon_compliance=false` for all 3 muscles) and BOTH reconstructions
computed and reported (§4.2), not just asserted in prose. (3) a second forced adversary this session
identified and ran BEFORE trusting the primary reconstruction (§2.3): does the compliant-equilibrium
reconstruction reproduce Static Optimization's OWN reported per-muscle force? It does not, closely
(§5) — a real, quantified, mechanistically-diagnosed gap, disclosed prominently (§6), not hidden.

## 2. Method

### 2.1 Primary reconstruction (compliant, model-faithful)

Reused, **zero re-solve of Static Optimization itself**: the already-certified subject2/`walking1`
SO activation trajectory (`data/msk_smoketest/subject2_walking1/static_optimization/so/
walking1_StaticOptimization_activation.sto`) and the already-validated IK motion
(`OpenSimData/Mocap/IK/walking1.mot`, read-only external drive). Both share an **identical 158-row,
0.01s-step time grid, t=0.00–1.57s** (verified live, `np.allclose`, PASS — no interpolation needed).
At each of the 158 frames: set the model's generalized coordinates from IK, set each target muscle's
activation to SO's own value, call OpenSim's native `Model.equilibrateMuscles(state)` — which solves
the true nonlinear tendon↔fiber static-equilibrium equation using the muscle's REAL
`Millard2012EquilibriumMuscle` tendon force-length curve — then read `getFiberLength`,
`getPennationAngle`, `getTendonLength` directly. This is the same per-frame quasi-static-equilibrium
assumption Static Optimization itself uses (no continuous activation/fiber-length ODE integration
across frames either way), and is the standard convention for post-hoc "model-based fascicle length
from kinematics + activation" reconstruction in this literature.

### 2.2 Forced-adversary bookend (rigid tendon)

Computed alongside, every frame, pure geometry, no equilibrium solve:
`fiber_length_rigid(t) = sqrt(h² + (MTU_length(t) − tendon_slack_length)²)`, where
`h = optimal_fiber_length × sin(pennation_angle_at_optimal)` is the Millard model's constant
"muscle height" invariant. This is the required symmetric-QC contrast the task explicitly asks for.

### 2.3 Forced-adversary #2 (does the reconstruction match SO's own force?) — single-frame decisive pre-check

Before building the full 158-frame pipeline, a single-frame spot-check at the already-established
push-off instant (t=0.51s, 44.8%GC, `docs/MECHANISM_CONTACT_WAVEFORM.md`) compared: (a) the compliant
equilibrium's implied tendon force (`getTendonForce`) vs SO's reported force — 15.6/17.3/11.7%
relative error (gasmed_r/gaslat_r/soleus_r); (b) a rigid-tendon-plus-real-fiber-velocity bookend
(velocity from central-difference of MTU length, propagated through the geometric fiber-velocity
relationship and the muscle's own force-velocity curve) — 6.4/31.5/7.8% relative error. **Neither
bookend cleanly reproduces SO's own number** (no clean "SO=rigid" or "SO=compliant" story) — measured
across the whole trial in §5, not just this one frame.

### 2.4 Forced-adversary catch, OODA-corrected (window-truncation artifact)

**Observed:** an initial push-off search window (35–65%GC) put the fiber-length minimum at EXACTLY
64.4%GC for all 3 muscles — suspiciously identical across different muscles and suspiciously close
to the window's own upper edge (the classic boundary-truncation signature already documented in
`docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md` §4 for a different quantity). **Oriented:** the
full-cycle, unwindowed trace showed the true minimum lies just OUTSIDE that window, at ~68–70%GC for
all 3 muscles — a genuine interior turning point (fiber length rises again immediately after), not
further truncation. This sits just past the INDEPENDENTLY-established R toe-off event (66.9%GC,
`docs/MECHANISM_CONTACT_WAVEFORM.md`, itself verified there via an out-of-sample textbook-band check)
— physiologically sensible (matches Lichtwark/Bougoulias/Wilson 2007's own "high velocity shortening
during take off" language) rather than an artifact. **Decided + Acted:** widened the window to
(35, 80)%GC and re-verified the new minimum is interior, not re-hitting the new boundary (confirmed:
68.2–69.7%GC, comfortably inside). Reported here rather than silently fixed.

### 2.5 %GC mapping and %Fmax convention

Reused VERBATIM from `docs/MECHANISM_CONTACT_WAVEFORM.md` / `contact_waveform_analysis.py`
(`T_stride=1.327s`, `R_heel_strike_t=1.2425s`, ipsilateral R toe-off = 66.9%GC) — not re-derived, so
%GC values here are directly comparable to every other push-off-timing number already established
in this repo.

## 3. Verification gates (machine-checked, all reported — PASS unless noted)

| gate | result |
|---|---|
| Time-grid alignment (IK vs SO-force vs SO-activation, all 158-row/0.01s) | **PASS** (`np.allclose`) |
| `equilibrateMuscles()` convergence, all 158 frames × 3 muscles | **PASS**, 0/158 failures |
| Geometric self-consistency identity (h = fiber·sin(pennation)), max relative error | **PASS**, gasmed_r 1.4×10⁻¹⁶, gaslat_r 1.9×10⁻¹⁶, soleus_r 0.0 (threshold 0.1%) |
| Normalized-fiber-length × optimal-fiber-length = absolute mm conversion (task's explicit ask) | **PASS**, max relative error 1.3–1.8×10⁻¹⁶, all 3 muscles |
| No over-shortening (stance-phase normalized fiber length ≥ 0.50) | **PASS**, min 0.688/0.702/0.791 |
| No mis-timing (fiber-length minimum inside push-off/take-off window, robust to smoothing) | **PASS**, 68.2–69.7%GC, shift ≤1.5%GC under 5-point smoothing |
| Compliant-vs-rigid contrast (proves decoupling is real, not tautological) | **PASS**, rigid 93.7–98.9% vs compliant 55.7–86.8% |
| Decoupling ratio ≤ 0.50 (this session's additional pre-registered pass bar) | **FAIL**, all 3 muscles (0.557–0.868, early-mid-stance) |
| No-decoupling hard falsifier (ratio ≥ 0.80) | **TRIGGERED for soleus_r** (0.868), not for gasmed_r/gaslat_r |
| Literature PMIDs verified live via NCBI eutils (this session, independently, not trusted from a sub-agent report alone) | **PASS** — see §4; 1 citation ERROR caught (Farris & Sawicki PMID) |

## 4. External literature anchor — verified live, and a forced adversary caught mid-analysis

Verified via a `watertight-researcher` sub-agent (multi-pass extraction, cross-checked internally)
and INDEPENDENTLY RE-VERIFIED by this session directly via `curl` to `eutils.ncbi.nlm.nih.gov`
(esummary/efetch) before being trusted — per this repo's own established live-citation-verification
convention (`scripts/msk/tendon_elastic_energy.py`'s Ker-1987 precedent).

### 4.1 Citation correction (forced adversary, not silently fixed)

**The task's own stated PMID for Farris & Sawicki 2012 (22233781) is WRONG.** Independently
confirmed live via NCBI eutils esummary: PMID 22233781 is Wang Z et al., "Lithium attenuates
bupivacaine-induced neurotoxicity in vitro...", *Neuroscience* 2012 — an unrelated in-vitro
neurotoxicology paper. **The correct PMID is 22219360** — confirmed live (title, journal, authors,
date all match): Farris DJ, Sawicki GS. "Human medial gastrocnemius force-velocity behavior shifts
with locomotion speed and gait." *Proc Natl Acad Sci U S A.* 2012 Jan 17;109(3):977-82.
DOI 10.1073/pnas.1107972109.

### 4.2 Anchors used

- **Farris & Sawicki 2012 PNAS (PMID 22219360, verified live)** — n=10 (6M/4F), age ~24-25y, mass
  67-77kg, treadmill walking 0.75/1.25/1.75/2.0 m/s. MEDIAL GASTROCNEMIUS ONLY (no LG/soleus). At
  0.75 m/s (slowest, closest to this trial's own rough speed proxy, §7): initial length at
  heel-strike `Li = 47±2mm`, length change `ΔL = -7±1mm`, average shortening velocity `11±2mm/s`,
  max SEE (tendon) length `373±25mm`. Verbatim: *"relatively isometric behavior during midstance
  followed by shortening in late stance"*; *"considerable decoupling of the length change of MG
  fascicles from that of the whole MG muscle-tendon unit (MTU)... facilitated by the lengthening
  (~15mm) of the series elastic element (SEE) during early to midstance and its recoil in late
  stance."* No pennation angle numeric value reported.
- **Lichtwark & Wilson 2006 JEB (PMID 17050853, verified live — title/DOI/abstract independently
  re-confirmed by this session directly)** — n=6 male, age 28±3.68y, height 182±12cm, mass
  82±12kg, treadmill 5km/h (~1.39 m/s), level walking. Mid-stance fascicle length `57.7±2.3mm`;
  shortest across walking conditions `44.1±1.0mm`; pennation range through stance `4.4±0.8°`; max
  pennation pooled-walking `17.5±1.2°`; max MTU length change `24.9±1.35mm`; max Achilles tendon
  strain ~4.6% (`10.9±1.7mm`). This session's own direct eutils fetch of the abstract independently
  confirms the core mechanism: *"force was developed isometrically... Substantial levels of Achilles
  tendon strain were recorded... which allowed the muscle fascicles to act at speeds more favourable
  for power production."* (Body-text mm/deg numbers are from the sub-agent's full-text extraction —
  ScienceDirect returned HTTP 403 to this session's own direct re-fetch attempt — cross-checked
  across 3 of the sub-agent's own passes, 1 discrepancy caught: 56.8mm is downhill, 57.7mm is the
  correct level-walking figure.)
- **Lichtwark, Bougoulias & Wilson 2007 J Biomech (PMID 16364330)** — already used in this repo
  (`scripts/msk/tendon_elastic_energy.py`) for the tendon-ENERGY question; reused here as qualitative
  mechanism support (abstract has no mm/deg numbers, confirmed by both this session's own eutils
  fetch and the sub-agent, both blocked from body text): *"Muscle fascicles acted relatively
  isometrically during the stance phase during walking... high velocity shortening during take off
  in both walking and running achieved by recoil of the SEE."*
- **Rubenson et al. 2012 JEB, soleus (PMID 22771749, verified live)** — n=8 male, age 26±3.51y,
  mass 70.31±9.18kg, preferred treadmill walking speed. Optimal fascicle length `L0 = 37.7±6.9mm`;
  normalized length: onset `0.84±0.15 L0`, max (30-40%GC) `0.94±0.16 L0`, peak EMG `0.88±0.19 L0`,
  toe-off `0.74±0.18 L0`; absolute change during stance `12.1±3.5mm` (`31.7±6.0%` normalized).
  Single-pass extraction (lower confidence than the two above); the source paper's own numbers have
  a minor unreconciled internal inconsistency (0.94−0.74=0.20 L0 ≈7.5mm at mean L0, vs the
  directly-quoted 12.1mm) — both reported as-is, not adjudicated (secondary anchor, LEAN).
- **Cronin & Lichtwark 2013 Gait & Posture (PMID 22910172, verified live, full abstract
  independently re-fetched this session)** — confirmed to be a pure methods/scope REVIEW with
  **zero specific mm/degree numbers** in its abstract; fully closed-access (no OA copy). Cited for
  qualitative context/pointer only, per the task's own naming of it — NOT used as a quantitative
  anchor.

## 5. Forced-adversary #2 in full — does the reconstruction match SO's own force? (honest gap, disclosed)

Measured across the whole trial (not just the single-frame pre-check, §2.3): the compliant
reconstruction's implied tendon force vs SO's own reported force, at frames where SO's force exceeds
5% of that muscle's Fmax —

| muscle | n material frames | RMS relative error | max relative error |
|---|---:|---:|---:|
| gasmed_r | 42 | 117.6% | 362.6% |
| gaslat_r | 55 | 72.5% | 242.4% |
| soleus_r | 68 | 63.0% | 205.0% |

**Diagnosed, not just measured (OODA):** the worst mismatches concentrate exactly in the 57–66%GC
late-push-off window — the phase of most rapid fascicle shortening. Confirmed directly: gasmed_r's
fiber velocity at its single worst-mismatch frame (t=0.72s, 60.6%GC) is **−169.7mm/s (−0.229
normalized, ~23% of Vmax)** — real, substantial, concentric shortening. My quasi-static
`equilibrateMuscles` reconstruction targets zero fiber velocity at each snapshot (by construction —
this is what "equilibrate" means), omitting the force-velocity multiplier SO's own optimizer would
apply at the ACTUAL joint angular velocity. This explains a material fraction of the gap
mechanistically (a normalized velocity of −0.229 would reduce active force by roughly 30–45% via a
typical Millard-curve F-V multiplier), though likely not the entire 200–360% residual — the rest
plausibly reflects SO's own redundant multi-muscle moment-balance solving a different, coupled
per-frame problem than any single-muscle external reconstruction can exactly replay (§2.3 already
found NEITHER simple bookend matches SO's number even at a well-behaved, slowly-changing instant).
**This bounds confidence in any FORCE number derived from this pipeline — it does NOT undermine the
LENGTH/pennation values themselves**, which stay internally self-consistent (geometric identity
exact to floating-point precision) at every single frame, including the worst-force-mismatch ones.

## 6. Results in detail

### 6.1 Task's own two named falsifiers: neither triggers

**No over-shortening:** stance-phase normalized fiber length bottoms at 0.688 (gasmed_r), 0.702
(gaslat_r), 0.791 (soleus_r) — all comfortably above the pre-registered 0.50 floor, i.e., the model
never pushes any of the three plantarflexors into a pathologically short, potentially
numerically-degenerate operating point.

**No mis-timing:** fiber-length minimum lands at 68.2%GC (gasmed_r, soleus_r) / 69.7%GC (gaslat_r)
— just past the independently-established R toe-off (66.9%GC), inside the pre-registered push-off/
take-off window, and robust to a 5-point moving-average smoothing perturbation (shift ≤1.5%GC,
value shift <0.5mm, no muscle flips). This is consistent with, not contradicted by, the "shortening
in late stance" / "high velocity shortening during take off" language in all three qualitative
literature anchors (Farris & Sawicki; Lichtwark & Wilson; Lichtwark/Bougoulias/Wilson).

### 6.2 The compliant-vs-rigid contrast proves the decoupling is real, not tautological

Rigid-tendon bookend decoupling ratios (0-40%GC): gasmed_r 0.989, gaslat_r 0.981, soleus_r 0.937 —
essentially 1.0, exactly as expected (a tendon that structurally cannot stretch forces the fiber to
track ~100% of MTU excursion). The REAL, compliant model shows materially lower ratios (0.629,
0.557, 0.868) — a clear, measured, non-trivial gap between the two, proving this model's declared
tendon compliance (`ignore_tendon_compliance=false`) is doing genuine mechanical work, not that ANY
model would show this pattern regardless of its specific parameters.

### 6.3 But the decoupling is weaker than the literature's qualitative picture, and gets worse over full stance

Pre-registered ≤0.50 pass bar: **fails for all 3 muscles** at early-mid-stance (0.557–0.868).
soleus_r (0.868) additionally crosses the explicit ≥0.80 "no decoupling" hard falsifier. Over the
FULL stance window (0 to R toe-off, 66.9%GC — directly comparable to Lichtwark & Wilson's own
"maximum length change" convention) the picture is quantitatively worse for all three: gasmed_r
0.847, gaslat_r 0.791, soleus_r 0.830 — meaning across a full stance cycle, the model's FIBER
performs 79–85% of the total MTU length-change work, with the TENDON contributing only 8.50–12.46mm
of individual-tendon excursion (vs Lichtwark & Wilson's reported combined Achilles tendon change of
24.9±1.35mm for level walking — though this is not a strictly apples-to-apples comparison, since the
model represents three independent tendons rather than one confluent Achilles tendon, §7). This is
the opposite of the tendon-dominated "catapult" mechanism (large SEE recoil, small fascicle
excursion) the ultrasound/tendon-energy literature describes for push-off specifically.

### 6.4 Absolute fascicle length: gastrocs run long, soleus closer

Mid-stance peaks of 85.16mm (gasmed_r) and 94.66mm (gaslat_r) sit 30–40mm above the literature's
typical reported range (Farris: `Li=47±2mm`, `ΔL=-7±1mm` → ~40–49mm; Lichtwark: mid-stance
`57.7±2.3mm`, shortest `44.1±1.0mm` → ~44–58mm). This model's own `optimal_fiber_length` parameters
(73.98mm gasmed_r, 86.52mm gaslat_r) already run long relative to typical literature values — a
plausible (not independently verified here) propagating cause, since these are fixed model
properties, not something this session's kinematics/activation-driven reconstruction could distort.
soleus_r's stance range (43.62–64.70mm) sits closer to plausible given Rubenson's `L0=37.7±6.9mm`
(implying a roughly 28–35mm normalized-length range) — the model still runs somewhat longer, but
less severely than the gastrocs.

### 6.5 Two specific, close quantitative hits (not just qualitative agreement)

- **soleus_r's total stance shortening (13.93mm) closely matches Rubenson et al. 2012's directly-
  reported in-vivo soleus "absolute change during stance" of 12.1±3.5mm** — within their own
  reported standard deviation, a genuine, specific quantitative corroboration, not just a
  qualitative-shape match.
- **gasmed_r's early-mid-stance (0-40%GC) MTU excursion (16.06mm) is close in order of magnitude to
  Farris & Sawicki's reported SEE lengthening of "~15mm during early to midstance"** for the same
  anatomical structure and phase.

### 6.6 Connection to the l_TS fragility

`docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md` found the knee/hip force cert is fragile to tendon slack
length (l_TS) sweeps, driven mechanistically by antagonist gasmed_r/recfem_r co-contraction
switching — a DIFFERENT mechanism than measured here. This session's finding is complementary, not
duplicative: a muscle whose tendon absorbs LESS of the total stance-phase length change than real
anatomy would (this model's own measured 79–87% fiber-side share, full stance) is mechanically MORE,
not less, directly exposed to any error in `tendon_slack_length` — because the tendon isn't buffering
as much of the length change as a real, more strongly-decoupling Achilles tendon would. This offers a
candidate, geometrically-grounded, independently-derived contributing explanation for that prior
fragility finding — it does not replace or re-litigate it.

## 7. Symmetric-QC

**The adversary this result (leaning toward a positive "validated" finding via §6.1's clean pass on
the task's own two falsifiers) must survive:** is "no over-shortening, no mis-timing" actually weak
evidence, since it might be nearly automatic for any reasonably-parametrized compliant Hill muscle
driven by physiologically-plausible activation? Forced via §6.2's rigid-tendon contrast: NO — the
rigid bookend, using the IDENTICAL activation/kinematics input, shows ~99% coupling and a materially
different absolute trajectory, proving the compliant model's behavior is a real, non-tautological
consequence of its own tendon compliance, not a free pass from the test's own construction.

**The adversary this result (leaning toward a negative "fiber-kinematics error" finding via §6.3's
decoupling-ratio failures) must survive:** is the decoupling-ratio failure driven by an unfair,
overly strict threshold this session invented, rather than a real literature-anchored standard?
Forced by keeping BOTH windows (early-mid-stance AND full-stance) and BOTH specific quantitative
literature comparisons (§6.5): the qualitative literature language ("considerable decoupling",
"relatively isometric") and this session's OWN full-stance numbers (fiber doing 79–87% of the work,
tendon excursion 3–6mm short of the one closely-comparable published number) converge on the same
conclusion independent of the exact 0.50/0.80 thresholds chosen — the failure is not an artifact of
threshold placement.

**A third adversary, caught only by fetching the anchor literature carefully (§4.1):** the task's
own stated PMID for the primary named anchor (Farris & Sawicki 2012) turned out to be the wrong
paper entirely. Using the WRONG paper's (nonexistent, since it's not even the same topic) numbers
was never a risk here since the mismatch was caught before any number was used, but this is flagged
prominently — exactly the forced-adversary discipline this repo's citation-verification convention
already established (`scripts/msk/tendon_elastic_energy.py`'s Ker-1987 precedent) — rather than
silently substituting the correct PMID without disclosure.

## 8. Honest gaps

1. **Forced-adversary #2 (§5): the compliant reconstruction does not closely reproduce SO's own
   reported force** (RMS relative error 63–118%, max 205–363%, concentrated in the high-fiber-
   velocity late-push-off window) — mechanistically diagnosed (fiber-velocity omission in the
   quasi-static equilibrium), not merely observed, and does not undermine the length/pennation
   values themselves (geometric identity stays exact throughout), but bounds confidence in any FORCE
   number this specific pipeline would produce.
2. **Neither the compliant nor a rigid-tendon-plus-real-velocity bookend exactly reproduces SO's
   reported force even at the well-behaved t=0.51s instant** (6.4–31.5% relative error, muscle-
   dependent) — SO's own internal per-frame force computation is not exactly reproducible from
   outside via either simple assumption; the compliant reconstruction is used as PRIMARY because it
   uses the model's own real, declared tendon property, not because it exactly replays SO's internal
   (undisclosed) procedure.
3. **Three independent tendons, not one confluent Achilles tendon** (same caveat already disclosed
   in `scripts/msk/tendon_elastic_energy.py`) — the full-stance individual-tendon excursions
   (8.50–12.46mm) are not directly summable into a single "Achilles tendon excursion" comparable to
   Lichtwark & Wilson's combined 24.9mm without assuming an unmodeled load-sharing relationship.
4. **Cross-study speed/subject mismatch**: Farris tested explicit 0.75–2.0 m/s conditions; Lichtwark
   fixed 5km/h; Rubenson "preferred" speed (not stated in m/s). Subject2/walking1's own rough
   whole-trial pelvis-translation speed proxy is **~1.06 m/s** (a crude COM-displacement/duration
   estimate, NOT a validated stride-based measurement) — closest to/between Farris's 0.75 and 1.25
   m/s conditions. No claim of a speed-matched comparison.
5. **Single subject (subject2), single trial (walking1), right leg only** — same scope caveat as
   every other cert in this family (`docs/MECHANISM_CROSS_SUBJECT.md`); whether the same
   weak-decoupling/long-absolute-length pattern recurs in other subjects is untested here.
6. **Lichtwark & Wilson 2006's body-text mm/degree numbers were not independently re-fetched by this
   session directly** (ScienceDirect returned HTTP 403) — sourced from a sub-agent's full-text
   extraction, internally cross-checked across 3 of its own passes (1 discrepancy caught and
   corrected), while the PMID/DOI/abstract identity and core qualitative mechanism claims WERE
   independently re-confirmed live by this session via NCBI eutils.
7. **Rubenson et al. 2012 (soleus anchor) is single-pass extraction only**, with an unreconciled
   minor internal inconsistency in the source paper's own reported numbers (0.94−0.74=0.20 L0≈7.5mm
   vs the directly-quoted 12.1mm) — both numbers reported as-is, not adjudicated, per LEAN discipline
   on a secondary anchor.
8. **No claim that the absolute fiber-length-scale discrepancy is CAUSED specifically by
   optimal_fiber_length mis-scaling** during subject-specific OpenSim scaling — a plausible, named,
   but not independently verified mechanism here; would require comparing subject2's own imaging/
   cadaver data against the generic LaiArnoldModified2017 template's scaling factors, out of scope
   for this session.
9. **Cronin & Lichtwark 2013 contributed no quantitative content** (confirmed live: a pure methods/
   scope review abstract, fully closed-access body text) — used for qualitative context/pointer only,
   per the task's own naming of it, exactly as disclosed rather than silently dropped or
   silently treated as if it had numbers.

## 9. Files

- `scripts/msk/gastroc_fascicle_kinematics.py` (new) — the full, re-runnable, self-contained
  pipeline (OpenSim API: `model.equilibrateMuscles`, `getFiberLength`, `getPennationAngle`,
  `getTendonLength`, `getGeometryPath().getLength`). Run with
  `.venv-msk/bin/python3 scripts/msk/gastroc_fascicle_kinematics.py` (needs OpenSim 4.6, confirmed
  present: `opensim 4.6-2026-06-22-85aaf64`).
- `data/msk_smoketest/subject2_walking1/gastroc_fascicle_kinematics/gastroc_fascicle_kinematics_results.json`
  (new) — every number in this document, machine-written: muscle properties, per-frame time series
  (activation, MTU/fiber/tendon length, pennation, compliant AND rigid reconstructions, equilibrium-
  solve success flags), the geometric identity gate, derived per-muscle statistics, full-stance
  excursions, the forced-adversary-#2 force-match diagnostic, the smoothing-robustness check, the
  full literature-anchor block (PMIDs, DOIs, verbatim quotes/numbers, extraction-confidence labels),
  the pre-registered-threshold verdict computation, and honest gaps.
- Reused, unedited, zero new OpenSim solve: `data/msk_smoketest/subject2_walking1/
  static_optimization/so/walking1_StaticOptimization_{force,activation}.sto`;
  `OpenSimData/Mocap/IK/walking1.mot` (read-only external drive);
  `OpenSimData/Mocap/Model/LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` (read-only
  external drive, loaded fresh, never modified — confirmed via a `find -newer` self-audit that only
  OpenSim's own `opensim.log` (its C++ logger duplicating stdout) changed on the external drive, the
  same pre-existing, already-disclosed behavior other docs in this repo have found).
- External literature verified live via `curl` to `eutils.ncbi.nlm.nih.gov` (esummary/efetch), both
  by a `watertight-researcher` sub-agent and independently re-checked directly by this session:
  PMID 22219360 (Farris & Sawicki 2012, CORRECTED from the task's stated 22233781), PMID 17050853
  (Lichtwark & Wilson 2006), PMID 16364330 (Lichtwark/Bougoulias/Wilson 2007, already used elsewhere
  in this repo), PMID 22771749 (Rubenson et al. 2012, soleus), PMID 22910172 (Cronin & Lichtwark
  2013, confirmed non-quantitative).
- **Not touched**: `scripts/msk/validate_joint_force.py`, `scripts/msk/static_opt_knee.py`,
  `scripts/msk/tendon_elastic_energy.py` (all reused unchanged via import, per this repo's
  wrap-don't-edit convention), `scripts/msk/contact_waveform_analysis.py`,
  `scripts/msk/contact_muscle_decomp.py`, `scripts/msk/pushoff_plantarflexor_achilles.py`.

No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task). All new
files are untracked, for the coordinator to commit.
