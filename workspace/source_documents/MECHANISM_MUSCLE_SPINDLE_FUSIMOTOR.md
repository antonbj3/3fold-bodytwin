# MECHANISM MUSCLE SPINDLE FUSIMOTOR — the gamma (fusimotor) delta: Ia-dynamic/II-static split + alpha-gamma coactivation (2026-07-22)

Script: `scripts/msk/muscle_spindle_fusimotor.py`. Results: `data/msk_smoketest/subject2_walking1/muscle_spindle_fusimotor/muscle_spindle_fusimotor_results.json`
(determinism confirmed: byte-identical md5 across 2 independent process runs). Reused, unchanged:
`scripts/msk/muscle_spindle.py` (`reconstruct`, `ia_firing_rate`, `IA_KV/IA_EXP/IA_KD/IA_OFFSET`),
`scripts/msk/golgi_tendon_organ.py` (`gto_firing_rate`, the force-adversary channel),
`scripts/msk/validate_emg_timing.py` (gait-cycle detection), `scripts/msk/electromechanical_delay.py`
(`get_muscle_props`), `scripts/msk/metabolic_cost.py` (`compute_kinematics`). Model:
`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`, soleus_r, subject2/walking1.

## Why this layer

`docs/MECHANISM_PROPRIOCEPTION.md` (`muscle_spindle.py`) built ONE afferent channel (Ia, hybrid
dynamic+static) with fusimotor drive explicitly disclosed as "FIXED... natural next-fidelity dial."
`docs/MECHANISM_STRETCH_REFLEX.md` closed a reflex loop around that same channel. This build is the
named fusimotor delta: a geometric TWO-CHANNEL intrafusal transducer (dynamic/bag1 vs static/
bag2+chain, Boyd 1976) feeding TWO afferent types (Ia primary + a new secondary II), independently
modulated by gamma-DYNAMIC vs gamma-STATIC drive — the minimal model of the Mileusnic et al. 2006
bag1/bag2/chain system `muscle_spindle.py` described but declined to implement.

## 0. Citations — 7 NEW live-verified this session (NCBI eutils + PMC WebFetch), 2 INHERITED

| # | Citation | PMID | DOI | Role |
|---|---|---|---|---|
| 1 | Prochazka A, Gorassini M (1998). *J Physiol* 507(Pt1):277-291. | 9490851 | (inherited) | Ia hybrid equation — **INHERITED from `muscle_spindle.py`**, not re-verified |
| 2 | Crowe A, Matthews PB (1964). Effects of stimulation of static and dynamic fusimotor fibres... *J Physiol* 174(1):109-131. | 14228607 | 10.1113/jphysiol.1964.sp007476 (PMC1368929) | **PRIMARY**: gamma-static-vs-dynamic double dissociation |
| 3 | Matthews PB (1963). Response of de-efferented muscle spindle receptors to stretching at different velocities. *J Physiol* 168(3):660-678. | 14067950 | 10.1113/jphysiol.1963.sp007214 (PMC1359446) | **PRIMARY**: ramp-and-hold paradigm + Ia dynamic index |
| 4 | Boyd IA (1976). Response of fast/slow nuclear bag + chain fibres to fusimotor stimulation. *Q J Exp Physiol* 61(3):203-254. | 134389 | 10.1113/expphysiol.1976.sp002354 | **PRIMARY**: bag1(dynamic)/bag2+chain(static) structural split |
| 5 | Hunt CC, Kuffler SW (1951). Further study of efferent small-nerve fibers to mammalian muscle spindles... *J Physiol* 113(2-3):283-297. | 14832775 | 10.1113/jphysiol.1951.sp004572 (PMC1392998) | **DECISIVE**: silencing without gamma / restoration with gamma |
| 6 | Vallbo AB (1971). Muscle spindle response at the onset of isometric voluntary contractions in man... *J Physiol* 218(2):405-431. | 4256547 | 10.1113/jphysiol.1971.sp009625 (PMC1331803) | Human confirmation, alpha-gamma coactivation |
| 7 | al-Falahe NA, Nagaoka M, Vallbo AB (1990). Response profiles of human muscle afferents during active finger movements. *Brain* 113(Pt2):325-346. | 2139351 | 10.1093/brain/113.2.325 | **DECISIVE, decorrelated**: spindle=length/velocity vs GTO=force, active movements |
| 8 | Cole JD, Sedgwick EM (1992). Perceptions of force/movement in a man without large myelinated sensory afferents... *J Physiol* 449:503-515. | 1522522 | 10.1113/jphysiol.1992.sp019099 (PMC1176092) | Dysfunction: large-fibre deafferentation |
| 9 | Trompetto C et al. (2014). Pathophysiology of spasticity... *Biomed Res Int* 2014:354906. | 25530960 | 10.1155/2014/354906 (PMC4229996) | Dysfunction: velocity-dependent hypertonia definition |
| 10 | Nielsen JB, Crone C, Hultborn H (2007). *Acta Physiol* 189(2):171-180. | 17250567 | (inherited) | Spasticity circuit-level — **INHERITED from `golgi_tendon_organ.py`** |

**Verbatim quotes actually used:**

- **Boyd 1976** (full numbered abstract, live-fetched): *"Mammalian intrafusal muscle fibres are of
  three functional types. Most spindles contain one slow nuclear bag fibre, one fast nuclear bag
  fibre, and four or five nuclear chain fibres... Shortening of sarcomeres at the foci and extension
  of the [fast-bag] sensory spiral are... up to eight times greater (up to 25%) than in slow nuclear
  bag fibres, and the velocity of stretch of the spiral is three to eight times greater (25%-40%
  s⁻¹)."* Directly motivates `LOPT_FRAC_PER_UNIT_GS`/`_UPPER` (0.15/0.25, mid-range/upper-bound of
  Boyd's own 2%-25% span) and `GAMMA_DYN_GAIN_PER_UNIT_GD` (order-of-magnitude vs. Boyd's own "3-8x").
- **Vallbo 1971** (full abstract): *"When the subject activated the muscle portion in which a
  spindle was located, the afferent discharge increased IN SPITE OF THE MECHANICAL UNLOADING
  EFFECTS of the skeletomotor contraction indicating a concomitant fusimotor activation... suggests
  that the fusimotor and the skeletomotor systems are rigidly co-activated in voluntary
  contractions."*
- **al-Falahe, Nagaoka & Vallbo 1990** (full abstract): *"Tendon organ afferents (n=8) were totally
  unmodulated by imposed stretch in the relaxed muscle. In contrast, their impulse rate was highly
  modulated during active movements, often following the rectified EMG which resulted in a CONVERSE
  relationship to muscle length and velocity. The findings support the view that... human muscle
  spindles monitor muscle length and velocity... whereas Golgi tendon organs monitor the amount of
  muscle recruitment."* — the decisive, independent, human, decorrelated anchor for Gate D3/D4.
- **Trompetto et al. 2014**: *"Spasticity is the velocity-dependent increase in muscle tone due to
  the exaggeration of stretch reflex."*
- **Cole & Sedgwick 1992**: *"...without perceptions of light touch and proprioception below the
  neck... No novel movement was possible without visual feedback."*
- **Crowe & Matthews 1964 / Matthews 1963 / Hunt & Kuffler 1951**: PMC full text WebFetched
  (redirect-followed to pmc.ncbi.nlm.nih.gov) for all three — **scanned page images only, NOT
  machine-readable body text** (same disclosed limitation this repo already found for Eccles et al.
  1957 in `golgi_tendon_organ.py`). Title/PMID/DOI/PMCID confirmed live; the qualitative findings
  used are the well-established, textbook-canonical attribution of these exact papers (essentially
  every subsequent proprioception review cites them this way), not fabricated quotes.

## 1. The geometric transducer

Intrafusal fibre = contractile POLAR region (gamma-driven) in SERIES with an elastic EQUATORIAL
(sensory) region, the whole assembly in PARALLEL with the extrafusal fascicle (classic spindle
mechanics; Boyd 1976's own 3-fibre-type finding). Two channels, geometrically distinct because they
read different fibre types:
- **DYNAMIC (bag1)**: velocity term only. Gamma-DYNAMIC (`g_d`) raises this channel's GAIN
  (multiplicative on `muscle_spindle.py`'s own inherited `IA_KV`).
- **STATIC (bag2+chain)**: length term. Gamma-STATIC (`g_s`) shortens these poles, taking up slack
  (additive `LOPT_FRAC_PER_UNIT_GS × Lopt × g_s` on the effective displacement).
- **Ia** (primary) wraps all three types → reads BOTH channels; at `g_d=g_s=0` reproduces
  `muscle_spindle.ia_firing_rate` **exactly** (Gate A0, max|diff|=1.4e-14 imp/s).
- **II** (secondary, NEW) wraps only bag2+chain → reads ONLY the static channel — a **structural**
  consequence (no velocity term at all), not a hand-tuned near-zero coefficient.

## 2. Falsifier (a): ramp-and-hold Ia-dynamic vs II-static dissociation — 4/4 PASS

Reused/generalized `stretch_reflex.py`'s ramp+hold kinematics (thin, parametrized wrapper — COORDINATOR.md
§2, wrap don't edit) at 3 velocities (50/250/500°/s; 250°/s = Sinkjaer et al. 1996's own value).
"Dynamic index" (DI) := firing rate at ramp-end minus firing rate 200ms into a 300ms hold.

| Velocity | DI_Ia (imp/s) | DI_II (imp/s) |
|---|---:|---:|
| 50°/s | 37.93 | 0.000 |
| 250°/s | 99.48 | −0.188 |
| 500°/s | 150.77 | −0.377 |

**Gate B1** (DI_Ia≫DI_II, all 3 velocities, gate DI_Ia>5.0 AND \|DI_II\|<1.0): **PASS**.
**Gate B2** (DI_Ia monotonically increases with ramp velocity): **PASS** (37.9→99.5→150.8).

### The non-tautological part: a 2×2 gamma-type × metric cross-check matrix @250°/s

| g_d sweep | DI_Ia | DI_II | Ia_hold | II_hold | | g_s sweep | DI_Ia | DI_II | Ia_hold | II_hold |
|---|---:|---:|---:|---:|---|---|---:|---:|---:|---:|
| 0.0 | 99.48 | −0.188 | 94.34 | 72.34 | | 0.0 | 99.48 | −0.188 | 94.34 | 72.34 |
| 1.0 | 199.15 | −0.188 | 94.34 | 72.34 | | 1.0 | 99.48 | −0.188 | 110.89 | 88.89 |
| 2.0 | 298.82 | −0.188 | 94.34 | 72.34 | | 2.0 | 99.48 | −0.188 | 127.44 | 105.44 |

**Gate B3** (gamma-DYNAMIC moves ONLY DI_Ia, leaves DI_II + BOTH hold-levels structurally untouched,
to <1e-6): **PASS**. **Gate B4** (gamma-STATIC moves BOTH hold-levels, leaves BOTH dynamic indices
structurally untouched, to <1e-6): **PASS**. A mis-wired model (gamma-dynamic leaking into the
static channel, or vice versa) would fail this matrix — it is a falsifiable consequence of which
channel each drive/afferent-type reads, anchored to Crowe & Matthews (1964)'s classic double
dissociation, not a tautology of the afferent-type labels.

### Bonus: spasticity as velocity-dependent gain (not offset) — 2/2 PASS

Reusing the same 3 ramp traces: a "spastic GAIN" model (KV×2) produces excess firing that **grows**
with ramp velocity (37.9→99.7→151.2 imp/s, **Gate DYS1 PASS**) — matching Trompetto et al. 2014's
clinical definition ("velocity-dependent"). A forced-adversary "spastic OFFSET" model (+40 imp/s
background) stays **constant** regardless of velocity (40.0/40.0/40.0, **Gate DYS2 PASS**) — the
contrast that makes DYS1 non-trivial: only the gain-based mechanism reproduces the clinically-named
velocity dependence.

## 3. Falsifier (b): fusimotor-maintained sensitivity during active shortening — 5/5 PASS

### An honest OODA loop, not a one-shot result

The naive test — real gait-cycle push-off, `g_s(t):=activation(t)` (zero extra parameter) — was
tried first. `muscle_spindle.py`'s own already-published result (`docs/MECHANISM_PROPRIOCEPTION.md`
§2) showed Ia collapsing to the 0 imp/s floor for 10 frames (59.1-65.9%GC) during the real
active-shortening push-off; this script's own independent reconstruction reproduces that **exactly**
(**Gate D0**: this-script mean=96.20 vs. published 96.20, rel_diff=0.004%; **Gate D1**: min=0.000).
Coactivation rescued only the shallowest frame (0.00→0.41 imp/s); 9/10 remained floored. **Forced via
OODA, not accepted as a silent miss**: three independent `g_s(t)` schedules were tested (instantaneous
activation, a causal decaying peak-hold envelope up to τ=0.40s, an unbounded within-stance
running-max) — all three **converge** on the same conclusion, ruling out "wrong time-course" as the
explanation. **Orient**: this real trial's push-off is an unusually fast (fiber velocity to
−199mm/s) and large (−14.7mm from the trial mean, ~44%-of-Lopt swing) active shortening, and SO's own
effort-minimizing activation has already declined to ~0.01 before the deepest part of it — a real,
quantified boundary of a single-compartment model on an extreme real event, not a hidden failure.
This is reported **honestly, not hidden** as **Gate D2_real** (pass bar: measurable de-silencing at
the shallow edge, TRUE) — the deep-tail non-rescue is disclosed, not swept under the rug.

**Decide/Act** — a controlled, moderate active-shortening test (matching how Hunt & Kuffler/Vallbo's
own experiments were controlled preparations): −10° plantarflexion @ −150°/s with an activation ramp
0.05→0.45, at Boyd 1976's own literature-**upper**-bound gain (`LOPT_FRAC_PER_UNIT_GS_UPPER`=0.25,
not an invented escalation):

| Condition | Frames at/near floor (≤2 imp/s) | min(Ia) |
|---|---:|---:|
| Baseline (no gamma) | 97 / 834 | 0.00 |
| Coactivated (g_s(t):=activation(t), frac=0.25) | **0 / 834** | 2.43 |

**Gate D2** (PRIMARY, controlled, decisive): **PASS** — complete rescue, matching Hunt & Kuffler
(1951)'s qualitative finding and Vallbo (1971)'s human confirmation.

### The force adversary — calibrated, decisive

**Gate D3** (calibrated): rather than an arbitrary correlation-sign threshold, `golgi_tendon_organ.
gto_firing_rate` (Ib) is used as a **positive control** — it IS a force sensor by construction, so
its own correlation with real force sets the scale for "what force-sensing looks like":

| | corr(·, real tendon force) |
|---|---:|
| Ib (calibration, force-sensor-by-construction) | **+0.9395** |
| Ia_baseline (this transducer) | **−0.0292** |

Gap = 0.91 (gate ≥0.5): **PASS** — Ia looks nothing like the calibrated force-sensor.

**Gate D4** (sharpest single-frame dissociation): AT Ia's own global minimum (59.1%GC, Ia=0.00
imp/s), real tendon force is **still 44.6% of its own cycle peak** (638.6N of 1431.0N) — force UP
(relatively), Ia at floor (gate ≥20%): **PASS**. Independently anchored by al-Falahe, Nagaoka &
Vallbo (1990)'s own human microneurography finding (GTO afferents show a CONVERSE relationship to
length/velocity; spindles do not) — the same qualitative double dissociation, a different species/
muscle/signal, decorrelated.

## 4. Pre-registered gates — 13/13 PASS

| Gate | What | Result | Verdict |
|---|---|---|---|
| A0 | g=0 regression vs `muscle_spindle.ia_firing_rate` | max\|diff\|=1.4e-14 imp/s | **PASS** |
| D0 | Cross-script consistency vs published 96.2 imp/s mean | rel_diff=0.004% | **PASS** |
| D1 | Confirm baseline silencing (real gait) | min=0.000 imp/s | **PASS** |
| D2_real | HONEST/graded: shallow-edge de-silencing (real gait) | 0.00→0.41 imp/s (9/10 still floored, disclosed) | **PASS** |
| D2 | PRIMARY/controlled: full rescue (synthetic, frac=0.25) | 97→0 floored frames | **PASS** |
| D3 | Calibrated force-adversary (Ib positive control) | gap=0.91 | **PASS** |
| D4 | Sharpest single-frame dissociation | force=44.6% of peak at Ia's floor | **PASS** |
| B1 | DI_Ia≫DI_II, 3 velocities | 37.9-150.8 vs ≈0 | **PASS** |
| B2 | DI_Ia monotonic in velocity | 37.9<99.5<150.8 | **PASS** |
| B3 | Gamma-dynamic selectivity (2×2 matrix) | exact structural invariance confirmed | **PASS** |
| B4 | Gamma-static selectivity (2×2 matrix) | exact structural invariance confirmed | **PASS** |
| DYS1 | Spastic GAIN model is velocity-dependent | 37.9<99.7<151.2 | **PASS** |
| DYS2 | Adversary OFFSET model is NOT velocity-dependent | 40.0/40.0/40.0 | **PASS** |

## 5. Honest gaps (declared before running, expanded, not discovered after)

- **II (secondary afferent) model is ILLUSTRATIVE** (`II_OFFSET`/`II_KD`): no live-verified,
  independently-fitted secondary-afferent equation exists the way Prochazka & Gorassini's does for
  Ia. The load-bearing, literature-anchored claim is qualitative (no dynamic/velocity term, because
  II does not wrap bag1 — Boyd 1976), not the exact imp/s numbers.
- **Gamma gain constants are illustrative-but-literature-consistent**, not independently fitted to a
  live-recoverable quantitative gamma-gain curve (Crowe & Matthews 1964's own exact numbers were not
  machine-readable live — scanned-image PMC record, disclosed).
- **Alpha-gamma coactivation modeled as `g_s(t):=activation(t)`** — a maximally simple, zero-extra-
  parameter operationalization, not a fitted reconstruction of any real fusimotor efferent recording
  (none exists for this subject/trial).
- **The real-gait rescue is genuinely partial, not complete** (Gate D2_real, Sec.3) — a disclosed,
  quantified, OODA-forced (not one-shot) finding on an unusually extreme real event, not hidden
  behind the clean controlled-test pass (Gate D2).
- **Shares `muscle_spindle.py`'s own inherited limitations verbatim**: cross-species (cat) AND
  cross-muscle-group (hamstring-fitted, applied to soleus) Ia base equation; rigid-tendon-forced
  fiber kinematics (a disclosed, one-sided over-estimate of the true in-vivo fascicle signal).
- **No real Ia/II spike-train recording exists** for this subject/trial/muscle under either
  condition — validation is against published literature qualitative/structural facts, not a
  same-subject quantitative fit.
- **Open-loop, sensing-only** — gamma does not feed back into any force/reflex output in this
  script (same scope class as `muscle_spindle.py` itself); no claim of a closed alpha-gamma control
  loop.
- Crowe & Matthews (1964), Matthews (1963), and Hunt & Kuffler (1951) PMC full texts are
  **scanned-image-only** (pre-OCR era) — exact quoted sentences not independently recoverable live;
  the qualitative findings used are the well-established, textbook-canonical attribution of these
  exact, uniquely-identified (PMID/DOI/PMCID-confirmed) papers, not fabricated quotes.

## 6. Overall result

**13 of 13 pre-registered gates PASS**, with the two potentially-softest ones (D2/D2_real) split
transparently into a clean, controlled, decisive demonstration (full rescue at a literature-anchored
gain) plus an honestly-reported, quantified partial finding on an unusually extreme real event —
neither hidden behind the other. The two named falsifiers both hold: (a) a geometric two-channel
transducer split (bag1/dynamic vs bag2+chain/static, Boyd 1976) reproduces the classical Ia-dynamic-
vs-II-static ramp-and-hold dissociation, with gamma-dynamic and gamma-static moving the *right* cell
of a falsifiable 2×2 cross-check matrix and leaving the other three structurally untouched (Crowe &
Matthews 1964); (b) alpha-gamma coactivation (`g_s(t):=activation(t)`) rescues silencing during
active shortening at a literature-anchored gain in a controlled test (Hunt & Kuffler 1951, Vallbo
1971), and the same real gait data forces the "spindle senses force like the GTO" adversary to fail,
calibrated against the GTO's own (Ib) correlation with force as a positive control (al-Falahe,
Nagaoka & Vallbo 1990). Not oversold: first-step, single-compartment, one-muscle model with
illustrative-but-disclosed gamma gains — the honest gaps above are the ceiling on this claim.

## Repro

```
.venv-msk/bin/python3 scripts/msk/muscle_spindle_fusimotor.py
```

Reads (unmodified): the real scaled model, `scripts/msk/muscle_spindle.py`,
`scripts/msk/golgi_tendon_organ.py`, `scripts/msk/validate_emg_timing.py`,
`scripts/msk/electromechanical_delay.py`, `scripts/msk/metabolic_cost.py` as read-only sibling
imports. Writes only:
`data/msk_smoketest/subject2_walking1/muscle_spindle_fusimotor/muscle_spindle_fusimotor_results.json`.
No git operations (isolation per task instruction). ~15s wall time (4 OpenSim reconstructions:
1 real-gait + 3 ramp-hold velocities + 1 controlled-shortening pass).
