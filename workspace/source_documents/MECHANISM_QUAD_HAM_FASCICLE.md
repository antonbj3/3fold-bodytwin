# MECHANISM QUADRICEPS/HAMSTRING FASCICLE KINEMATICS — extending muscle-INTERNAL validation from the ankle to the knee prime movers (2026-07-21)

Extends `docs/MECHANISM_GASTROC_FASCICLE.md` (gastrocnemius/soleus fascicle kinematics vs in-vivo
B-mode ultrasound) to the **knee prime movers**: vastus lateralis (`vaslat_r`, quadriceps) and
biceps femoris (`bflh_r` long head + `bfsh_r` short head, hamstrings). Reuses the IDENTICAL method
(compliant `Model.equilibrateMuscles` reconstruction + rigid-tendon forced-adversary bookend +
geometric self-consistency identity) on a mechanically DIFFERENT muscle-tendon architecture (shorter,
stiffer patellar/hamstring tendons vs the long, highly-compliant Achilles) to test whether gastroc's
already-documented systematic pattern (decoupling weaker than the literature's qualitative picture;
absolute fiber length running long) is Achilles/ankle-specific or a model-wide issue.

**Headline: the pattern RECURS, and this session found its externally-anchored, machine-verified
ROOT CAUSE.** A newly-verified single-source cadaveric-architecture reference (Ward, Eng, Smallwood &
Lieber 2009, *Clin Orthop Relat Res*, PMID 18972175, open-access PMC full text) that happens to cover
all SIX muscles examined across both docs shows this model's own `optimal_fiber_length` parameter
runs **+0.82 to +2.92 SD above** a modern, high-resolution (21-specimen) cadaver reference for EVERY
SINGLE muscle checked (vaslat_r +2.39, bflh_r +0.82, bfsh_r +1.11, gasmed_r +2.34, gaslat_r +2.92,
soleus_r +1.13 — reusing the sibling gastroc cert's own already-certified values, zero re-solve) while
**pennation angle stays within ±1 SD for all six** (−0.64 to +0.79) — a precise, geometrically clean
diagnosis that this model's fiber-length SCALING specifically (not its pennation) runs systematically
long, with severity that tracks the muscle groups already flagged (gastrocs+VL worst, BF+soleus
milder). Dynamically, this propagates exactly as predicted: `vaslat_r`'s stance-phase fiber range
(101.46–141.03mm) sits clearly above the task's own stated ~80–100mm envelope, and its
loading-response decoupling ratio (**0.777** compliant vs 0.946 rigid-bookend — a real, non-tautological
~17-point tendon-compliance effect, but far from Chleboun et al. 2007's own reported *statistically
flat* fascicle length, P=.22, in the exact same phase) directly **contradicts** the literature's
qualitative mechanism. `bflh_r` shows an even WEAKER tendon-buffering effect (0.931 compliant vs 0.981
rigid, only a ~5-point gap) and **crosses this cert family's own pre-registered ≥0.80 hard
"no-decoupling" falsifier** — matching soleus_r's own hard-fail in the gastroc doc. The ONE clean
qualitative TIMING mechanism tested for the knee muscles — the well-established biarticular-hamstring
"stretched-at-terminal-swing" signature — DOES reproduce correctly (global fiber-length maximum at
96.8%GC, robust to smoothing, falls inside the pre-registered wrapped window) — this is not a
wholesale failure, timing mechanisms work, and the compliant-vs-SO-force match is actually BETTER for
the knee muscles (13.6–32.5% RMS relerr) than gastroc/soleus (63–118%) — but the **fiber-length-SCALE
and tendon-BUFFERING-STRENGTH gap is now confirmed across all 5–6 muscles checked in this cert
family**, with a real, quantified, externally-anchored root cause, not just a repeated symptom. No
in-vivo biceps-femoris-during-WALKING measurement exists anywhere in the literature (confirmed by an
exhaustive dual search — this session directly, plus an independently-dispatched sub-agent, plus a
cross-confirming unrelated third-party 2022 paper hitting the identical gap, §4.3) — but a REAL
in-vivo B-mode ultrasound dataset for the same muscle during RUNNING (Mao et al. 2024, independently
re-verified word-for-word by this session from open-access full text) reports fascicle lengthening
is only **~30%** of MTU length change in the analogous eccentric late-swing phase — a real,
externally-measured decoupling ratio roughly **3× lower** than this model's own 0.931 for `bflh_r`,
upgrading the decoupling finding from a model-internal threshold comparison to an externally
corroborated one, in the same direction the Chleboun-anchored `vaslat_r` finding already showed.

## Headline numbers

| muscle | role | key window | fiber_pp (mm) | MTU_pp (mm) | decoupling (compliant/rigid) | full-cycle fiber range (mm) | OFL z-score vs Ward&Lieber 2009 |
|---|---|---|---:|---:|---:|---:|---:|
| vaslat_r | quad (VL) | loading-response 0–23%GC | 10.79 | 13.89 | **0.777** / 0.946 | 101.46–150.59 | **+2.39 SD** |
| bflh_r | hamstring (BF long head) | terminal-swing/early-stance 75–20%GC (wrapped) | 40.38 | 43.37 | **0.931** / 0.981 | 80.40–127.59 | +0.82 SD |
| bfsh_r | hamstring (BF short head, mono-articular, bonus/descriptive) | full cycle | 27.06 | 28.43 | 0.952 / 0.967 | 123.51–150.57 | +1.11 SD |
| *gasmed_r (context)* | *plantarflexor* | *0–40%GC (gastroc doc)* | — | — | *0.629 / 0.989* | *50.92–85.16* | *+2.34 SD* |
| *gaslat_r (context)* | *plantarflexor* | *0–40%GC (gastroc doc)* | — | — | *0.557 / 0.981* | *60.78–94.66* | *+2.92 SD* |
| *soleus_r (context)* | *plantarflexor* | *0–40%GC (gastroc doc)* | — | — | *0.868 / 0.937* | *43.62–64.70* | *+1.13 SD* |

All new numbers (vaslat_r/bflh_r/bfsh_r rows + all 6 OFL z-scores) machine-computed by
`scripts/msk/quad_ham_fascicle_kinematics.py`, written to
`data/msk_smoketest/subject2_walking1/quad_ham_fascicle_kinematics/quad_ham_fascicle_kinematics_results.json`.
Italicized rows reused verbatim (zero re-solve) from the sibling
`gastroc_fascicle_kinematics_results.json` for direct comparison. Model:
`Millard2012EquilibriumMuscle`, `ignore_tendon_compliance=false` for all 3 new muscles, confirmed
live from the OpenSim API (Step 0).

## 1. Pre-registration (stated before this session's OODA-forced smoothing-bug correction, §2.4, changed any number)

**C / ¬C, thresholds fixed in advance, directly reusing this task's own framing ("do VL/BF absolute
fascicle lengths + their stance-phase length-change fall within the measured ultrasound envelope, or
do they run long / mis-time like gastroc did"):**

- **C (fascicle kinematics validated):** absolute fiber length stays within a geometry-derived outer
  band (`[0.25,1.90]×optimal_fiber_length`, the model's OWN `ActiveForceLengthCurve` support,
  verified live, identical for every muscle checked) AND a tighter soft pass bar
  (`[0.50,1.50]×optimal_fiber_length`, mirroring gastroc's own 0.50 floor with a symmetric 1.50
  ceiling since VL/BF have real lengthening-dominant events gastroc's push-off-only test never
  needed) AND the muscle-specific key-event timing (vaslat_r: near-isometric loading-response
  decoupling ≤0.50; bflh_r: peak-stretch lands inside the terminal-swing/early-stance window,
  robust to smoothing) AND passes the geometric self-consistency identity.
- **¬C (fiber-kinematics ERROR, the SAME systematic issue gastroc found, now in a different
  muscle-tendon architecture):** over-shortening/over-lengthening past the soft bar, OR mis-timing
  (bflh_r's peak-stretch lands outside its window), OR — reusing gastroc's own additional
  pre-registered axis VERBATIM, not re-derived, for direct cross-muscle-family comparability — no
  decoupling (fiber/MTU peak-to-peak ratio ≥0.80 in the muscle's own key window).
- **Forced adversary (the one this task is tempted to skip):** near-isometric-then-lengthen/shorten
  behavior could be trivial for ANY reasonably-parametrized compliant Hill muscle. Forced via the
  IDENTICAL rigid-tendon geometric bookend gastroc used (§2.2): if the rigid bookend showed the SAME
  qualitative pattern, the compliant result would be a test-construction artifact. It does not
  (rigid decoupling 0.946–0.981 vs compliant 0.777–0.952, all 3 muscles) — the adversary FALLS,
  proving the (weak) decoupling that does exist is caused by real tendon compliance, not tautology.
- **A second forced adversary, identified independently this session (not present in gastroc,
  because gastroc's stance-only scope never needed it):** is "the model's fiber length runs long"
  just this session's own impression, or does it hold up against an independent, high-quality,
  EXTERNAL reference covering ALL the muscles in this cert family? Forced via Ward, Eng, Smallwood &
  Lieber 2009's cadaveric architecture data (§4.2) — it holds up: all 6 muscles' `optimal_fiber_length`
  parameters exceed the reference mean, by an amount (+0.82 to +2.92 SD) that is NOT uniform across
  muscles and correlates with which muscles show the worse dynamic "runs long" symptom.

**Symmetric-QC checks forced, not assumed:** (1) the conversion-identity check (normalized fiber
length × optimal-fiber-length = absolute mm) at all 158 frames × 3 muscles: max relative error
2.19×10⁻¹⁶ (bfsh_r) to 2.20×10⁻¹⁶ — exact to floating-point precision, PASS (§3). (2) the
compliant-vs-rigid-tendon bookend — BOTH reconstructions computed and reported for all 3 muscles,
not asserted in prose (§6.2). (3) a forced-adversary #2 (does the reconstruction reproduce SO's own
per-muscle force?), reused verbatim from gastroc, run and disclosed (§5) — here it is BETTER
(13.6–32.5% RMS relerr) than gastroc's own (63–118%), an honest counter-example against a falsely
uniform "everything gets worse" narrative. (4) A bonus same-%GC cross-cycle consistency check this
session identified from understanding this specific trial's own %GC-wraparound geometry (§2.5,
§6.5) — vaslat_r borderline FAILS (5.04mm vs the 5.0mm pre-registered threshold), traced (not just
observed) to a genuine small stride-to-stride difference at the recording's own tail frame, not a
pipeline bug (verified by direct frame-pair inspection).

## 2. Method

### 2.1 Primary reconstruction (compliant, model-faithful) — identical to gastroc

Reused, **zero re-solve of Static Optimization itself**: the already-certified subject2/`walking1`
SO activation trajectory and the already-validated IK motion (same files gastroc used). At each of
the 158 frames (now the FULL gait cycle, 0.57–99.81%GC — unlike gastroc's stance-only scope, since
VL/BF's literature-anchored key events include swing phase): set the model's coordinates from IK, set
`vaslat_r`/`bflh_r`/`bfsh_r`'s activation to SO's own value, call `Model.equilibrateMuscles(state)`,
read `getFiberLength`/`getPennationAngle`/`getTendonLength`. 0/158 equilibrium failures.

### 2.2 Forced-adversary bookend (rigid tendon) — identical method to gastroc

`fiber_length_rigid(t) = sqrt(h² + (MTU_length(t) − tendon_slack_length)²)`, computed alongside every
frame, no equilibrium solve.

### 2.3 Muscle-appropriate key-event windows — DERIVED, not copied from gastroc's push-off window

Push-off/take-off is a plantarflexor-specific concept; it does not apply to knee extensors/flexors.
This session derived each muscle's key window from (a) this trial's OWN already-certified SO
activation trace (a legitimate, already-established input, not an outcome being tested) and (b)
textbook gait-phase/EMG-timing convention (Perry & Burnfield gait phases; Winter's EMG-phase
tables) — **before** looking at any fiber-length/literature-comparison outcome:

- **vaslat_r**: SO's own activation peaks at **11.6%GC** (this session's own measurement) — textbook
  loading response (0–15%GC), while the knee is FLEXING (lengthening the knee-extensor MTU) — the
  classic eccentric-quadriceps-control event. Window pre-registered generously as **(0.0, 23.0)%GC**.
  A SECOND VL burst is textbook-expected (and IS reported by Chleboun 2007, §4.1) at terminal swing;
  this model's SO solution does **not** reproduce it — `vaslat_r`/`vasmed_r`/`vasint_r` activation
  stays at the 0.01 floor across 70–100%GC (checked at 0.75%GC resolution, §6.1), while `recfem_r`
  (rectus femoris, a redundant biarticular alternative) peaks at 0.362 @ 24.7%GC and separately at
  0.159 @ 72.7%GC — a real, machine-verified, disclosed redundant-actuator strategy difference from
  Chleboun's real subjects, not silently worked around.
- **bflh_r**: SO's own activation burst spans **~79–100%GC wrapping into ~0–16%GC** (peak 0.099 @
  0.6%GC) — textbook terminal swing (knee extending + hip flexing, BOTH actions simultaneously
  STRETCH this biarticular knee-flexor/hip-extensor — the well-known "hamstrings maximally stretched
  at terminal swing" mechanism from sprint-injury biomechanics, present here at walking-appropriate
  lower magnitude). Independently cross-checked via PURE KINEMATICS (Step-0 exploration, before any
  equilibrium solve): `bflh_r`'s own global MTU-length maximum across the full cycle lands at
  **98.3%GC**, coincident with the activation burst — a geometry-only fact, not assumed. Window
  pre-registered as the WRAPPED range **(75.0, 20.0)%GC** (i.e., ≥75 OR ≤20).
- **bfsh_r** (mono-articular knee flexor): activation is broad/near-continuous (elevated ~30–71%GC
  and again ~79–100%GC wrapping to ~16%GC, one brief dip 72–77%GC) — a materially different temporal
  pattern from `bflh_r`'s sharp single burst, a real functional dissociation this model's SO solution
  produces between the two "biceps femoris" heads. Reported descriptively (full-cycle only, no
  forced literature falsifier) since no literature anchor isolating the short head during gait was
  found this session, consistent with LEAN discipline.
- **Geometric cross-check** (internal-consistency confirmation, not an independent finding):
  `vaslat_r` (knee extensor) and `bfsh_r` (knee flexor) are BOTH pure mono-articular knee muscles;
  their MTU-length global extrema land at the SAME %GC (`vaslat_r` max / `bfsh_r` min at 74.2%GC;
  `vaslat_r` min / `bfsh_r` max at 41.0%GC) — exactly the mirror-image behavior single-joint
  antagonist geometry requires.

### 2.4 Forced-adversary catch, OODA-corrected (smoothing edge-padding artifact)

**Observed:** a first run of the smoothing-robustness check (§6.6) reported `vaslat_r`'s
loading-response fiber peak-to-peak SHIFTING from 10.79mm (raw) to 47.59mm (smoothed) — a suspicious
+36.8mm jump, an order of magnitude larger than any robustness shift this cert family has ever seen
(gastroc's own shifts were ≤1.5%GC / <0.5mm). **Oriented:** direct inspection of the smoothed array
found `smoothed[0] = 65.50mm` vs `raw[0] = 108.55mm` — a ~40% corruption at the array's true first
index, exactly consistent with `np.convolve(y, kernel, mode='same')`'s implicit ZERO-padding at the
boundary (108.55×3/5 = 65.1, matches to 3 significant figures). This is a genuine bug this session
introduced, not present in gastroc (whose push-off window, 35–80%GC, sits far from the true array
edges) — `vaslat_r`'s loading-response window (0–23%GC) uniquely includes array index 0 itself
(t=0 is already at 6.37%GC). **Decided + Acted:** replaced `np.convolve(...,mode='same')` with an
edge-REPLICATING pad + `mode='valid'` (`smooth_edge_safe()`, `scripts/msk/quad_ham_fascicle_kinematics.py`);
re-ran — the shift dropped to a sane **−0.33mm**, consistent with the rest of this cert family.
Caught, diagnosed, fixed, and disclosed here rather than silently patched.

### 2.5 %GC mapping — reused verbatim; a bonus cross-cycle check this specific mapping enables

Reused VERBATIM from `docs/MECHANISM_CONTACT_WAVEFORM.md` (`T_stride=1.327s`,
`R_heel_strike_t=1.2425s`). Because this trial's recording (t=0–1.57s) spans MORE than one stride
(1.327s) starting MID-cycle (heel-strike at t=1.2425s falls inside the window), the %GC range
[~6.4%,~24.7%] is sampled at two distinct real times: the tail of the recording (t<0.5s, wrapping to
low %GC) and the head of the fully-captured stride (t>1.0s). This session exploited that as a free
test-retest internal-consistency check (§6.5) — not present in gastroc, which never needed to leave
its stance-only window near this seam.

## 3. Verification gates (machine-checked, all reported — PASS unless noted)

| gate | result |
|---|---|
| Time-grid alignment (IK vs SO-force vs SO-activation) | **PASS** (`np.allclose`) |
| `equilibrateMuscles()` convergence, all 158 frames × 3 muscles | **PASS**, 0/158 failures |
| Geometric self-consistency identity (h = fiber·sin(pennation)) | **PASS**, max relerr 0.0% (all 3, floating-point exact) |
| Conversion-identity (normalized fiber × OFL = absolute mm) | **PASS**, max relerr 2.19–2.20×10⁻¹⁶ |
| AFLC curve-support hard degeneracy gate `[0.25,1.90]×OFL` | **PASS**, all 3 muscles, all 158 frames |
| Soft pass bar `[0.50,1.50]×OFL` (non-degenerate operating point) | **PASS**, all 3 muscles (norm range 0.675–1.131) |
| Compliant-vs-rigid contrast (proves decoupling is real, not tautological) | **PASS**, rigid 0.946–0.981 vs compliant 0.777–0.952 |
| Decoupling ratio ≤0.50 pass bar (reused verbatim from gastroc) | **FAIL**, all 3 muscles (0.777–0.952) |
| No-decoupling hard falsifier (ratio ≥0.80) | **TRIGGERED** for bflh_r (0.931) and bfsh_r (0.952, descriptive); NOT triggered for vaslat_r (0.777, close but under) |
| bflh_r mis-timing check (peak-stretch inside wrapped window, robust to smoothing) | **PASS**, 96.8%GC raw → 97.6%GC smoothed, still inside |
| Cross-cycle same-%GC consistency (bonus check, <5.0mm) | **PASS** bflh_r (3.87mm)/bfsh_r (2.60mm); **borderline FAIL** vaslat_r (5.04mm, traced to the recording's own last frame) |
| Ward & Lieber 2009 architecture cross-check, pennation (\|z\|<1 SD sanity) | **PASS**, all 6 muscles (−0.64 to +0.79 SD) |
| Ward & Lieber 2009 architecture cross-check, fiber length (\|z\|<1 SD would be the "well-calibrated" bar) | **FAIL**, all 6 muscles (+0.82 to +2.92 SD) |

## 4. External literature anchors — verified live this session

### 4.1 Vastus lateralis — Chleboun, Busic, Graham & Stuckey 2007, JOSPT (PMID 17710906, verified live)

Confirmed via direct `curl` to `eutils.ncbi.nlm.nih.gov` (esearch → esummary → full efetch
abstract), independently, this session — not sourced from a sub-agent. **"Fascicle length change of
the human tibialis anterior and vastus lateralis during walking."** *J Orthop Sports Phys Ther.* 2007
Jul;37(7):372-9. DOI 10.2519/jospt.2007.2440. n=9 healthy subjects, treadmill walking, fascicle
length via real-time B-mode ultrasound + EMG + knee/ankle joint angle, repeated-measures ANOVA.

**Verbatim findings:** *"During the initial portion of stance, when the TA and VL muscles were
active, the ankle plantar flexed and the knee joint flexed, suggesting muscle-tendon complex
lengthening, but the fascicle length of both muscles remained constant (TA, P=.93; VL, P=.22)."*
*"The VL muscle became active again at the end of swing as the knee extended, and the fascicle
length decreased (P<.05)."* *"The lack of change in fascicle length during the initial portions of
stance phase suggests a nearly isometric muscle action of the TA and VL. There is a possible
interaction occurring between the fascicle and tendon... such that the tendon lengthens to allow
joint motion and potentially to store elastic energy."*

**This abstract reports P-values/qualitative direction, NOT absolute mm figures** — full text is
paywalled (confirmed via Unpaywall API this session: `is_oa:false`, no repository copy, no legitimate
open-access route found). The falsifiable claim tested against this model (§6.1) is therefore
QUALITATIVE/DIRECTIONAL (near-isometric fascicle length during the loading-response window while the
MTU lengthens, i.e. a LOW decoupling ratio), matching exactly what this abstract reports — not an
absolute-mm comparison. A companion paper by the same group, PMID 18827324 ("Vastus lateralis
fascicle length changes during stair ascent and descent," 2008, same journal), was also verified live
but is a stair task, not level walking, and was not used as the primary anchor.

### 4.2 Ward, Eng, Smallwood & Lieber 2009, Clin Orthop Relat Res (PMID 18972175, verified live, OPEN ACCESS)

Found via NCBI esearch, confirmed via esummary, and — critically, unlike most anchors in this cert
family — its **full text was directly fetched and read this session** (PMCID PMC2650051, no
paywall). **"Are current measurements of lower extremity muscle architecture accurate?"** 21 human
lower-extremity cadaveric specimens, 27 muscles dissected, fiber length normalized to a reference
sarcomere length of 2.7μm (removing cross-study fixation-length confounds), n=18–20 specimens per
muscle. This is explicitly a MODERN, higher-quality alternative to the older, much smaller,
"most widely used" Wickiewicz et al. dataset the paper's own discussion says many musculoskeletal
models are calibrated against — and the paper's own verbatim finding is that its NEWER data runs
**LONGER**, not shorter, for exactly the muscle groups this doc examines: *"The key differences in
our data are the **longer** fiber lengths of the knee extensors, knee flexors, and ankle
plantarflexors and shorter fiber lengths of the ankle dorsiflexors. Thus, musculoskeletal models
based on previous architectural data are likely in error."* This makes the comparison below
conservative (harder to trigger, not a cherry-picked low bar): this model's parameters are compared
against the HIGHER of two competing literature conventions, and still run long.

**Table 3 values used (normalized fiber length Lf, mean±SD; pennation angle, mean±SD), hand-transcribed
from the verified full-text fetch:**

| muscle | Lf mean±SD (mm) | pennation mean±SD (deg) | n |
|---|---:|---:|---:|
| Vastus lateralis | 99.4 ± 17.6 | 18.4 ± 6.8 | 19 |
| Biceps femoris long head | 97.6 ± 26.2 | 11.6 ± 5.5 | 18 |
| Biceps femoris short head | 110.3 ± 20.6 | 12.3 ± 3.6 | 19 |
| Gastrocnemius medial head | 51.0 ± 9.8 | 9.9 ± 4.4 | 20 |
| Gastrocnemius lateral head | 58.8 ± 9.5 | 12.0 ± 3.1 | 20 |
| Soleus | 44.0 ± 9.9 | 28.3 ± 10.1 | 19 |

**This model vs Ward & Lieber, z-scores (machine-computed, §8 of the pipeline):**

| muscle | model OFL (mm) | z (fiber length) | model pennation (deg) | z (pennation) |
|---|---:|---:|---:|---:|
| vaslat_r | 141.43 | **+2.39** | 14.49 | −0.58 |
| bflh_r | 119.18 | +0.82 | 10.08 | −0.28 |
| bfsh_r | 133.18 | +1.11 | 15.14 | +0.79 |
| gasmed_r *(context)* | 73.98 | *+2.34* | 9.49 | *−0.09* |
| gaslat_r *(context)* | 86.52 | *+2.92* | 12.05 | *+0.01* |
| soleus_r *(context)* | 55.16 | *+1.13* | 21.85 | *−0.64* |

**Interpretation, geometrically precise:** pennation is well-calibrated for every muscle in this
cert family (\|z\|≤0.79); fiber LENGTH specifically runs long for every muscle, by an amount that
tracks the already-observed severity ordering (gastrocs + vaslat_r worst at +2.3 to +2.9 SD; bflh_r/
bfsh_r/soleus_r milder at +0.8 to +1.1 SD) — a precise, single-source, externally-anchored root-cause
diagnosis the gastroc doc itself could only call "plausible... not independently verified" (its §6.4,
§8.8). Framed cautiously: with only reported summary statistics (not raw per-specimen data) from an
n=18–20 sample of unknown distributional shape, "SD above the mean" is a plausibility/percentile
framing, not a formal hypothesis-test p-value.

### 4.3 Biceps femoris during WALKING — an extensive, exhaustive, independently-cross-confirmed negative search, and a real decorrelated fallback anchor

The task's own starting hypothesis named "Kellis" and ~90–110mm for biceps femoris. This session
searched hard, directly, live (NCBI esearch, ~10 query variants combining `Kellis`/`biceps femoris`+
`fascicle`+`walking`/`gait`/`treadmill`/`ultrasound`, and author-agnostic variants substituting
`Silder`/`Schache`/`Thelen`), then dispatched an independent `watertight-researcher` sub-agent for a
SECOND, broader pass (~20 further query variants plus author-targeted sweeps), and found **NO clean
in-vivo B-mode ultrasound (or any other in-vivo modality) biceps-femoris-during-WALKING paper**, by
Kellis or anyone else. This negative is independently corroborated by a THIRD, unrelated party:
**Luis, Afschrift, De Groote & Gutierrez-Farewik 2022** (PMID 36277379, "Evaluation of
musculoskeletal models, scaling methods, and performance criteria for estimating muscle excitations
and fiber lengths across walking speeds," *Front Bioeng Biotechnol* — PMID/title/journal/authors
independently confirmed live via NCBI esummary by this session), doing exactly this kind of
muscle-fiber-length model-validation exercise for a different purpose, whose own Figure 6 lists
in-vivo experimental normalized fiber-length sources for SO/GM/GL/VL (citing Farris & Sawicki 2012,
Chleboun 2007 — the SAME anchors already established in this cert family) but conspicuously has
**no BF or ST entry**, falling back to the SAME Ward & Lieber 2009 cadaveric dataset used here (§4.2)
for BF/ST normalization — independent evidence this is a known, field-wide data gap, not a search
failure specific to this session.

**Symmetric-QC — near-misses checked and explicitly killed, not just successes reported:**
Thelen, Lenz, Francis, Lenhart & Hernández 2013 (PMID 23540723, *J Biomech*, "Empirical assessment
of dynamic hamstring function during human walking") looked promising by title but is
electrically-stimulated MEDIAL hamstrings (semimembranosus/semitendinosus), not BF, and measures
induced joint motion, not fascicle length — **not a match**. Alexander, Schwameder, Baker & Trinler
2021 (PMID 34509042, *Gait Posture*) reports BF long head muscle FORCE model-agreement (the reference model vs
OpenSim), no length value — **not a match**. Riley, Franz, Dicharry & Kerrigan 2010 (PMID 20022251,
*Gait Posture*, verified live — title/journal/authors match) reports BF long head MUSCLE-TENDON-UNIT
length (not fascicle length — a different quantity) during walking vs running in n=26 runners — a
genuine, walking-relevant partial lead, but not a fascicle-length measurement, and its abstract gives
no absolute value (full text closed access).

**A real fallback anchor: Mao, Ren, Huang, Wu & Ruan 2024** (PMID 39228786, *J Sports Sci Med*,
open access, PMCID PMC11366839 — PMID/title/journal/authors independently confirmed live via NCBI
esummary AND full efetch abstract BY THIS SESSION DIRECTLY; the two headline numbers below were ALSO
independently re-fetched and grep-verified by this session directly from the open-access PMC full
text, not merely trusted from the sub-agent that first surfaced this paper). "Fascicle Behavior and
Muscle Activity of The Biceps Femoris Long Head during Running at Increasing Speeds." n=7 college
male sprinters, treadmill running 4/5/6 m/s, real-time B-mode ultrasound of BFlh (middle region) +
EMG + 3D kinematics + 2 force plates. Not walking — but a REAL in-vivo dynamic B-mode ultrasound
measurement of the exact target muscle, independently verified word-for-word from the full text:

> "Peak BFlh fascicle lengths during running are 6.13 ± 0.85, 5.94 ± 0.81, and 5.73 ± 0.82cm at 4,
> 5, and 6m/s, respectively." … "BFlh fascicle lengthening accounted for **about 30%** of MTU length
> change during the late swing phase."

This "~30% of MTU length change" figure is a REAL, in-vivo-measured decoupling ratio for the
analogous eccentric/lengthening late-swing phase — directly comparable in construction to this
session's own computed metric (§6.2). Used as a decorrelated, cross-activity (not cross-muscle)
external check in §6.2.

**Resolving the task's own ~90–110mm figure (a category-error diagnosis, not force-fit):** it does
not match the Kellis in-vivo-ultrasound-validated cadaver numbers (63.5–71.2mm, elderly cadavers) or
Mao 2024's real in-vivo running peak (57.3–61.3mm) — but matches Ward & Lieber 2009's cadaveric
NORMALIZED fiber length almost exactly (97.6mm long head / 110.3mm short head, §4.2 — independently
cross-derived twice, by this session's own manual table transcription AND the sub-agent's
independent raw-HTML-table parse, landing on identical values both times). The most defensible
explanation: Ward & Lieber 2009 is the dataset underlying the generic OpenSim/Hill-type
`optimal_fiber_length` MODEL PARAMETER for BF across the Hamner/Rajagopal/Lai-Arnold model lineage
this session's own model belongs to (confirmed live, Step 0) — so the task's ~90–110mm figure most
plausibly reflects a STATIC/GENERIC MODEL-PARAMETER convention, not an in-vivo WALKING measurement.
Structurally the same kind of category error as this cert family's own precedent (the gastroc doc's
Farris/Sawicki PMID correction) — here a wrong TYPE of number rather than a wrong PMID, caught and
disclosed rather than silently used as if it were walking-gait data. **This causal explanation is
this session's own plausible reconstruction, not independently confirmed against whatever original
source produced the task's figure** — reported as the most defensible available explanation, not an
independently proven fact.

## 5. Forced-adversary #2 — does the reconstruction match SO's own force? (better here than gastroc)

| muscle | n material frames (>5% Fmax) | RMS relative error | max relative error |
|---|---:|---:|---:|
| vaslat_r | 35 | 13.63% | 20.82% |
| bflh_r | 21 | 32.49% | 44.65% |
| bfsh_r | 61 | 21.06% | 51.74% |
| *gasmed_r (gastroc doc)* | *42* | *117.6%* | *362.6%* |
| *gaslat_r (gastroc doc)* | *55* | *72.5%* | *242.4%* |
| *soleus_r (gastroc doc)* | *68* | *63.0%* | *205.0%* |

**Symmetric-QC point:** the knee muscles' compliant reconstruction matches SO's own reported force
MATERIALLY better than the ankle plantarflexors did — a real, disclosed counter-example to a
falsely-uniform "the pattern only gets worse" narrative. The gastroc doc's own mechanistic diagnosis
(quasi-static equilibrium omits the force-velocity multiplier at high fiber-shortening-velocity
frames, concentrated in gastroc's specific high-velocity push-off phase) plausibly explains part of
why: this cert's key windows (VL loading-response, BF terminal-swing) are comparatively lower fiber
velocity than gastroc's push-off, not independently re-verified here (an honest, named, not
re-litigated point — out of scope to re-derive).

## 6. Results in detail

### 6.1 vaslat_r: the loading-response decoupling test — Chleboun's qualitative claim does NOT hold up

Chleboun 2007 reports VL fascicle length statistically UNCHANGED (P=.22, i.e., not distinguishable
from zero) during initial stance while the knee flexes under load. This model shows a REAL, material
fiber-length change instead: fiber_pp=10.79mm against an MTU_pp of 13.89mm in the identical
activation-confirmed window (0–23%GC), a decoupling ratio of **0.777** — far from "near isometric,"
and failing this cert family's own ≤0.50 pass bar (though not crossing the ≥0.80 hard-fail bar,
sitting 0.023 below it). The rigid-tendon bookend (0.946) proves this IS a real tendon-compliance
effect (not a construction tautology) — just a much weaker one than Chleboun's subjects showed. This
model's own SO solution additionally does NOT reproduce Chleboun's second finding (VL reactivating
at terminal swing, fascicle shortening as the knee extends) — `vaslat_r` stays at the activation
floor across 70–100%GC in this trial, with `recfem_r` (rectus femoris) carrying that functional role
instead — a real, disclosed, redundant-actuator difference from Chleboun's real subjects, not the
same muscle-recruitment strategy, though a mechanically valid alternative given real muscle
redundancy.

### 6.2 bflh_r: the terminal-swing timing test PASSES; the decoupling test fails harder than anywhere else in this cert family — now externally corroborated, not just model-internal

The biarticular-hamstring-stretch-at-terminal-swing mechanism (well-established in the sprint-injury
literature, and independently confirmed for THIS trial via pure kinematics before any equilibrium
solve, §2.3) reproduces correctly: global fiber-length maximum at **96.8%GC**, robust to 5-point
smoothing (shift +0.8%GC, still inside the pre-registered wrapped window). But the decoupling ratio
in that same window is **0.931** (compliant) vs 0.981 (rigid) — only a ~5-percentage-point
tendon-compliance effect, the SMALLEST of any muscle examined across both docs, and it clearly
**crosses the ≥0.80 hard "no decoupling" falsifier** — matching soleus_r's own hard-fail in the
gastroc doc (0.868). Absolute tendon excursion is correspondingly tiny: bflh_r's tendon length range
across the full cycle is only 406.82–409.98mm (3.16mm total, Step-1 timeseries), compared to
gastroc's own 8.5–12.5mm individual-tendon excursions.

**This is no longer only a model-internal threshold comparison.** Mao et al. 2024's real in-vivo
B-mode ultrasound measurement of BFlh during running (§4.3, independently re-verified word-for-word
by this session from the open-access full text) reports fascicle lengthening accounts for **"about
30%"** of MTU length change during the analogous late-swing eccentric phase — a real decoupling
ratio of ~0.30. This model's own computed 0.931 for the same construction (fiber_pp/MTU_pp in the
biarticular-stretch window) is roughly **3× higher** — a large, externally-anchored, directionally
unambiguous gap, even granting the honest cross-activity mismatch (running vs walking, and Mao's
probe measured the muscle's MIDDLE region specifically, which this model's lumped-fiber
representation does not distinguish). Both the real in-vivo running data (Mao) and the real in-vivo
walking data (Chleboun, for VL) point the SAME direction this model diverges from: real muscle
fascicles change LESS than the MTU during its own eccentric/lengthening phase; this model's
fascicles track the MTU almost 1:1 regardless of which muscle or which real dataset is used as the
comparison.

### 6.3 Absolute fiber length: vaslat_r runs clearly long; bflh_r runs closer to plausible

`vaslat_r`'s stance-phase range (101.46–141.03mm) sits well above the task's own stated ~80–100mm
envelope, consistent with (and now root-caused by, §4.2) its own `optimal_fiber_length` (141.43mm)
running +2.39 SD above Ward & Lieber's cadaveric mean. `bflh_r`'s full-cycle range (80.40–127.59mm)
straddles the task's stated ~90–110mm envelope more evenly (its own OFL running a milder +0.82 SD
over Ward & Lieber) — a materially SMALLER absolute-scale mismatch than either vaslat_r or gastroc/
gaslat_r show, directly tracking the Ward & Lieber root-cause severity ordering.

### 6.4 bfsh_r: bonus context, negligible tendon buffering, no literature falsifier forced

Full-cycle fiber range 123.51–150.57mm, decoupling 0.952 (compliant) vs 0.967 (rigid) — essentially
no measurable tendon-compliance effect, consistent with its geometrically short tendon
(127.77–128.97mm across the whole cycle, tendon_slack_length 127.77mm) and mono-articular role. No
literature anchor isolating the short head during gait was found this session; reported descriptively
only, per LEAN discipline.

### 6.5 Bonus cross-cycle consistency check

vaslat_r borderline FAILS (max diff 5.04mm vs the pre-registered <5.0mm threshold, only ~1% over) —
traced directly (not just observed) to the pairing between t=0.25s (tail of the prior, only
partially-recorded stride) and t=1.57s (the LAST recorded frame of the whole trial) — a genuine,
small, real stride-to-stride difference plus a plausible recording-edge effect, confirmed by
inspecting the actual frame pair (§2.5), not a pipeline artifact. bflh_r (3.87mm) and bfsh_r (2.60mm)
pass comfortably.

### 6.6 Connection to `docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md` (a plausible, not over-claimed, link)

That doc found the knee/hip force cert is fragile to tendon-slack-length (l_TS) sweeps, mechanistically
driven THERE by a specific antagonist pair, `gasmed_r`+`recfem_r`, discretely switching dominance — a
different specific pair than this doc's own target muscles. The GENERAL mechanism gastroc's own §6.6
identified (a muscle whose tendon absorbs LESS of the total length change than real anatomy is
mechanically MORE exposed to any l_TS mis-calibration) plausibly extends to `vaslat_r`/`bflh_r` too,
given this doc's own measured weak decoupling — but that doc's own root-cause muscle pair was
`gasmed_r`+`recfem_r` specifically, not vaslat_r/bflh_r, so this is offered as a plausible
generalization, not a re-verified finding of that specific document.

## 7. Symmetric-QC

**The adversary this result (leaning toward "the model is broken across the board") must survive:**
is the pattern actually uniform, or is this session cherry-picking the muscles/windows that fail?
Forced by reporting EVERY dimension measured, including favorable ones: bflh_r's terminal-swing
TIMING mechanism passes cleanly and robustly; the compliant-vs-SO-force match is BETTER for all 3
new muscles than for any gastroc/soleus muscle (§5); vaslat_r's decoupling ratio (0.777), while
failing the pass bar, does NOT cross the hard-fail bar (unlike bflh_r/soleus_r); bfsh_r/bflh_r's
absolute-length OFL z-scores (+0.82 to +1.11 SD) are materially milder than vaslat_r/gastrocs'
(+2.34 to +2.92 SD). The picture is genuinely mixed, muscle-dependent, and reported as such — not
flattened into either "everything fails" or "everything passes."

**The adversary this result (leaning toward "the Ward & Lieber comparison is unfair/cherry-picked")
must survive:** is comparing this model's parameters against a modern cadaveric dataset a fair test,
given the model may have been deliberately calibrated to an older convention? Forced by reading the
paper's own discussion (§4.2): Ward & Lieber's own explicit finding is that their (newer) data runs
LONGER, not shorter, than the older Wickiewicz-type convention for exactly these muscle groups — so
this comparison uses the HIGHER, more model-favorable, of the two literature conventions, and this
model's parameters still exceed it by up to +2.92 SD. The comparison is conservative, not
adversarially unfair to the model.

**The adversary this result (leaning toward "the BF anchor gap invalidates the whole BF analysis")
must survive:** does the absence of a clean walking-specific BF ultrasound paper undermine the
bflh_r timing/decoupling findings? No — the TIMING finding (§6.2) rests on independently-confirmed
pure kinematics (this trial's own MTU-length extremum) plus textbook biarticular-muscle mechanics,
not on the missing mm-matched study; only the ABSOLUTE-LENGTH comparison (§6.3) relies on the
fallback Ward & Lieber anchor, disclosed explicitly as such (§4.3), not silently substituted.

## 8. Honest gaps

1. **No in-vivo (any modality) biceps-femoris-during-WALKING measurement was found**, despite an
   exhaustive search: this session's own ~10 direct NCBI query variants, an independently-dispatched
   `watertight-researcher` sub-agent's ~20 further variants + author-targeted sweeps (both negative),
   cross-confirmed by an unrelated third party (Luis et al. 2022, PMID 36277379, hit the identical
   gap for a different purpose and used the same Ward & Lieber 2009 fallback) — this is disclosed as
   a genuine, field-wide data gap (§4.3), not a search failure. The task's own "Kellis" + ~90–110mm
   attribution most plausibly reflects a category error (a static OpenSim model-parameter convention,
   not an in-vivo walking measurement, §4.3) — this session's own reconstruction of why, not an
   independently confirmed fact about the figure's original source. Two real fallback anchors are
   used instead, each disclosed for exactly what it is: Ward & Lieber 2009 cadaveric architecture
   (97.6±26.2mm BFlh, absolute-length comparison, §6.3) and Mao et al. 2024's real in-vivo BFlh
   ultrasound during RUNNING (57.3–61.3mm peak, ~30% decoupling ratio, §6.2) — both independently
   re-verified word-for-word by this session directly from open-access full text, not merely trusted
   from the sub-agent that first located them.
2. **Chleboun 2007's full text (mm-level VL fascicle-length figures) was not accessible** — confirmed
   paywalled via a live Unpaywall API check this session (`is_oa:false`, no repository copy); only
   the abstract's qualitative/P-value claims were used, matching what this session could actually
   verify (no unverified full-text numbers were used, unlike some prior anchors in this cert family
   that relied on a sub-agent's paywalled-source extraction).
3. **vaslat_r's SO-solution does not reproduce Chleboun's second (terminal-swing VL reactivation)
   finding** (§6.1) — a real, disclosed, redundant-actuator-strategy difference (rectus femoris
   carries that role instead in this model's solve), not evidence the model's fiber-length KINEMATICS
   are wrong at that phase (vaslat_r is barely active there, so its fiber trajectory in that window is
   not a meaningfully activation-driven test).
4. **The Ward & Lieber comparison uses summary statistics (mean±SD, n=18–20), not raw per-specimen
   data** — "SD above the mean" is reported as a plausibility/percentile framing, not a formal
   hypothesis-test p-value, since the true shape of each muscle's fiber-length distribution across
   specimens is unknown.
5. **Forced-adversary #2 (§5): the compliant reconstruction does not closely reproduce SO's own
   force** (RMS 13.6–32.5%, max 20.8–51.7%) — better than gastroc/soleus but still a real,
   undiagnosed-in-detail-here gap (the gastroc doc's own fiber-velocity mechanistic explanation is
   referenced, not re-derived for these 3 muscles specifically).
6. **bfsh_r has no literature falsifier** — reported descriptively only; no in-vivo study isolating
   the short head from the long head during gait was found or expected to exist (both heads
   typically image together on standard postero-lateral B-mode ultrasound).
7. **Single subject (subject2), single trial (walking1), right leg only** — same scope caveat as
   every other cert in this family (`docs/MECHANISM_CROSS_SUBJECT.md`).
8. **The §6.6 connection to `docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md` is a plausible generalization,
   not a re-verified finding** — that doc's own root-cause antagonist pair was `gasmed_r`+`recfem_r`,
   not this doc's target muscles.
9. **The Mao et al. 2024 running-vs-walking and probe-region (middle-BFlh-only) mismatches** mean
   §6.2's ~3× decoupling-ratio gap, while large and directionally unambiguous, is not a strictly
   apples-to-apples comparison — no claim is made that a walking-specific measurement would show
   exactly 0.30, only that real in-vivo BFlh behavior (in the one dynamic dataset available) is
   materially more tendon-buffered than this model's reconstruction, in both of the two real datasets
   checked (Mao running; Chleboun walking, for the analogous VL question).
10. **The category-error explanation for the task's ~90–110mm figure (§4.3) is a plausible
    reconstruction, not an independently verified fact** — this session does not have visibility into
    what actually produced that starting figure; the Ward & Lieber 2009 match is machine-verified,
    the causal attribution is this session's own best inference.

## 9. Files

- `scripts/msk/quad_ham_fascicle_kinematics.py` (new) — the full, re-runnable, self-contained
  pipeline, closely modeled on `scripts/msk/gastroc_fascicle_kinematics.py` (same OpenSim API calls,
  same wrap-don't-edit reuse of `validate_joint_force.py`), extended to `vaslat_r`/`bflh_r`/`bfsh_r`
  with muscle-appropriate (not copied) phase windows, full-gait-cycle scope, the geometric
  `ActiveForceLengthCurve`-derived degeneracy gate, the Ward & Lieber 2009 cross-check (Step 8/10),
  and the edge-safe smoothing fix (`smooth_edge_safe`, §2.4). Run with
  `.venv-msk/bin/python3 scripts/msk/quad_ham_fascicle_kinematics.py`.
- `data/msk_smoketest/subject2_walking1/quad_ham_fascicle_kinematics/quad_ham_fascicle_kinematics_results.json`
  (new) — every number in this document, machine-written: muscle properties (incl. AFLC curve
  support), per-frame time series (activation, MTU/fiber/tendon length, pennation, compliant AND
  rigid reconstructions), the geometric identity gate, per-muscle derived statistics (full-cycle +
  muscle-specific windows), the cross-cycle consistency bonus check, the forced-adversary-#2
  force-match diagnostic, the smoothing-robustness check, the Ward & Lieber 2009 architecture
  cross-check (all 6 muscles), the full literature-anchor block, the pre-registered-threshold
  verdict computation, and honest gaps.
- Reused, unedited, zero new OpenSim solve: `data/msk_smoketest/subject2_walking1/
  static_optimization/so/walking1_StaticOptimization_{force,activation}.sto`;
  `OpenSimData/Mocap/IK/walking1.mot` (read-only external drive);
  `OpenSimData/Mocap/Model/LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` (read-only
  external drive); `data/msk_smoketest/subject2_walking1/gastroc_fascicle_kinematics/
  gastroc_fascicle_kinematics_results.json` (sibling cert, read for the gasmed_r/gaslat_r/soleus_r
  context rows in §4.2's Ward & Lieber table only — its own OpenSim solve not touched or re-run).
- External literature verified live via direct `curl` to `eutils.ncbi.nlm.nih.gov` (esearch/
  esummary/efetch) by this session: PMID 17710906 (Chleboun et al. 2007, VL walking, primary),
  PMID 18827324 (Chleboun et al. 2008, VL stairs, companion, not primary), PMID 18972175 (Ward,
  Eng, Smallwood & Lieber 2009, cadaveric architecture cross-check, full text fetched via PMC —
  PMCID PMC2650051), PMID 19646698 / 20727788 / 22564790 (Kellis cadaveric hamstring architecture,
  checked and found NOT to be the walking-gait anchor), PMID 36064431 (Kellis & Blazevich 2022,
  confirmed non-primary review), PMID 39228786 (Mao et al. 2024, BFlh running ultrasound —
  independently confirmed live AND its two headline numbers independently re-fetched/grep-verified
  by this session directly from the open-access PMC full text, PMCID PMC11366839), PMID 36277379
  (Luis et al. 2022, independent third-party corroboration of the BF/ST in-vivo data gap), PMID
  20022251 (Riley et al. 2010, genuine partial lead, MTU not fascicle length), PMID 23540723 and
  PMID 34509042 (near-misses, checked and killed, §4.3). Chleboun 2007's open-access status was
  checked directly by this session via the Unpaywall API (`is_oa:false`, confirmed); a
  `watertight-researcher` sub-agent separately, independently checked BOTH Chleboun papers via 4
  routes (PMC elink, direct jospt.org fetch -> HTTP 403/Cloudflare, Semantic Scholar API, Unpaywall
  API), all agreeing closed -- this session directly re-confirmed only the Unpaywall leg itself, not
  all 4 of the sub-agent's routes.
- A `watertight-researcher` sub-agent was also dispatched this session to independently search for
  the VL/BF literature anchors in parallel with this session's own direct search; its findings
  (the Mao 2024, Riley 2010, and Luis 2022 papers, the near-miss kills, and the Chleboun-paywall
  4-route check) are incorporated above. The load-bearing NEW numeric claims the sub-agent
  introduced were independently re-verified by this session directly before being trusted (per this
  repo's "subagent output is a hypothesis, not truth" discipline): PMID 39228786's existence
  (esummary) AND its two headline numbers (efetch abstract + a direct PMC full-text fetch + grep,
  matching the sub-agent's extraction exactly), and PMIDs 36277379 / 20022251's existence
  (esummary). This session did NOT independently re-derive the sub-agent's own Ward & Lieber 2009
  raw-HTML-table parse (§4.2's table was transcribed by this session directly from its OWN separate
  PMC fetch, before the sub-agent's report arrived, and the two independently matched).
- **Not touched**: `scripts/msk/gastroc_fascicle_kinematics.py`, `scripts/msk/validate_joint_force.py`,
  `scripts/msk/contact_waveform_analysis.py`, `scripts/msk/tendon_elastic_energy.py`, and every
  other pre-existing script in this cert family (all reused unchanged via import or direct JSON
  read, per this repo's wrap-don't-edit convention).

No git commit, no git push performed (isolation respected). All new files are untracked, for the
coordinator to commit.
