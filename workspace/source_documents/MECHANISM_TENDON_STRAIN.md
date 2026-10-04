# MECHANISM ACHILLES + PATELLAR TENDON STRAIN — a decorrelated ultrasound-anchored view of the elastic tendon (2026-07-21)

Computes Achilles (triceps-surae) and patellar (quadriceps) tendon **strain (%)** across the
subject2/`walking1` gait cycle from the model's own Static-Optimization tendon force and its
Millard tendon force-strain curve, anchored against DECORRELATED in-vivo ultrasound literature —
a new, independent observable (strain, not force or fascicle length) probing the same elastic
tendon that `docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md` (l_TS-fragility) and
`docs/MECHANISM_GASTROC_FASCICLE.md` (weak fiber/tendon decoupling) already implicated.

**CONFIDENCE TIER: in-vivo-anchored (ultrasound/ex-vivo mechanical-testing literature, all 5
citations independently verified live this session), single-subject (subject2), single-trial
(`walking1`), right leg only.** Achilles anchor is activity-matched (walking-vs-walking, medium
confidence on the exact 4.6% figure, Sec. 8.2); patellar anchor is force-level-matched but NOT
activity-matched (MVC, not walking — no walking-specific patellar ultrasound anchor exists in the
literature, Sec. 8.5/12.7).

**Headline: the model's tendon strain is REAL, MEASURABLE, and physiologically bounded (nowhere
near the ~8% rupture-risk zone) — but it is MATERIALLY BELOW, not above, the ~4-6% in-vivo walking
band: the falsifier resolves to the "implausible" branch, in the LOW direction the task's own
framing did not name.** Achilles-equivalent peak strain (highest of gasmed_r/gaslat_r/soleus_r,
the fair conservative envelope) is **2.501%** at t=0.51s (44.8%GC, the same push-off instant every
other force cert in this repo already flags), reproduced independently by a genuinely different
OpenSim computation path (2.37–2.79%, Pathway 2 below) and confirmed non-boundary-contaminated.
Patellar-equivalent (quadriceps-lumped-tendon convention) peak is **2.084%** at 43.3%GC
(interior-corrected; the naive whole-trial value hit the recording boundary, a now-familiar
artifact class in this repo, and was corrected, not silently kept). **Both fall clearly below the
task's own 4-6% walking band and even below this session's own generous ±1-percentage-point outer
tolerance (3-7%)** — and both sit at roughly a THIRD of the ~8% rupture-risk/plastic
threshold (ratio 0.25–0.33 vs the 7.5–9.9% Wren 2001 failure-strain anchor, Sec. 11), so
injury-risk is cleanly ruled OUT. A geometric decomposition (Sec. 7) shows this is
overwhelmingly a **force-allocation** gap, not a **curve-shape** gap: the model's generic tendon
curve is broadly realistic when compared at a MATCHED relative force to the literature (Sec. 7.1),
but Static Optimization's own redundant recruitment never asks these muscles for anywhere near the
~92% of Fmax this same curve would require to reach 4.6% strain — the peak relative force SO
actually delivers is only 25-40% of that requirement (Sec. 7.2). This connects concretely, but only
PARTIALLY and non-uniformly (mixed real-EMG evidence, Sec. 8), to the already-found weak tendon/
fiber decoupling and l_TS fragility.

## Headline numbers

| quantity | value | timing | pathway | confidence |
|---|---:|---|---|---|
| **Achilles-equivalent peak strain** (max of gasmed_r/gaslat_r/soleus_r) | **2.501%** | t=0.51s, 44.8%GC | Pathway 1 (primary, force-curve inversion) | high (bit-reproduces prior cert) |
| — cross-check, compliant-equilibrium reconstruction | 2.37–2.79% (per muscle) | t=0.53–0.67s | Pathway 2 (decorrelated) | high (independent method, same conclusion) |
| **Patellar-equivalent peak strain** (max of recfem_r/vasmed_r/vaslat_r/vasint_r; quadriceps-lumped tendon) | **2.084%** (interior-corrected) | t=0.49s, 43.3%GC | Pathway 1 (primary) | high |
| — patellar tendon, anatomically isolated (Blankevoort PT1-6 ligament) | 4.50–6.63% (settled) | varies, 13.9–61.4%GC | Pathway 3 (secondary, heavily caveated) | LOW (known flat-activation defect in source data) |
| Walking anchor band (task-given, Lichtwark & Wilson) | 4-6% | — | literature | medium (see Sec. 8 citation correction) |
| Rupture-risk/plastic threshold (task-given) | ~8% | — | literature, independently verified | high (Sec. 8.4, two convergent sources) |
| F/Fmax this model's curve needs for 4.6% strain | **91.6%** | — | geometric (Sec. 7.2) | high (exact, curve-derived) |
| SO's actual peak F/Fmax, these 7 muscles | 4.2–36.6% | — | Pathway 1 | high |

All numbers machine-computed by `scripts/msk/tendon_strain_gait_cycle.py`, written to
`data/msk_smoketest/subject2_walking1/tendon_strain_gait_cycle/tendon_strain_gait_cycle_evidence.json`
— none transcribed from console prose.

---

## 1. Pre-registration (thresholds fixed in the script BEFORE any strain number was computed)

**C / ¬C, thresholds fixed in advance** (reusing the already-verified-live anchors from
`docs/MECHANISM_GASTROC_FASCICLE.md`/`docs/MECHANISM_TENDON_ELASTIC.md` to set the bands, then
independently re-verifying those specific citations live THIS session, Sec. 8 — not circular: the
threshold values pre-date this session's fresh verification, the fresh verification confirms/
corrects them):

- **PASS (physiological):** peak strain ∈ [3%, 7%] (task's 4-6% band + 1 percentage-point outer
  tolerance for anchor uncertainty) AND clearly below 8%.
- **FAIL_HIGH (implausible/injury-risk):** peak strain > 8%.
- **FAIL_LOW (materially under-strained):** peak strain < 3%.
- **Forced adversary (the one this task is tempted to skip):** a clean "under 4-6%" reading could
  be a construction artifact (wrong muscles, wrong instant, a boundary-search bug, or a
  single-pathway quirk) rather than a real finding. Forced three ways: (a) re-derive the FULL
  per-frame strain trajectory independently from the SAME already-validated
  `tendon_elastic_energy.py` functions and prove bit-exact reproduction of its already-published
  peak_strain scalars (Sec. 2, gate PASS, relerr 0.00e+00 all 7 muscles) before trusting anything
  new; (b) cross-check via a SECOND, genuinely different OpenSim computation route
  (`equilibrateMuscles()`-based reconstruction, Pathway 2) that a prior session already showed
  disagrees with Pathway 1 by 63-118% RMS relative error IN FORCE — so agreement/disagreement in
  STRAIN is a real test, not assumed (Sec. 5); (c) explicitly check every muscle's whole-trial peak
  for the SAME recording-boundary artifact class this repo has repeatedly found elsewhere
  (`docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md` §4, `docs/MECHANISM_TENDON_ELASTIC.md` §5) — caught
  and corrected for recfem_r (Sec. 3).

## 2. Method — reuse, no re-solve (three decorrelated pathways)

**Pathway 1 (PRIMARY, all 7 muscles: gasmed_r/gaslat_r/soleus_r + recfem_r/vasmed_r/vaslat_r/
vasint_r).** Imports `scripts/msk/tendon_elastic_energy.py`'s own already-validated functions
(`get_muscle_props`, `build_inversion_table`, `invert_force_to_strain`) UNCHANGED and applies them
to the full 158-frame SO `force.sto` (already computed, zero re-solve) to recover the FULL
per-frame strain trajectory — that prior script's own JSON stored only the scalar peak, not the
trajectory, so this is a genuine extension, not a duplicate. **Self-audit, run before trusting
anything downstream:** re-derived peak_strain for all 7 muscles bit-exactly reproduces
`tendon_elastic_energy_results.json`'s already-published values (relerr 0.00e+00 for all 7 —
Sec. "reproduction_check" in the evidence JSON) — proves this reuse is faithful, not a subtly
different re-implementation.

**Pathway 2 (cross-check, Achilles-only: gasmed_r/gaslat_r/soleus_r).** Reads the ALREADY-COMPUTED
`gastroc_fascicle_kinematics_results.json` per-frame `tendon_length_compliant_m` (from OpenSim's own
`equilibrateMuscles()`, driven by SO **activation**, a genuinely different computational route than
Pathway 1's force-curve inversion) and derives `strain = L_tendon/L_slack − 1`. Time grids verified
`np.allclose` before comparison (both are the identical 158-row, 0.01s IK/SO grid).

**Pathway 3 (patellar tendon only, anatomically precise): JAM/lenhart2015 native ligament strain.**
`docs/MECHANISM_PATELLOFEMORAL_FORCE.md` already established that the LaiArnold model used by
Pathways 1-2 has **zero** separate patellar-tendon element — the quadriceps' "tendon" in Pathway 1
is the WHOLE lumped quad-tendon+patellar-tendon segment of one continuous `PathActuator` (verified
there: 0 Ligament/PathSpring elements in the model's 93-force `ForceSet`). The JAM/lenhart2015 model
(`data/msk_smoketest/jam_contact_decorr/results/joint-mechanics/gait_driven.h5`, already built for
that same document) DOES have an explicit `Blankevoort1991Ligament` "PT1".."PT6" bundle group
(patella→tibia only) with a **native per-frame `strain` output** already computed by OpenSim's own
`ForsimTool` — read directly here, zero derivation, zero re-solve. Confirmed from the live XML
(`scratchpad/knee_lig/lenhart2015.osim`, PT1: `slack_length=0.046792m`, `linear_stiffness=3000`
N/strain, no custom `reference_strain` override) that "strain" is the anatomically standard
`(length−slack_length)/slack_length` — directly comparable to Pathways 1-2's definition.

## 3. Forced-adversary catch #1 — the recfem_r/patellar boundary artifact (OODA, not silently fixed)

**Observed:** recfem_r's whole-trial peak strain (2.587%) lands at **t=1.57s — the literal last
recorded frame.** **Oriented:** this is the EXACT signature `docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md`
§4 and `docs/MECHANISM_TENDON_ELASTIC.md` §5 already diagnosed for this SAME muscle in this SAME
trial (a genuine, still-rising recruitment burst truncated by the recording boundary, not a real
interior local maximum — `docs/MECHANISM_TENDON_ELASTIC.md` §5 already showed recfem_r's force is
"still rising steeply" at that instant with no sign of turning over). **Decided + Acted:** applied
this repo's own established `EDGE_MARGIN_FRAMES=5` interior-only convention (`vjf.SAVGOL_WINDOW_
PRIMARY//2`, the SAME constant used throughout this cert family) — the corrected, non-boundary
peak is **2.084% at t=0.49s (43.3%GC)**, comfortably interior. The other 6 muscles' whole-trial
peaks are all comfortably interior already (10.9–52.3%GC, nowhere near either edge) — this artifact
affects only recfem_r/the patellar-equivalent envelope it drives, not the Achilles group.

## 4. Forced-adversary catch #2 — the JAM h5's known startup transient (Pathway 3)

A first pass on Pathway 3 found the naive whole-trial peak for ALL 6 PT bundles landing at **t=0.00s
(frame 0)**, with strongly NEGATIVE mean/min strain elsewhere (mean −2.1 to −4.8%, min −38 to −50%)
— implausible on its face. **Oriented:** this matches an ALREADY-DOCUMENTED "violent disengage/
re-engage transient" in this EXACT h5 file (`docs/MECHANISM_PATELLOFEMORAL_FORCE.md` §3: contact
area/force ringing non-monotonically before settling, fixed there via a force-based
`detect_settle_index` giving `sidx=9`, `t=0.09s`). **Decided + Acted:** reused that ALREADY-
COMPUTED, already-published `settle_index=9` verbatim (read directly from
`patellofemoral_force_results.json`, not re-derived) to exclude the first 9 frames. Settled
(t≥0.09s) peaks: PT1 6.63%, PT2 4.89%, PT3 4.50%, PT4 4.64%, PT5 4.89%, PT6 5.78% — now falling
INSIDE or just above the 4-6% walking band. **But** mean/min strain remain strongly negative even
after this correction (these bundles spend most of the trial well below their own slack length) —
consistent with, and very likely explained by, the ALREADY-DISCLOSED defect in this same h5 (§6 of
`docs/MECHANISM_PATELLOFEMORAL_FORCE.md`: quadriceps muscle activation is frozen at a flat 0.05
default for the entire trial, never consuming real time-varying control input) — a patella held at
a constant, low, non-physiological drive level is not expected to track a realistic tension
trajectory. **This pathway is reported as a secondary, heavily-caveated cross-check, not the primary
patellar-tendon number** — its brief-window peaks landing near/above the walking band, while its
dominant (most-of-trial) behavior is clearly non-physiological, is exactly why it is not promoted
above Pathway 1.

## 5. Cross-pathway agreement (Achilles) — the forced adversary FALLS, strengthening the finding

| muscle | Pathway 1 peak (%) | Pathway 2 peak (%) | RMS diff (pp) | max diff (pp) |
|---|---:|---:|---:|---:|
| gasmed_r | 2.501 | 2.790 | 0.416 | 1.426 |
| gaslat_r | 1.983 | 2.367 | 0.324 | 1.033 |
| soleus_r | 1.856 | 2.617 | 0.295 | 1.185 |

Despite `docs/MECHANISM_GASTROC_FASCICLE.md` already finding these two pipelines' IMPLIED **FORCE**
disagrees by 63–118% RMS relative error, their implied **STRAIN** agrees far more closely (RMS
0.3–0.4 percentage points, both pathways landing in the SAME 1.9–2.8% band, both clearly below the
3% floor). This is not a coincidence but a **geometric** consequence of the curve's own shape: the
tendon force-strain curve is concave (soft "toe" region at low force, stiff linear region near
F=Fmax, Sec. 7), so `d(strain)/d(F)` is LARGE at the low force levels these muscles actually operate
at during gait — a given fractional force error maps to a much smaller absolute strain error. The
adversary ("maybe Pathway 1 is a single-method artifact") was forced via a genuinely independent
computation and FELL: both agree the muscle is under-strained relative to the walking anchor.

## 6. Gait-cycle trajectory (10%GC bins, envelope = max across the group's muscles per frame)

Strain is **intensive** (a length ratio), not extensive like force/energy — it must NOT be summed
across parallel tendons the way `docs/MECHANISM_TENDON_ELASTIC.md` correctly sums stored ENERGY or
`docs/MECHANISM_PUSHOFF_PLANTARFLEXOR.md` correctly sums FORCE across gasmed_r+gaslat_r+soleus_r.
The fair, non-cherry-picked summary is the per-frame MAX across the group (the single most-strained
tendon at each instant — using the LOWEST-strain muscle would understate; using the highest is the
conservative, adversarial-to-my-own-under-strain-finding choice):

| %GC bin | Achilles-envelope (%) | Patellar-envelope (%) |
|---|---:|---:|
| [0,10) | 0.26 | 0.52 |
| [10,20) | 0.58 | 1.20 |
| [20,30) | 1.18 | 1.57 |
| [30,40) | 1.56 | 1.66 |
| **[40,50)** | **2.34** | **1.92** |
| [50,60) | 1.65 | 1.74 |
| [60,70) | 0.46 | 1.20 |
| [70,80) | 0.25 | 0.63 |
| [80,90) | 0.26 | 0.18 |
| [90,100) | 0.35 | 0.19 |

Both groups peak in the same 40-60%GC push-off/loading-response window already established as this
model's own force/contact peak region (`docs/MECHANISM_CONTACT_WAVEFORM.md`) — the timing is
physiologically sensible (matches "isometric-ish stance, tendon loads at push-off," the same
qualitative mechanism `docs/MECHANISM_GASTROC_FASCICLE.md` already confirmed) even though the
MAGNITUDE falls short of the literature band.

## 7. Geometric decomposition — curve-SHAPE realism (Q1) vs force-ALLOCATION realism (Q2)

All 11 of this model's muscles share one generic, model-wide default tendon curve (`strain_at_one_
norm_force` e0=4.9% at F/Fmax=1.0, already disclosed in `docs/MECHANISM_TENDON_ELASTIC.md` §7.1 as
the "dominant uncertainty"). Two DIFFERENT, decorrelated questions, kept separate (geometric
derivation, not a heuristic comparison of mismatched activities):

### 7.1 Q1 — is the curve SHAPE itself realistic, at a matched relative force?

| comparison point | this model's curve predicts | literature (verified live, Sec. 8) | ratio |
|---|---:|---|---:|
| F/Fmax=0.50 (≈50% MVC) | 3.086% | Farris 2013: free-tendon 5.2±1.7%, proximal AT 2.6±2.0% (50% MVC, isokinetic dynamometer) | 0.59–1.19× (inside the reported range) |
| F/Fmax=1.00 (≈MVC) | 4.900% (=e0, by construction) | Hansen 2006: patellar tendon 6.8–6.9% (MVC isometric ramp) | 0.71–0.72× |

At a FORCE-MATCHED comparison, this model's generic default curve is **not wildly unrealistic** —
it sits inside or just below the real, directly-measured human range at the SAME relative loading
level. The curve-shape/compliance parameter is a real, disclosed uncertainty (per
`docs/MECHANISM_TENDON_ELASTIC.md` §7.1) but is NOT, by itself, sufficient to explain the size of the
gait under-strain finding.

### 7.2 Q2 — does gait ever ask enough of these muscles to reach the walking anchor?

Inverting the SAME curve the other direction: reaching 4.6% strain (the walking-specific literature
figure) requires **F/Fmax = 91.6%** on this model's curve — for EVERY one of the 7 muscles (shared
generic curve). Static Optimization's own actual peak relative force, this trial:

| muscle | SO peak F/Fmax | ratio to the 91.6% needed |
|---|---:|---:|
| gasmed_r | 36.6% | 0.40 |
| recfem_r | 27.8% | 0.30 |
| gaslat_r | 25.8% | 0.28 |
| soleus_r | 23.4% | 0.26 |
| vaslat_r | 13.5% | 0.15 |
| vasmed_r | 6.7% | 0.07 |
| vasint_r | 4.2% | 0.05 |

**No muscle gets closer than 40% of the way to the relative force this model's own curve would need
to reach the literature's walking strain figure.** This is the dominant, geometrically clean
explanation: the gap is overwhelmingly a FORCE-ALLOCATION story (SO's redundant multi-muscle
recruitment keeps every one of these muscles well under half the force fraction needed), not
primarily a curve-compliance story (Sec. 7.1 shows the curve itself is roughly in the right
ballpark at matched force).

## 8. External literature anchors — verified live (this session, directly + via one sub-agent cross-check)

**A forced adversary the task's own citation invited, caught and corrected:** the task's own
"Achilles peak strain in walking ~4-6% (Lichtwark & Wilson 2005...)" attribution is WRONG on the
paper. Independently re-verified live, twice — once by a `watertight-researcher` sub-agent and once
directly by this session via `curl` to `eutils.ncbi.nlm.nih.gov` (both esummary and full efetch
abstract text, matching this repo's established citation-verification convention):

1. **Lichtwark GA, Wilson AM. "In vivo mechanical properties of the human Achilles tendon during
   one-legged hopping." J Exp Biol. 2005 Dec;208(Pt 24):4715-25. PMID 16326953, DOI
   10.1242/jeb.01950.** This IS "Lichtwark & Wilson 2005" — but it is a **ONE-LEGGED HOPPING**
   study, not walking. Verbatim (fetched directly, this session): *"the average peak strain was
   8.3%... prolonged hopping may well cause tendon damage."* Used here as the verified anchor for
   the task's OWN "~8% rupture-risk" framing (the source paper itself links high strain to damage
   risk), not for the 4-6% walking figure.
2. **Lichtwark GA, Wilson AM. "Interactions between the human gastrocnemius muscle and the Achilles
   tendon during incline, level and decline locomotion." J Exp Biol. 2006 Nov;209(Pt 21):4379-88.
   PMID 17050853, DOI 10.1242/jeb.02434.** The ACTUAL walking(/incline/decline)-specific paper.
   Verbatim abstract (fetched directly): *"Substantial levels of Achilles tendon strain were
   recorded in both walking and running conditions."* The specific "**4.6% (10.9±1.7mm)**" figure
   already used in `docs/MECHANISM_GASTROC_FASCICLE.md` comes from a PRIOR session's full-text
   extraction — this session independently attempted to re-verify it directly (PMC: no deposit;
   Europe PMC: confirmed non-open-access) and **could not access the full text either** — the
   PMID/topic/qualitative claim is now triply-confirmed (2 prior docs + this session, all live), but
   the exact "4.6%" number remains at medium confidence, not independently re-derived from primary
   text this session. Disclosed, not silently upgraded.
3. **Farris DJ, Trewartha G, McGuigan MP, Lichtwark GA. "Differential strain patterns of the human
   Achilles tendon determined in vivo with freehand three-dimensional ultrasound imaging." J Exp
   Biol. 2013 Feb 15;216(Pt 4):594-600. PMID 23125339, DOI 10.1242/jeb.077131.** THE "Farris" paper
   the task's framing pointed to (confirmed: this is a direct Achilles-tendon-STRAIN paper, unlike
   Farris & Sawicki 2012 PNAS which reports MTU/fascicle length in mm but no strain %, already
   checked in `docs/MECHANISM_GASTROC_FASCICLE.md`). Verbatim (fetched directly, this session):
   *"The free-tendon exhibited significantly (P=0.03) greater longitudinal strain (5.2±1.7%) than
   the proximal AT (2.6±2.0%)"* — measured during **50% MVC, fixed-end isokinetic-dynamometer
   plantar-flexor contractions**, NOT walking. Used here for the force-matched Q1 comparison
   (Sec. 7.1), correctly labeled by its real condition, not implied to be a walking measurement.
4. **Wren TA, Yerby SA, Beaupré GS, Carter DR. "Mechanical properties of the human achilles
   tendon." Clin Biomech (Bristol). 2001 Mar;16(3):245-51. PMID 11240060, DOI
   10.1016/s0268-0033(00)00089-9.** Ex-vivo excised-tendon mechanical testing. Verbatim (fetched
   directly): *"The failure strain was 12.8% (SD 1.7) for the bone-tendon complex and 7.5% (SD 1.1)
   for the tendon substance"* at 1%/s; **9.9% (SD 1.9) at the faster, more physiologically relevant
   10%/s rate.** The primary, independently-verified anchor for the task's "~8% rupture-risk/plastic
   threshold" (an ULTIMATE/failure-strain framing specifically, not a separate "yield-onset"
   number — disclosed precisely, not blurred).
5. **Hansen P, Bojsen-Moller J, Aagaard P, Kjaer M, Magnusson SP. "Mechanical properties of the
   human patellar tendon, in vivo." Clin Biomech (Bristol). 2006 Jan;21(1):54-8. PMID 16183183, DOI
   10.1016/j.clinbiomech.2005.07.008.** Verbatim (fetched directly): *"strain [trial a, 6.9 (SEM
   0.6)%; trial b, 6.8 (SEM 0.7)%]"* during **maximal 10-s ramp isometric knee extension** (MVC), not
   walking. No walking-specific in-vivo patellar tendon strain measurement was found to exist in
   the literature (consistent with `docs/MECHANISM_PATELLOFEMORAL_FORCE.md` §8's own already-
   verified finding that no in-vivo instrumented PFJ walking measurement exists either) — used here
   as the best available, correctly-labeled anchor (Sec. 7.1's Q1 check), not a walking figure.

## 9. Connecting to the l_TS + weak-decoupling findings — real, but partial and non-uniform

**Consistent, convergent signal, three independent angles, same underlying mechanism candidate:**
1. `docs/MECHANISM_GASTROC_FASCICLE.md`: this model's Achilles-equivalent tendon absorbs LESS of the
   total stance-phase MTU length change than real anatomy (fiber does 79–87% of the work vs the
   literature's tendon-dominated "catapult" picture; individual-tendon excursion only 8.5–12.5mm vs
   Lichtwark & Wilson's reported 24.9mm combined). A tendon stretching less than real anatomy is,
   geometrically, a tendon carrying less relative force than real anatomy — directly consistent
   with (not duplicative of) this session's strain finding.
2. `docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md`: the knee/hip force cert is fragile to l_TS because
   Static Optimization's redundant recruitment discretely switches an antagonist pair (gasmed_r/
   recfem_r) between near-silent and dominant across a narrow l_TS window — the SAME "SO's
   force-allocation among these exact muscles is unstable/non-physiological" character this
   session's Sec. 7.2 finding also points to (SO never commits enough force to any of the 7 muscles
   individually).
3. `docs/MECHANISM_PUSHOFF_PLANTARFLEXOR.md`: real EMG cross-check gives a MIXED, not uniform,
   picture — **soleus_r's real EMG implies 2.81× MORE force than SO chose** at push-off (strongly
   consistent with an SO-under-recruitment explanation for soleus_r's own under-strain, the LOWEST
   of the three Achilles muscles here, 1.86%) — but **gasmed_r's real EMG implies LESS force than SO
   chose** (ratio 0.805, the OPPOSITE direction), and gaslat_r has no EMG channel for this subject.
   **This session's own dominant Achilles-envelope driver is gasmed_r** (2.501%, the LEAST
   under-strained of the three) — for that specific muscle, the available real-EMG evidence does
   NOT support "SO under-recruited it relative to real physiology"; if anything the reverse. **Even
   applying soleus_r's own 2.81× EMG correction does not close the full gap to the 91.6% F/Fmax
   this curve needs** (2.81×23.4%≈66%, still short of 91.6%) — so the connection to prior findings
   is real and geometrically forced (Sec. 7.2 is an exact fact about this model, independent of any
   EMG evidence), but its FULL attribution to "SO mis-recruits vs real physiology" is only partially
   and non-uniformly evidenced, not a clean single-mechanism resolution. Reported honestly, not
   oversold into a tidier story than the data supports.

## 10. Symmetric QC

**The adversary this result (leaning toward "confirms prior findings," a positive-feeling
narrative) must survive:** is the convergence with the l_TS/decoupling findings just narrative
pattern-matching? Forced via Sec. 9.3's own real-EMG check: NO clean uniform story — gasmed_r's own
EMG evidence points the OPPOSITE direction from soleus_r's, and even soleus_r's own correction
doesn't fully close the gap. The convergence is real at the LEVEL OF THE GEOMETRIC FACT (SO never
approaches 91.6% F/Fmax for any of these 7 muscles — Sec. 7.2, unconditionally true, no EMG
needed) but not at the level of a single, fully-resolving causal story.

**The adversary this result (leaning toward "under-strain finding is real, not an artifact") must
survive:** forced three ways in Sec. 1/3/4/5 (bit-exact reproduction of the prior cert, explicit
boundary-artifact correction for recfem_r, a genuinely independent second computation pathway) —
all three checks left the FAIL_LOW conclusion intact or slightly reinforced it (Pathway 2's numbers
are if anything marginally HIGHER, not lower, than Pathway 1's).

## 11. Machine cross-checks (gates, all PASS unless noted)

| gate | result |
|---|---|
| Re-derived peak_strain bit-reproduces `tendon_elastic_energy_results.json`, all 7 muscles | **PASS**, relerr 0.00e+00 |
| Round-trip force→strain→force inversion, all 7 muscles | **PASS**, max relerr 8.00e-06 |
| Pathway 1/Pathway 2 time-grid alignment (Achilles) | **PASS**, `np.allclose` |
| recfem_r/patellar-envelope boundary-artifact check + correction | **PASS** (caught, corrected, disclosed) |
| Pathway 3 startup-transient check (reused sidx=9 from sibling doc) | **PASS** (caught, corrected, disclosed) |
| Achilles-envelope peak is non-boundary (interior) | **PASS** (frame 51 of 158, well inside [5,153)) |
| All 5 literature PMIDs independently verified live (esummary + full efetch abstract, this session directly) | **PASS**, 1 citation-attribution correction made (Sec. 8.1) |
| Achilles peak strain (2.50%) vs rupture-risk threshold (7.5–9.9%, Wren 2001) | **PASS** (well below, ratio 0.25–0.33) |
| Curve-shape check (Q1) vs Farris 2013 (force-matched) | **PASS** (0.59–1.19× the reported range) |

Full detail, every PMID/DOI/verbatim quote, every per-frame array: `data/msk_smoketest/
subject2_walking1/tendon_strain_gait_cycle/tendon_strain_gait_cycle_evidence.json`.

## 12. Honest gaps

1. **The exact "4.6%" Lichtwark & Wilson 2006 walking figure remains at medium confidence** — the
   PMID/topic/qualitative claim is independently, live, triply-confirmed (2 prior sessions + this
   one), but the specific number could not be re-derived from primary text this session either
   (non-open-access, no PMC deposit, Europe PMC confirmed non-OA) — inherited gap, not newly
   introduced, disclosed rather than silently upgraded to "verified."
2. **Pathway 3 (patellar-tendon-specific, anatomically precise) is low-confidence** — its settled
   peaks (4.5–6.6%) are intriguing (closer to the walking band than Pathway 1's lumped-tendon
   number) but its trial-long mean/min strain is strongly negative, most plausibly an artifact of
   the ALREADY-DISCLOSED flat/constant 0.05 muscle-activation defect in that h5 (not newly found
   here, carried forward from `docs/MECHANISM_PATELLOFEMORAL_FORCE.md` §6). Not treated as primary.
3. **Farris 2013 and Hansen 2006 are NOT walking measurements** (50% MVC and MVC respectively) —
   used deliberately for the FORCE-MATCHED Q1 comparison (Sec. 7.1), correctly labeled throughout,
   never implied to be walking-condition numbers.
4. **The l_TS/decoupling connection (Sec. 9) is geometrically forced but causally partial** — real
   EMG evidence is mixed across the three Achilles muscles (strongly supports under-recruitment for
   soleus_r, contradicts it for gasmed_r); no uniform single-mechanism resolution is claimed.
5. **Single subject (subject2), single trial (`walking1`), right leg only** — same scope caveat as
   every other cert in this family.
6. **The generic, model-wide default tendon curve (e0=4.9%, all 11 muscles identical) is a real,
   disclosed uncertainty** (inherited from `docs/MECHANISM_TENDON_ELASTIC.md` §7.1) — Sec. 7.1 shows
   it is roughly in the right ballpark at matched relative force, but it is still not a
   subject-measured, ultrasound-calibrated compliance value.
7. **No walking-specific in-vivo patellar tendon strain anchor exists in the literature** (checked
   live by the sub-agent and consistent with this repo's own prior PFJ-force finding) — the patellar
   side of this report is necessarily anchored at MVC-condition literature (Hansen 2006), a weaker
   match than the Achilles side's activity-matched anchor (Lichtwark & Wilson 2006).
8. **Pathway 2 exists only for the 3 Achilles muscles**, not the 4 quadriceps muscles (the
   `gastroc_fascicle_kinematics.py` reconstruction was built Achilles-only) — the patellar-
   equivalent Pathway-1 number has no equivalent independent cross-check beyond Pathway 3's
   heavily-caveated one.

## 13. Files

- `scripts/msk/tendon_strain_gait_cycle.py` (new) — the full, re-runnable, self-contained pipeline
  (imports `tendon_elastic_energy.py`, `gastroc_fascicle_kinematics.py`, `validate_joint_force.py`
  unchanged; reads `gait_driven.h5` directly via `h5py`). Run with
  `.venv-msk/bin/python3 scripts/msk/tendon_strain_gait_cycle.py` (needs OpenSim 4.6 + h5py).
- `data/msk_smoketest/subject2_walking1/tendon_strain_gait_cycle/tendon_strain_gait_cycle_evidence.json`
  (new) — every number in this document: full per-frame strain trajectories (all 3 pathways),
  the reproduction self-audit, the boundary-artifact corrections, the geometric Q1/Q2 decomposition,
  the 10%GC-binned trajectory, the verdict computation, and all literature anchors with verbatim
  quotes/PMIDs/DOIs.
- Reused, unedited, directly imported: `scripts/msk/tendon_elastic_energy.py`, `scripts/msk/
  gastroc_fascicle_kinematics.py`, `scripts/msk/validate_joint_force.py`. NOT imported, but its
  already-computed OUTPUT reused verbatim (the `settle_index=9` value, read directly from
  `patellofemoral_force_results.json`, itself produced by `scripts/msk/jam_contact_decorr_analyze.
  py`'s `detect_settle_index`): the settle-window correction in Sec. 4 was reused, not re-derived.
- Reused, unedited data: `data/msk_smoketest/subject2_walking1/static_optimization/so/
  walking1_StaticOptimization_force.sto`; `data/msk_smoketest/subject2_walking1/
  gastroc_fascicle_kinematics/gastroc_fascicle_kinematics_results.json`; `data/msk_smoketest/
  subject2_walking1/tendon_elastic_energy/tendon_elastic_energy_results.json`; `data/msk_smoketest/
  jam_contact_decorr/results/joint-mechanics/gait_driven.h5`; `data/msk_smoketest/
  patellofemoral_force/patellofemoral_force_results.json` (settle_index only); the scaled
  `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` (read-only, external drive).
- External literature verified live via `curl` to `eutils.ncbi.nlm.nih.gov` (esummary + efetch),
  both by a `watertight-researcher` sub-agent and independently re-checked directly by this session:
  PMID 16326953 (Lichtwark & Wilson 2005, hopping), PMID 17050853 (Lichtwark & Wilson 2006,
  walking), PMID 23125339 (Farris et al. 2013), PMID 11240060 (Wren et al. 2001), PMID 16183183
  (Hansen et al. 2006).

No git commit, no git push performed (isolation respected). All new files are untracked.
