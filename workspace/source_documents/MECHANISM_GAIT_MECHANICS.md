# MECHANISM GAIT MECHANICS — the inverted-pendulum energy-exchange model, tested against the twin's OWN measured walk (2026-07-22)

The function baseline the operator's dysfunction-as-deviation reasoning depends on: normal walking is an
**inverted pendulum**, not a spring. The stiff stance leg vaults the body's center of mass (CoM) over the
foot; gravitational potential energy (PE, peak at mid-stance, CoM highest) and kinetic energy (KE, peak at
double-support, CoM fastest) trade **out of phase**, recovering a large share of mechanical energy without
muscular cost — pendulum-style, not spring-mass (running). Script: `scripts/msk/gait_mechanics.py`. Raw
machine-checked results: `reports/probes/gait_mechanics_results.json` (7 pre-registered gates, **7/7 PASS**,
plus descriptive/robustness detail below). Citation ledger: `docs/MECHANISM_GAIT_MECHANICS_evidence.json`.

**Scope / duplicate-check** (`docs/` grepped for pendulum/Cavagna/Froude/double-support/M-shape/spring-mass;
`data/MECHANISM_ANCHOR_GRAPH.json`, 1057 nodes, grepped live for the same terms — not assumed): no existing
doc computes CoM energy-recovery, double-support fraction, GRF M-shape, or a Froude-number check for this
twin's own walking trial. The one graph hit, `MSK-GAIT-LOCOMOTION` (OPEN, SEED-DESIGN only, alias
`AUTO-TEST-THE-CORE-REDUCED-REP-GENERALIZES-CL`), proposes parametrizing gait *control* in the inverted-
pendulum/SLIP model's low-dim latent space to test a **different, unrelated** hypothesis (whether reduced-
rank representations generalize outside CAD/topology-optimization) — it needs a validated forward pendulum
model as an input; this doc supplies exactly that (not that node's own ambition), the same relationship
`docs/MECHANISM_SPRINT_SPRING_MASS.md` has to it on the running side. `docs/MECHANISM_GAIT_KINEMATICS_FIDELITY.md`
already certified (16/17 pre-registered gates PASS) that subject2_walking1's joint kinematics **are** normal,
healthy, textbook-consistent gait — the precondition this doc's "function baseline" framing leans on rather
than re-derives.

## Data: the twin's OWN measured walk, not literature alone

- **Raw force-plate GRF** (2000 Hz, right+left, 1.579 s): `/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/ForceData/walking1_forces.mot`
  — this is the **exact file** every `*_setup_externalLoads_patched.xml` under
  `data/msk_smoketest/subject2_walking1/` already points `<datafile>` at (verified by reading one directly);
  every existing ID/JR/CMC cert on this subject is already built on it. Read-only, never modified.
- **Static standing calibration trial** (2000 Hz, 8 s, same subject/session):
  `.../subject2/ForceData/static1_forces.mot` — an independent, decorrelated mass measurement via F=mg,
  not an assumed anthropometric value.
- **Twin's own IK kinematics** (100 Hz, 1.57 s), inside the read-only smoketest tree:
  `data/msk_smoketest/subject2_walking1/walking1_smoketest.mot`.
- **Twin's own subject-scaled OpenSim model**:
  `data/msk_smoketest/subject2_walking1/rra/pass1/subject2_walking1_rra_pass1_adjusted.osim` — used for a
  direct geometric leg-length measurement (joint-offset XML, not an anthropometric ratio).

Method is geometric, not curve-fit: Newton's second law integrated twice on the **raw, measured** ground
reaction force gives the CoM's vertical position and forward velocity — the classic Cavagna force-plate
method — with kinematics used only as an *independent* cross-check of mean speed and leg length, so two
decorrelated data streams (force-plate dynamics vs. marker/IK kinematics) validate each other.

## Part 1 — Claim (a): anti-phase KE↔PE exchange, forced against the "muscles do it all" adversary

**C** (pre-registered): the twin's own walking1 stance cycle shows KE and PE trading **out of phase**
(r(KE,PE) < 0) with a substantial, Cavagna-consistent %recovery, in **real, unfabricated** half-step windows
bounded by directly-observed gait events (heel-strike/toe-off), at **two independent mass estimates**.
**¬C / adversary, steelmanned**: "no obligatory pendulum-like exchange — muscles independently drive KE and
PE with no privileged phase relationship." Forced to its strongest fair form as a **circular-shift null**
(house rule for this kind of decorrelation control: preserves each signal's own autocorrelation/amplitude/
frequency content, destroys only the cross-phase relationship — harder to beat than a naive full reshuffle,
which would also destroy each signal's own smoothness and so under-sell the adversary).

**Measured** (raw GRF, `Fy_tot=R_vy+L_vy`, `Fx_tot=R_vx+L_vx`; `a=F/m`, double-integrated with two
physically-motivated calibration constants — zero net vertical drift over the unit, mean forward speed
matched to the independently-measured kinematic value):

| window (real gait events, no reconstruction) | mass | recovery % | r(KE,PE) | null percentile rank |
|---|---|---|---|---|
| L-heelstrike→R-heelstrike (0.579–1.2425 s) | 78.2 kg (established) | **80.0%** | **−0.954** | 99.5 |
| L-heelstrike→R-heelstrike | 79.24 kg (static-trial) | **75.2%** | **−0.901** | 95.1 |
| R-toeoff→L-toeoff (0.8025–1.468 s) | 78.2 kg | **52.2%** | **−0.599** | 70.4 |
| R-toeoff→L-toeoff | 79.24 kg | **61.2%** | **−0.746** | 79.3 |

**Pre-registered gates (all machine-checked, all PASS):** r(KE,PE)<0 in 4/4; recovery% inside a
[30,85]% plausibility band in 4/4; real value beats >50% of the circular-shift null in 4/4 (adversary falls
on direction in all four; 2/4 additionally clear a strict ≥95th-percentile bar). The anti-phase sign is
unanimous (4/4, r ranging −0.60 to −0.95) and recovery is always large and positive — never near the
~0% the adversary predicts.

### The OODA loop this doc forces on itself before trusting the spread

First pass showed 52–80% recovery across the two half-step windows — a real spread, not noise, and it
would have been premature to average-and-move-on. **Orient**: a half-step is not the classic Cavagna
periodic unit; it brackets the *same* central single-support vault with **mechanistically different**
double-support sub-phases — DS1 (heel-strike **collision**, net deceleration) vs. DS2 (**push-off**, active
plantarflexor work, the Kuo/Donelan step-to-step-transition distinction). A window ending at DS2 includes
active muscular energy *addition*, which necessarily depresses the *apparent passive* recovery fraction
relative to a window that only brackets DS1. This is a real physiological distinction, not an artifact —
and it directly explains (not hand-waves) why the heelstrike-anchored window (DS1-only) reads higher
(75–80%, right at Cavagna's own literature ceiling) than the toeoff-anchored window (DS2-included, 52–61%).

**Act (the fix, not a shrug):** reconstruct the classic full-stride unit. The file is truncated before the
next L heel-strike, so the fourth quarter-phase (single-R-support) is missing; it is completed by mirroring
the *measured* single-L-support segment — justified by this trial's own measured L-R symmetry (DS1 vs DS2
differ by only **0.89%**, see Part 2). Result: **60.1% / 68.7%** recovery (both masses), r=−0.58/−0.81,
beating the null at the 91.7–92.4th percentile — squarely inside Cavagna's own literature band (below).

**Specificity check on the fix itself** (forced, not assumed): splicing in a **time-reversed** copy of the
same segment — which has no physical justification, since gait's symmetry is a phase-shift + left-right
relabeling, not a time reversal — gives a markedly weaker, sometimes sign-flipping result (r=+0.07/−0.15,
recovery 38.6%/50.4%). The pipeline does not manufacture a pendulum signal from arbitrary input; only the
physically-motivated reconstruction produces the strong, consistent result. This is reported as a disclosed,
real-data-sensitive approximation (honest gap below), not hidden behind the headline number.

## Part 2 — Claim (b): the twin's OWN GRF shows the double-hump (M-shape) + a walking-consistent double-support fraction

**C:** the fully-captured single-limb stance (left foot, 0.579–1.468 s, 0.889 s duration) shows the
classic **two-peak, one-trough** vertical GRF shape, machine-detected (`scipy.signal.find_peaks`,
prominence≥20 N — never eyeballed), and double-support occupies a nonzero, walking-consistent (not
running, not degenerate) fraction of the gait cycle.

**Measured:** peak 1 (loading response) = 775.2 N = **1.011×BW** @ 0.803 s; trough (mid-stance unloading)
= 666.8 N = **0.869×BW**; peak 2 (push-off) = 825.5 N = **1.076×BW** @ 1.200 s. Trough sits 14–19% below
both peaks — **gate PASS** (exactly 2 peaks, 1 trough, trough <95% of the lower peak). Gait timing (all
directly observed from stance/liftoff events, 20 N standard footfall threshold): L heel-strike 0.579 s,
R toe-off 0.8025 s, R heel-strike 1.2425 s, L toe-off 1.468 s → DS1=0.2235 s, DS2=0.2255 s (**0.89%
left-right asymmetry** — good bilateral symmetry, itself a sanity check on data quality), stride estimate
(2× the directly-measured L-heelstrike→R-heelstrike step time) = **1.327 s**, double-support fraction
= 0.449/1.327 = **33.8%** of the cycle. **Gate PASS**: strictly inside the (0%, 50%) band that separates
walking (has double support) from running (0%, flight phase) and from a degenerate/mislabeled trial (>50%).

**Honest, non-hidden tension with the commonly-quoted "~20%" figure:** 33.8% is meaningfully above the
textbook headline. This is not swept past: the twin's own kinematics put this trial's forward speed at
**1.060 m/s** (measured from `pelvis_tx`, below), slower than the ~1.2–1.4 m/s typically associated with
the "~20%" figure, and double-support fraction is well-established to rise as speed falls (more time is
needed in bipedal support to redistribute momentum when less forward momentum is available per stride).
`Öberg, Karsznia & Öberg (1993, PMID 8035350)`, live-verified, confirms the existence of exactly this
speed-stratified normative structure (separate reference tables for slow/normal/fast gait) but its
abstract alone does not hand over the exact double-support-vs-speed slope — flagged in honest gaps rather
than asserted from memory.

**Axis/sign sanity check** (machine-checked, not assumed): mean horizontal GRF over the first quarter of
single-limb stance = **−56.5 N** (braking, correct sign) vs. the last quarter = **+65.1 N** (propulsive,
correct sign) — confirms the coordinate convention before any energy integration is trusted. **Gate PASS.**

## Part 3 — Claim (c): Froude number confirms this is ordinary walking, well clear of the walk-run transition

**C:** Fr = v²/(gL) for this trial sits comfortably below the empirically-established ~0.5 human walk-run
transition value — internal consistency of the "walking" label, derived geometrically from the twin's own
model, not assumed.

**Leg length, derived from the twin's own scaled model geometry** (not an anthropometric ratio): summing
the vertical joint-offset translations along the hip→knee→ankle kinematic chain in
`subject2_walking1_rra_pass1_adjusted.osim` gives thigh=0.4910 m + shank=0.5034 m + ankle-ground clearance
(≈0.07 m, the one coarse term) = **1.0645 m**. Independent cross-check: mean `pelvis_ty` over the gait
cycle (1.1443 m) minus the model's own hip-joint offset within the pelvis body (0.0819 m) = **1.0624 m** —
**0.19% agreement** between two independent geometric derivations from the same model.

**Measured:** v_avg (kinematic, `pelvis_tx` over the full 1.57 s trial) = **1.060 m/s**. **Fr = 1.060² /
(9.81 × 1.0645) = 0.108** — a **4.6× margin** below the ~0.5 transition. **Gate PASS.** Numeric anchor,
live-verified this session: `Carr & McGee (2009, PMID 19672305)`, PLoS ONE, states directly in its
abstract: *"Human unsuited Fr* is relatively constant (approximately 0.5) with gravity."* `Diedrich &
Warren (1995, PMID 7707029)` independently confirms (treadmill experiments) that the human walk-run
transition occurs "at a constant speed near the energy separatrix," the same phenomenon Fr=0.5
parametrizes; `Vaughan & O'Malley (2005, PMID 15760752)` confirms Fr=v²/gL (with L = leg length) as the
standard dynamic-similarity formalism applied to bipedal gait.

## Part 4 — anchor: does the pendulum-recovery claim reproduce the primary literature figure?

`Cavagna, Heglund & Taylor (1977, PMID 411381)`, Am J Physiol, live-verified, states directly: *"Walking...
involves an alternate transfer between gravitational-potential energy and kinetic energy within each stride
(as takes place in a pendulum). This transfer... can account for up to 70% of the total energy changes...
leaving only 30% to be supplied by muscles. No kinetic-gravitational energy transfer takes place during
running."* This doc's twin-measured, full-stride-reconstructed recovery (60.1%/68.7%, both mass hypotheses)
sits inside this exact band; three of the four real half-step measurements (52.2–75.2%) do too, and the
fourth (80.0%) exceeds the cited ceiling by 10 points — reported as-is (Part 1's OODA section explains the
mechanism, DS1-only bracketing, rather than hiding or rounding it down).

## Honest gaps

- **n=1 trial, n=1 (partially reconstructed) stride.** No repeated-stride ensemble is available in this
  trial to separate genuine cycle-to-cycle gait variability from this method's own measurement noise.
- **The full-stride reconstruction is a disclosed approximation, not raw data.** The file truncates before
  the next left heel-strike; the missing quarter-phase is completed by mirroring the measured single-limb-
  support segment (justified by this trial's own 0.89%-symmetric double-support timing), and shown (via
  the reversed-proxy specificity check) to be sensitive to getting the mirroring direction right — a real
  limitation of single-trial completion, not hidden behind the headline recovery number.
- **Double-support fraction (33.8%) vs. the commonly-quoted "~20%" figure**: the gap is explained by this
  trial's slower-than-typical 1.06 m/s speed via the well-established inverse speed/double-support
  relationship, but the exact quantitative slope of that relationship was not extracted from a live-fetched
  primary-source table this session — `Öberg et al. 1993 (PMID 8035350)` confirms such tables exist
  (abstract only quotes the study's scope, not its numeric tables); a full-text/table fetch would be needed
  to check the precise expected number at this exact speed rather than the directional argument given here.
- **Alexander (1989, PMID 2678167)** is live-verified as a real citation (title/DOI/PMID all confirmed via
  NCBI eutils) for the historical Froude-number-in-locomotion tradition, but MEDLINE carries no abstract
  text for it — the quantitative "~0.5" claim in this doc is instead anchored directly to `Carr & McGee
  (2009, PMID 19672305)`'s abstract, which states the number explicitly; Alexander 1989 is cited for
  intellectual lineage, not as the numeric source.
- **2D sagittal-plane simplification**: mediolateral (Z-axis) GRF and rotational/segmental internal kinetic
  energy are excluded from KE/PE, matching Cavagna's own classic method but still a real simplification.
- **Ankle-to-ground clearance (~0.07 m)** in the leg-length sum is a coarse anthropometric term (not
  extracted from the model's own foot mesh geometry) — a minor (~6.6%) contributor to the 1.0645 m total;
  the two independent leg-length derivations agreeing to 0.19% bounds how much this term could be wrong.
- **Mass estimates differ by 1.33%** (established repo value 78.2 kg vs. this session's static-trial F=mg
  measurement 79.24 kg) — small and physically plausible (footwear/force-plate calibration offset), not
  reconciled to sub-percent precision this session; both are carried through the analysis in parallel
  (Part 1's table) rather than one being discarded.
- **Dysfunction layer is conceptual, not measured this session.** No pathological-gait trial (antalgic,
  Trendelenburg, spastic/hemiparetic, Parkinsonian festination) exists in this repo's read-accessible data;
  the mapping below to this doc's own measured variables is a proposed_cell, not a claim.

## Function → dysfunction: the video-measurable deviation layer (conceptual mapping, not measured this session)

This doc's twin-measured variables are exactly the ones a markerless-video pipeline could extract and
compare against, per the operator's framing (low-functioning gait as deviation from this validated
baseline), each mapping to a specific, named pathology:

| variable this doc measured (function baseline) | dysfunction it would flag as a deviation |
|---|---|
| DS1≈DS2 (0.89% symmetric); stance duration equal L/R | **antalgic gait** — shortened stance time on the painful limb, DS-timing asymmetry rises |
| M-shape push-off peak (1.076×BW) driven by DS2 active work (Part 1) | **Parkinsonian festination** / weak push-off — the second (push-off) peak flattens or disappears |
| step width / frontal-plane stance stability (not measured this doc — sagittal-plane only) | **Trendelenburg gait** — hip-abductor weakness → contralateral pelvic drop during single support |
| stance-phase smoothness / GRF shape regularity | **spastic/hemiparetic gait** — irregular, non-two-peaked stance GRF, asymmetric stride time |

## couples_to

- `docs/MECHANISM_SPRINT_SPRING_MASS.md` (+`_evidence.json`) — the direct running/spring-mass counterpart;
  shares the Froude-number/dynamic-similarity vocabulary (`Vaughan 2005`, and now `Carr & McGee 2009`) and
  is joined at the walk-run transition this doc's Part 3 sits just below.
- `docs/MECHANISM_GAIT_KINEMATICS_FIDELITY.md` — already-certified (16/17 gates PASS) that subject2_walking1
  is normal, healthy, textbook-consistent gait; the precondition this doc's function-baseline framing relies
  on rather than re-derives.
- `docs/MECHANISM_CONTACT_WAVEFORM.md`, `docs/MECHANISM_PUSHOFF_PLANTARFLEXOR.md` — share the stance-phase
  GRF-shape vocabulary; this doc's DS1(collision)/DS2(push-off) mechanistic split (Part 1) is directly
  relevant context for `PUSHOFF_PLANTARFLEXOR`'s finding that contact-force over-prediction concentrates at
  the push-off hump specifically.
- `docs/MECHANISM_METABOLIC_SPEED_CURVE.md`, `docs/MECHANISM_METABOLIC_CALORIMETRY.md` — the metabolic side
  of the same walking-economy story: this doc quantifies the **mechanical, passive** share (52–80% of CoM
  energy fluctuation recovered pendulum-style); those docs quantify the **muscular/metabolic** remainder.
- `docs/MECHANISM_SPINE_GAIT_VBR.md` — general gait-force sibling, different mechanism layer (internal
  segment loading vs. this doc's whole-body CoM energetics).
- `MSK-GAIT-LOCOMOTION` (anchor graph node, OPEN/SEED-DESIGN, alias
  `AUTO-TEST-THE-CORE-REDUCED-REP-GENERALIZES-CL`) — proposes the inverted-pendulum/SLIP low-dim-latent
  parametrization for a different (reduced-rep-generalization) hypothesis; this doc supplies a
  twin's-own-data-validated forward pendulum model as an input, not that node's full ambition.
- The function↔dysfunction organizing axis: this doc's contribution is the validated **function** pole
  (Part 1–3) plus the explicit conceptual dysfunction mapping above.

## Proposed cell

Run the twin's markerless video pose pipeline on gait footage (starting with whatever is already indexed
in `data/video_index_tables/` and `data/video_sources/`) to extract stance-time asymmetry, step width, and
a GRF-shape proxy (e.g., vertical CoM/pelvis acceleration pattern) from video alone; machine-gate the
extracted values against this doc's own pre-registered numbers (33.8% double support, 1.011×/1.076×BW
M-shape peaks, Fr=0.108) rather than eyeballing footage, to detect antalgic/Trendelenburg/spastic/
festination deviations as quantified departures from this now-measured normal-function baseline. If any
pathological-gait trial ever becomes available in this repo's read-accessible data, re-run
`scripts/msk/gait_mechanics.py` against it directly for a real (not conceptual) function-vs-dysfunction
contrast.

## Repro

```
cd ~/projects/bodytwin
python3 scripts/msk/gait_mechanics.py   # writes reports/probes/gait_mechanics_results.json
```
Pure numpy/scipy (peak-finding only) — no OpenSim API needed at run time (`.mot`/`.osim` are parsed as
plain text). Reads the raw external-drive force-plate files read-only (verified mounted at
`/media/anton/8838D60F38D5FBDE/...`); if that drive is not mounted, the script fails loudly (file-not-found)
rather than silently falling back to synthetic data. Verified working under plain system `python3`
(numpy 2.2.6/scipy 1.15.3). Runtime: ~2 seconds. No git operations; no writes outside
`reports/probes/gait_mechanics_results.json` and this doc pair (this `.md` + the evidence JSON).
