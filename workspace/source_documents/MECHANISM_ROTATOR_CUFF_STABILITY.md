# MECHANISM ROTATOR CUFF STABILITY — concavity-compression + coronal force-couple (2026-07-22)

Complementary to `docs/MECHANISM_GLENOHUMERAL_FORCE.md` (joint reaction **magnitude** via Static
Optimization) and `docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md` (a crossing-muscle sign fix).
Neither computes joint-reaction magnitude — this cert answers the **stability/direction**
question: does the net muscle force vector on the humerus point INTO the glenoid concavity
(compressive, stable) or ACROSS its face (shear, unstable), and does adding the rotator cuff to a
deltoid-only pull change the answer. The underlying script,
`scripts/msk/validate_rotator_cuff_stability.py`, already existed and had already been run once —
but this session found that prior run **stale and, in its own most decisive leg, silently never
executed at all** (§1), root-caused why (§2-3), forced the adversary through on the corrected
substrate (§4), and anchored the result against two independently-verified cadaveric papers plus
a first-principles geometric derivation bridging them (§5-6). Every number below is either
reproduced bit-for-bit from the pre-existing repo JSON (unmodified, `.venv-msk`, read in place) or
newly machine-measured this session; no existing repo file was edited or overwritten.

## Headline result

| | | |
|---|---:|---|
| Existing on-disk result (`rotator_cuff_stability_results.json`) | **STALE** | predates 2 keys its own current script computes (§1) |
| Twin model (subject2) direct test, 4 poses: does cuff improve shear/compression ratio? | **0/4 — FAILS at every pose** | counter-anatomical (§2) |
| Root cause (machine-queried, not inferred) | **11/11** deltoid+cuff origins parented to `torso`, 0/11 to `scapula_r` | confirms an already-disclosed gap (§3) |
| Forced adversary on anatomically-valid (donor) substrate: deltoid-ALONE ratio | **3.25** (unstable, +4223.7 N superior shear) | worst point in the whole grid |
| Forced adversary: deltoid+CUFF ratio (same pose, max activation) | **1.01** | **3.21× improvement** |
| Monotonicity across the 4×5 activation grid (20 points) | **20/20 monotonic**, cuff strictly helps at every deltoid level | PASS |
| Concavity-compression stability ratio, labrum intact (Lippitt et al. 1993, PMID 22959294, live-verified) | **up to 60%** of compressive load | external anchor |
| Glenoid+labrum socket depth (Howell & Galinat 1989, PMID 2721051, live-verified) | **9 mm SI / 5 mm AP**, labrum ≈50% of depth | external anchor |
| Geometric bridge (this session, ratio = tan θ_max) cross-paper prediction | AP stability drops **32.9%** after Bankart-type detachment | vs. Lippitt's own independent ~20% (order-of-magnitude match) |
| Dysfunction pole: migration correlates with cuff-tear size (Keener et al. 2009, PMID 19487518, live-verified) | tear ≥175 mm², p=0.0002 | external anchor |

**CONFIDENCE TIER: the force-couple/concavity-compression mechanism is B-grade CONFIRMED on the
donor (anatomically-valid) substrate, externally anchored twice over (Lippitt ratio, Howell
depth) plus a derived geometric bridge between them. The twin model's OWN result is a DIAGNOSED
model-fidelity gap (not a refutation) — subject2's own digital twin cannot currently be used to
test this specific mechanism.**

## 0. Why this doc didn't already exist

`grep`-ing `docs/` for rotator-cuff/concavity/force-couple certs found `scripts/msk/
validate_rotator_cuff_stability.py` and a results JSON already on disk
(`data/msk_smoketest/subject2_rotator_cuff_stability/rotator_cuff_stability_results.json`), but no
`docs/MECHANISM_ROTATOR_CUFF_STABILITY.md` — the script had been written and run once, but never
written up, and (per §1) never actually finished producing its own most important output.

## 1. Staleness + a live crash bug caught before trusting anything

**mtime check**: the script (`08:11:02`) is newer than its own results JSON (`08:05:39`) — a
direct, dated signal the JSON might not reflect the script's current behavior, checked before
reusing anything. Confirmed by content: the on-disk JSON has exactly 4 top-level keys (`poses`,
`axis_stability_check`, `neutral_frame_axes`, `cuff_activation_sweep_at_30deg`) — **missing**
`twin_origin_clustering_check` and `donor_leg`, both of which the CURRENT script computes and
writes. The script's own docstring calls `donor_leg()` **"the decisive, diverse-instance-space
(activation-grid) test"** — i.e. the most important leg was never actually inspected.

Re-running is the obvious next step, but the isolation rule for this session is "touch only files
you create" (another process may depend on this exact JSON) — so rather than overwrite it, the
module was **imported read-only** (`sys.path` insertion, no edits) into a new wrapper script, and
its own functions were called fresh, writing only to a new scratch path.

**Reproducibility check first** (a machine cross-check, not blind trust): the fresh run's
`deltoid_only_ratio` at all 4 poses matched the stale JSON's own values to full float64 precision
(e.g. ~30deg: 2.1310257772545826 both runs) — confirming the model is unchanged and the stale
JSON's *existing* content is trustworthy as far as it goes; only the *missing* content was in
question.

**Calling `donor_leg()` directly reproduces a real crash**: `AssertionError:
DeltoideusClavicle_A: last point not on humerus_r`. Root cause: `muscle_pull_force`'s assertion
hardcodes the distal body name as `"humerus_r"`, but the donor model
(`ThoracoscapularShoulderModel.osim`, a single non-lateralized shoulder-only model — confirmed via
its own load-time geometry warnings, `humerus.vtp`/`ulna.vtp`, no `_r` suffix) names that body
`"humerus"`. **This function could never have produced the missing JSON key on this donor model —
it always crashes first.** This is why "stale" was not just a timing accident: the decisive leg
had a live bug, not merely an unrun one.

## 2. Twin-model (subject2) direct test — an apparent, counter-anatomical falsification

Reused (bit-for-bit, reproduced) from the existing JSON — deltoid = `DeltoideusClavicle_A,
DeltoideusScapula_{M,P}`; cuff = `Supraspinatus_{A,P}, Infraspinatus_{I,S}, Subscapularis_{I,M,S},
TeresMinor`; both at full (1.0) activation, "strongest fair form" per muscle's own real
`max_isometric_force`. Decomposed along the scapula-parent-offset-frame's own local axes
(compression = −Z/medial, shear = X-Y/anterior-superior plane):

| pose | measured total elevation | deltoid-only ratio | deltoid+cuff ratio | cuff made it... |
|---|---:|---:|---:|---|
| ~0deg | 0.00° | 5.694 | **+Infinity** (compression went **negative**, −511.7 N) | **worse** |
| ~30deg | 29.99° | 2.131 | 37.440 | **worse** |
| ~60deg | 59.99° | 1.068 | 5.115 | **worse** |
| ~90deg | 90.00° | 0.487 | 2.239 | **worse** |

**4/4 poses: adding the cuff makes the shear/compression ratio worse, not better** — the opposite
of the textbook force-couple/concavity-compression direction. A supplementary cuff-activation
sweep at the ~30deg pose sharpens this: as cuff activation rises 0→1.0, compression falls
monotonically (1906→183 N) and the ratio rises monotonically (2.13→37.44) — a *systematic*,
not noisy, wrong-direction signature.

Per this repo's own discipline, a one-shot counter-anatomical result is **not** accepted as an
honest negative without forcing the Orient step first.

## 3. OODA — root-causing the twin's counter-anatomical result (Orient, not hand-waved)

**Observe**: cuff addition systematically worsens the ratio at all 4 poses, sharpest at 0°.

**Orient — is the reference frame or the muscle geometry at fault?** Two independent, direct
machine checks, not narrated:

1. **Reference-frame check (already in the existing JSON, re-verified)**: the compression/shear
   axes (scapula-parent offset frame) rotate with the scapula's own `scapula_upward_rot_r` to
   within **0.173°** over a 30.2° scapular rotation (0.57% relative) — **PASS, the frame is not
   the problem.**
2. **Muscle-origin body-attachment ground truth (new this session, direct API query, not a
   proxy)**: for every one of the 11 deltoid+cuff muscles, `PathPoint.getBody().getName()` on the
   proximal ("origin") point:

   | | origin_body | insertion_body |
   |---|---|---|
   | DeltoideusClavicle_A, DeltoideusScapula_{M,P}, Supraspinatus_{A,P}, Infraspinatus_{I,S}, Subscapularis_{I,M,S}, TeresMinor (**11/11**) | **torso** | humerus_r |

   **0/11 attach to the new `scapula_r` body.** A positive control (3 of the 11 *scapular-
   stabilizer* muscles re-routed in a separate, later session,
   `docs/MECHANISM_SCAPULA_CLAVICLE.md`) shows mixed torso/scapula_r attachment at their own path
   points — confirming the query itself correctly discriminates the two bodies; it is not
   silently returning a constant.

**Decide/Act**: since `scapula_r` rotates ~30° relative to `torso` across this exact pose sweep
(the frame-stability check above), every deltoid/cuff muscle's geometrically-computed pull
direction is being evaluated from a **rigidly torso-fixed** origin that does **not** follow the
real scapula through scapulohumeral rhythm — at any non-neutral pose, the computed pull direction
is anatomically wrong. This is not a new mystery: `docs/MECHANISM_SCAPULA_CLAVICLE.md` §8 item 4
already disclosed "several of those 25 muscles' proximal points are still collapsed onto 'torso'
rather than re-attached to the new scapula_r/clavicle_r — an explicit, disclosed opportunity NOT
taken." This session supplies the direct, quantitative confirmation that disclosure predicted.
**Verdict: the twin-model negative is a diagnosed model-fidelity gap, not a refutation of
concavity-compression/force-couple stability** — the adversary must be forced on a substrate
where this specific gap does not exist.

## 4. Forced adversary on the anatomically-valid substrate — the decisive test

`ThoracoscapularShoulderModel.osim` (Seth, Dong, Matias, Delp 2019, PMID 31780916 — the same donor
`add_scapula_clavicle.py` already sourced): a real, separate scapula body, and (checked live) a
roughly isotropic origin spread across the 9 deltoid/cuff-relevant muscles (9.0 / 9.6 / 9.1 cm in
x/y/z — consistent with genuinely distinct scapular-fossa origins, unlike the twin's skewed
14.1/8.2/3.9 cm).

**Bug found and fixed before trusting any output** (§1's crash): `muscle_pull_force`'s
`assert ins.getBody().getName() == "humerus_r"` fails on this donor model outright. Fixed, in a
new script (repo script untouched), by resolving the distal body name **empirically** — querying
what body all 11 deltoid+cuff muscles' own last path point actually reports (unanimous, 11/11:
`"humerus"`) rather than hardcoding a side-specific name. Same algorithm, generalized, not a
different method.

At the donor's own default pose (`plane_elv≈0°, shoulder_elv=18.43°, axial_rot=58.95°` — a single
fixed pose, not swept, disclosed as a scope limit in §7):

| | ratio (shear/compression) | compression | superior shear | compression positive? |
|---|---:|---:|---:|:---:|
| **Deltoid ALONE**, max activation | **3.251** | 1305.1 N | **+4223.7 N** | yes |
| **Deltoid + CUFF**, max activation | **1.012** | 5237.6 N | +3453.6 N | yes |

**Deltoid-alone is the single worst point in the entire grid below** — large superior shear, weak
compression relative to it: the geometric signature of "the head rides up," the measured cuff-tear
phenotype this cert's own falsifier required as the adversary. Adding full cuff activation
improves the ratio **3.21×** and cuts superior shear by **18.2%** (it does not fully reverse it to
net-inferior at full bilateral co-activation — physiological in-vivo tasks do not usually demand
both muscle groups at 100% MVC simultaneously; the *dose-response direction* is the load-bearing
claim, not perfect cancellation at this one extreme).

**Diverse-instance-space test (the pre-registered gate)**: full 4×5 grid, deltoid activation ∈
{1.0, 0.75, 0.5, 0.3} × cuff activation ∈ {0.0, 0.25, 0.5, 0.75, 1.0} — **20/20 points**: for every
fixed deltoid-activation level, the ratio **strictly decreases** as cuff activation rises from 0
to 1.0 (e.g. at deltoid=1.0: 3.251 → 1.852 → 1.354 → 1.127 → 1.012). The forced adversary
("deltoid alone is fine, cuff adds nothing") **FALLS uniformly across the diverse activation
grid** — not a cherry-picked single comparison.

One shared, disclosed confound checked and ruled out as *differential*: `DeltoideusScapula_M`'s
`max_isometric_force` (2597.8 N — already flagged in `docs/MECHANISM_GLENOHUMERAL_FORCE.md` §4b as
large vs. typical 500-1200 N single-deltoid-head values) is **byte-identical between twin and
donor** (confirmed live) — i.e. it is the donor's own native, unscaled value, present equally in
both legs of this comparison, not a twin-specific distortion inflating one side only.

## 5. External anchors — live-verified via NCBI eutils this session, not recalled

**Lippitt SB, Vanderhooft JE, Harris SL, Sidles JA, Harryman DT 2nd, Matsen FA 3rd (1993).
"Glenohumeral stability from concavity-compression: A quantitative analysis." *J Shoulder Elbow
Surg* 2(1):27-35. PMID 22959294, doi:10.1016/S1058-2746(09)80134-1.** (10 cadaver shoulders, 50N
then 100N compressive load perpendicular to the glenoid, tangential force to dislocation measured
in 8 directions 45° apart, repeated after labral excision.) Quoted directly: **"the humeral head
resisted tangential forces of up to 60% of the compressive load"** with the labrum intact; **"the
degree of compression stabilization varied around the circumference... attributed to the greater
glenoid depth"** superiorly/inferiorly; **"resection of the glenoid labrum reduced the
effectiveness of compression stabilization by approximately 20%."**

**Howell SM, Galinat BJ (1989). "The glenoid-labral socket. A constrained articular surface."
*Clin Orthop Relat Res* (243):122-5. PMID 2721051.** (25 cadaver shoulders.) Quoted directly:
**"approximately 9 mm deep in the superoinferior (SI) direction and 5 mm deep in the
anteroposterior (AP) direction. The... labrum contributes approximately 50%... Detachment of the
labrum anteriorly, as in a Bankart lesion, may reduce the depth... from approximately 5.0 to
2.4 mm."**

An initial PMID for the Lippitt paper looked suspicious (22959294 — a number in a range typically
associated with 2012-era publications, for a nominally-1993 paper) and was **not** trusted at face
value; `esummary` resolved this: `epubdate: "2009 Feb 19"` — the paper was retro-digitized into
PubMed's electronic system in 2009, explaining the high PMID for an actual 1993 publication.
Confirmed genuine by the full author list, title, journal, volume/issue/pages all matching exactly
— not a coincidental keyword hit.

## 6. Geometric bridge — deriving *why* depth governs the ratio, then cross-checking two papers against each other through it

**Derivation (frictionless rigid-body contact equilibrium — geometric, not curve-fit)**: model the
glenoid+labrum concavity as a spherical cap of radius of curvature R_g and rim half-angle θ from
its own center of curvature. Concavity-compression is explicitly a *geometric*, friction-
independent mechanism (Lippitt & Matsen's own framing) — the socket can exert only a contact
*normal* force on the head, directed along the local radius through the current contact point. A
rigid body in equilibrium under exactly two force systems (the applied resultant of compression +
tangential load, and the single contact-normal reaction) requires those two to be exactly
colinear and opposite — so the applied resultant's angle-from-axis must equal whatever angle φ the
current contact point sits at (φ=0 at the pole, φ=θ at the rim). As the tangential/compressive
ratio rises, φ rises to match; the **critical (dislocating) ratio is reached exactly when φ = θ**,
i.e.:

**stability ratio = tan(θ), with cos(θ) = 1 − depth/R_g** (the standard spherical-cap sagitta
relation).

This single relation, with no fitting, predicts: (a) ratio rises monotonically with depth — Lippitt's
own "greatest [stability] magnitude... attributed to greater glenoid depth," confirmed; (b) the
ratio should be a pure function of geometry, **independent of the applied load magnitude** —
consistent with (not separately broken out as two numbers in) Lippitt's one reported "up to 60%"
figure holding across both the 50N and 100N test conditions; (c) reducing depth (labral loss)
reduces the ratio — both papers agree qualitatively.

**Cross-paper over-determination (the genuine test — neither paper's number was used to fit the
other)**: back-solve the ONE unknown, R_g, from Lippitt's SI ratio (0.60) + Howell's SI depth
(9 mm) alone: θ_max = arctan(0.60) = 30.96°, **R_g = 63.2 mm**. Then, **using that same R_g**,
predict what Howell's *independently-measured* AP-direction depths (5.0 mm intact → 2.4 mm
post-Bankart) imply: ratio 0.423 → 0.284, a **32.9% predicted AP-specific stability-ratio
reduction**. Lippitt's own, entirely independently-measured, all-directions-averaged reduction
after full labral excision was **~20%**. These are not expected to match exactly (Howell's number
is an anterior-only Bankart-type detachment; Lippitt's is a full circumferential excision averaged
over 8 directions) — but landing in the same order of magnitude, out-of-sample, through a
first-principles relation calibrated on a *different* direction's numbers, is a real,
non-tautological consistency check, not a coincidence manufactured by curve-fitting. (R_g=63.2mm
itself is back-solved, not independently measured from a third source — disclosed in §7, not
claimed as a validated anatomical radius-of-curvature.)

## 7. Dysfunction pole — cuff tear → loss of the inferior force-couple → superior migration

**Weiner DS, Macnab I (1970). "Superior migration of the humeral head. A radiological aid in the
diagnosis of tears of the rotator cuff." *J Bone Joint Surg Br* 52(3):524-7. PMID 5455085.**
Live-verified bibliographically (title/authors/journal exact match); no digitized abstract exists
(pre-abstract-database era) — the classic original description of the phenotype this cert's
adversary (§4, deltoid-alone) reproduces mechanistically (large superior shear, weak relative
compression).

**Keener JD, Wei AS, Kim HM, Steger-May K, Yamaguchi K (2009). "Proximal humeral migration in
shoulders with symptomatic and asymptomatic rotator cuff tears." *J Bone Joint Surg Am*
91(6):1405-13. PMID 19487518.** Live-verified with full abstract (n=160: 98 asymptomatic + 62
symptomatic cuff-tear shoulders). Quoted: migration **"correlates with rotator cuff tear
size"** (multivariate p=0.01); **"tears extending into the infraspinatus tendon are associated
with greater humeral migration than... isolated supraspinatus tears"** (p=0.01/0.03) — directly
consistent with this cert's own model, in which infraspinatus is one of the two largest
compression-restoring cuff contributors (§4); symptomatic > asymptomatic migration (p=0.03); a
critical tear-size threshold (≥175 mm²) below which neither pain nor size significantly predicts
migration.

**Inman VT, Saunders JB, Abbott LC. "Observations of the function of the shoulder joint." 1944.
*Clin Orthop Relat Res* 1996;(330):3-12 (Classic Article reprint of *J Bone Joint Surg*
1944;26:1-30). PMID 8804269.** Bibliographic fields (title/authors/journal/DOI) live-verified; **no
digitized abstract available via NCBI efetch** (empty attributes list) — cited for the
well-established, extensively secondary-attributed coronal deltoid/cuff force-couple and
scapulohumeral-rhythm concept this cert's own model (§4) numerically instantiates, not for a
specific quoted figure.

**Burkhart SS (1992). "Fluoroscopic comparison of kinematic patterns in massive rotator cuff
tears. A suspension bridge model." *Clin Orthop Relat Res* (284):144-52. PMID 1395284.**
Live-verified bibliographically; noted for completeness (transverse-couple/global-balance framing)
but not deeply used this session — see §8's honest gap on the transverse couple.

## 8. Honest gaps (full list)

1. **Twin model (subject2's own digital twin) cannot currently test this mechanism at all** —
   11/11 deltoid+cuff origins mis-parented to `torso` instead of `scapula_r` (root-caused, not
   fixed, this session — the fix is a graft-registration task on `add_arm_muscles.py`, out of a
   verification session's scope).
2. **Donor-model test uses a single fixed pose** (its own default/rest coordinates), not a
   multi-pose sweep like the twin's own 4-pose test.
3. **Transverse couple (subscapularis vs. infraspinatus/teres minor) — bonus check run,
   INCONCLUSIVE**, not chased further (LEAN scope call): at the donor's default pose, all 6
   muscles showed the SAME-sign anterior/posterior component (no opposition found), but two
   undiagnosed confounds are disclosed rather than either accepted as a kill or hidden: (a) the
   donor's own X/Y axes were never independently anatomically verified the way the twin's were;
   (b) the donor's default pose carries a large (58.9°) axial rotation whose provenance is
   unexamined, and the textbook transverse-couple description assumes near-neutral rotation.
4. **The back-solved glenoid radius of curvature (63.2 mm, §6)** is implied by the two anchor
   papers via this session's own geometric relation, not independently measured from a third
   source.
5. **`DeltoideusScapula_M`'s Fmax (2597.8 N)** is large vs. typical Delp-lineage single-head values
   (already disclosed elsewhere) — shared identically by both legs compared here, so not a
   differential confound for the *ratio* comparison, but it inflates absolute force magnitudes in
   both.
6. **No wrapping surfaces** on any deltoid/cuff muscle in either model (straight-line paths) —
   unchanged, previously-disclosed gap.
7. **Right arm only; static synthetic activation grid, not EMG-driven** — activation levels are a
   prescribed 0-1 grid, not measured muscle activity.
8. **The donor's own default-pose provenance** (why `axial_rot` defaults to 58.9°) was not
   investigated this session.

## Falsifier verdict

- **(a) Concavity-compression stability ratio vs. glenoid depth**: PASS, via external anchor +
  derived geometric bridge (§5-6) — genuine cross-paper over-determination (Lippitt's ratio,
  Howell's depth, jointly consistent through a first-principles relation, out-of-sample AP
  prediction lands in the right order of magnitude). Not a from-scratch re-derivation of the
  cadaveric experiment itself (no contact-mechanics substrate exists in this repo for the
  glenohumeral joint — checked, only muscle-driven kinematic models exist).
- **(b) Deltoid-cuff coronal couple keeping the resultant centered**: PASS on the anatomically
  valid (donor) substrate — deltoid-alone ratio 3.25 (the required-to-fail adversary, confirmed
  failed: worst point in the whole grid, large +superior shear); deltoid+cuff ratio 1.01 (3.21×
  improvement), monotonic across all 20 grid points. FAILS on the twin (subject2) model at all 4
  poses, but that negative was forced through OODA to a specific, machine-confirmed root cause
  (11/11 origins on the wrong body) rather than accepted as a refutation — a diagnosed
  model-fidelity gap, not a falsification of the mechanism.

## Files

- `scripts/msk/validate_rotator_cuff_stability.py` — pre-existing, **not modified**; imported
  read-only this session.
- `data/msk_smoketest/subject2_rotator_cuff_stability/rotator_cuff_stability_results.json` —
  pre-existing, **not modified**; its `poses`/`axis_stability_check`/`cuff_activation_sweep_at_30deg`
  content reused after a bit-for-bit reproducibility check; its absence of `donor_leg` explained
  (§1) and supplied fresh (§4), written only to this session's own new evidence file.
- `docs/MECHANISM_ROTATOR_CUFF_STABILITY_evidence.json` — **new this session**: twin-model
  reproduction, root-cause body-attachment query, fixed donor-leg activation grid, transverse
  bonus check, geometric-derivation numbers, all 6 literature anchors with PMIDs/quotes.
- `docs/MECHANISM_SCAPULA_CLAVICLE.md` §8 item 4 — the prior disclosure this session's root-cause
  finding directly confirms.
- `docs/MECHANISM_GLENOHUMERAL_FORCE.md` §0, §4b — the crossing-muscle set and the
  `DeltoideusScapula_M` oversized-Fmax disclosure this session cross-checked as a shared
  (non-differential) confound.
