# MECHANISM JOINT FORCE VALIDATION SCORECARD — synthesis across 6 joints (2026-07-21)

Integrates every corrected joint-force cert in this family (knee, hip, ankle, spine, shoulder,
elbow) into one honest scorecard. **Every number below was re-verified this session directly
against the on-disk `*_results.json` machine output of the cert that produced it** — not
transcribed from doc prose — before being placed in a table (files listed in each row's source;
full list in §7). One genuine cross-check finding surfaced by this process is reported in §6, not
hidden. Isolation respected: bodytwin only, read-only (no data written except this doc), no git
commit/push.

## 1. Master scorecard

| Joint | Tier-1 twin: pure-kinematic reaction (%BW) | Tier-2 twin: muscle-driven, corrected (%BW) | In-vivo / published anchor (%BW) | Ratio (Tier-2 / anchor) | Anchor TIER |
|---|---:|---:|---:|---:|---|
| **Knee** | 102.17 (ratio 0.396) | self-computed **391.10** / JointReaction **391.11** | OrthoLoad knee median **258.22** (n=72 trials, 9 subj, unassisted level walking) | **1.515×** | in-vivo instrumented — full contact-bearing implant |
| **Hip** | 88.59 (ratio 0.323) | self-computed **386.77** / JointReaction **387.04** | OrthoLoad hip median **273.93** (n=162, 18 subj) | **1.412×** | in-vivo instrumented — full contact-bearing implant |
| **Ankle** | 107.50 (ratio 0.226) | self-computed **484.49** / JointReaction **484.48**; Achilles (scalar sum) 320.26 | Cadaveric/FE literature band **390–540**, mid **476.7** (talocrural-only mean 445); Achilles anchor 390 (Giddings 2000) | **1.016×** | **cadaveric/model-literature plausibility — no in-vivo program has ever existed for this joint** |
| **Spine** | no_box 41.19 / with_box 52.78 | no_box **296.59** / with_box **389.17** | OrthoLoad spine_vbr stoop-lift-10kg median **170.94** (n=13, 4 subj) | **1.74× (no_box) – 2.28× (with_box)** | in-vivo instrumented, but a **partial-load-sharing implant** (VBR + posterior fixation shunts load away — structurally reads a fraction of true compression, not the whole joint) |
| **Shoulder** | native model: **5.306** (ratio 0.073) — **0/80 native muscles cross the joint; SO is structurally impossible (0 redundancy)** | non-native donor-grafted model: self-computed **19.065** / JointReaction **19.362** | OrthoLoad shoulder median **72.2** (n=23, 6 subj, unloaded 90° abduction) | **0.264×** (closes 20.6% of the pre-graft gap) | anchor is in-vivo instrumented (strong); twin side is a **diagnosed structural gap** (native) or a **non-native fidelity-limited stand-in** (grafted) |
| **Elbow** | **8.61** at 5 kg load (= bone-contact identically, 0 muscles to subtract); twin's own 10 N·m ideal actuator is **inadequate** at 5 kg (required 18.14 N·m, 1.81× over capacity) | **none possible** — 0/80 muscles cross the elbow at any of 3 tested poses (whole-arm structural gap) | **none exists** — no in-vivo elbow program anywhere in the literature. Non-twin illustrative demo (arm26/Holzbaur 2005): 455.2 N = 9.28× load weight = 59.4%BW (borrowed reference mass) | **N/A** — no anchor to ratio against | **diagnosed-gap on both sides** (no twin estimate, no anchor) — weakest joint in the family |

Source JSONs re-verified this session (machine-checked, not eyeballed): `data/msk_smoketest/subject2_walking1/joint_force_validation/joint_force_validation_results.json`,
`.../sign_bug_audit/sign_bug_audit_results.json`, `.../hip_force_validation/hip_force_validation_results.json`,
`.../ankle_force_validation/ankle_force_validation_results.json`,
`data/msk_smoketest/subject2_spine_stoop_lift/spine_force_validation_results.json`,
`data/msk_smoketest/subject2_arm_abduction/shoulder_force_validation/shoulder_force_validation_results.json`,
`.../shoulder_force_with_muscles/shoulder_force_with_muscles_results.json`,
`data/msk_smoketest/subject2_elbow_static/elbow_force_validation/elbow_force_results.json`.

### Anchor tier legend (what the operator asked to make explicit)

- **In-vivo instrumented, full contact-bearing** (knee, hip): implant telemetry measures the actual
  bone-on-bone/prosthesis contact force during real human movement — the strongest evidence class
  in this family.
- **In-vivo instrumented, partial-load-sharing** (spine): real implant telemetry, but the
  instrumented segment structurally shares axial load with adjacent hardware (posterior
  pedicle-screw/rod fixation) — a genuine in-vivo measurement of a *different* (smaller) physical
  quantity than "total joint compression," not a weaker measurement of the same quantity.
- **In-vivo instrumented anchor, twin-side gap** (shoulder): the *anchor* is full-tier quality
  (glenohumeral implant telemetry), but the *twin's own prediction* is what's degraded — either
  structurally absent (native model) or built on a non-native, fidelity-limited donor graft.
- **Cadaveric/model-literature plausibility** (ankle): published quasi-static/cadaveric/
  finite-element estimates (1977–2000), not a direct in-vivo instrumented measurement — no such
  program has ever existed for this joint.
- **Diagnosed-gap** (elbow, and shoulder's native side): the twin cannot produce a muscle-driven
  estimate at all — a structural absence of redundant actuation, not a convergence/tuning failure.
  Elbow additionally has no anchor of any kind.

## 2. The unifying finding — and precisely where it does and doesn't hold

**At the two joints with a full-contact-bearing in-vivo anchor and a fully-converged, corrected
muscle-driven estimate (knee, hip), the twin over-predicts in-vivo joint contact force by
~1.4–1.5× (knee 1.515×, hip 1.412×).** This is not a self-computed-vs-JointReaction disagreement —
the two structurally independent methods now agree with each other to <0.1% at both joints
(0.0039% knee, 0.070% hip), so the over-prediction is a property of the **method architecture**,
not noise or a residual code bug. Verified live (machine-fetched, verbatim-quoted, not recalled):

> Moissenet F, Chèze L, Dumas R (2014). "A 3D lower limb musculoskeletal model for simultaneous
> estimation of musculo-tendon, joint contact, ligament and bone forces during gait." *J Biomech*
> 47(1):50-58. DOI [10.1016/j.jbiomech.2013.10.015](https://doi.org/10.1016/j.jbiomech.2013.10.015),
> PMID [24210475](https://pubmed.ncbi.nlm.nih.gov/24210475/): *"Musculo-tendon forces and joint
> reaction forces are typically estimated using a two-step method, computing first the
> musculo-tendon forces by a static optimization procedure and then deducing the joint reaction
> forces from the force equilibrium... the joint reaction forces are usually overestimated."*

This twin's pipeline at every joint is exactly this two-step pattern (`R` computed from
kinematics+GRF first, SO muscle tensions computed separately, then `F_bone_contact = R −
Σmuscle-crossing`) — not a twin-specific defect, a documented tendency of the method class.

**This does NOT generalize uniformly to all 6 joints — that is the point of building a scorecard
rather than quoting one number:**

- **Ankle (1.016×, lands almost exactly on its anchor) is the exception, not a counter-example**
  (the remediation doc's own words, confirmed here): its anchor is the weaker cadaveric-literature
  tier, not in-vivo — there is no in-vivo signal to over-predict against in the first place.
- **Spine (1.74–2.28×, over-predicts by MORE than knee/hip) was pre-registered, before running, to
  read higher for a structurally different reason**: a VBR implant shares load with posterior
  fixation, so the anchor structurally under-reads true total compression — not (necessarily, or
  not only) the Moissenet mechanism. **Open confound this scorecard surfaces and no child doc
  resolves**: spine's Tier-2 uses the *identical* `R − Σmuscle-crossing` two-step architecture
  Moissenet flags as an over-predictor at knee/hip, so some fraction of the spine's 1.74–2.28× gap
  could be the *same* method-architecture over-prediction stacked on top of the load-sharing
  effect, not purely load-sharing. `docs/MECHANISM_SPINE_FORCE.md` §6 itself discloses the
  load-sharing correction is "qualitative... not a precisely-calibrated correction factor" but does
  not consider or rule out this second mechanism — nothing measured this session distinguishes the
  two. Flagged as open, not resolved.
- **Shoulder (0.264×) and elbow (no estimate possible) show the OPPOSITE finding** —
  under-prediction or no estimate at all — because the native twin arm has **zero** muscles: the
  Moissenet mechanism requires a redundant muscle set to over-recruit; a model with no redundancy
  to resolve cannot exhibit it. The grafted shoulder model (non-native) narrows but does not close
  this gap, and for reasons unrelated to Moissenet (missing scapular rhythm/wrapping, §4).

## 3. Preserved honest negatives (the operator's explicit checklist, verified against source docs)

1. **RRA didn't clear the residual.** The pelvis residual force peak (172.96 N, 22.6%BW) **FAILS**
   the pre-registered "good" band (<75 N) at knee, hip, and ankle — all three reuse the identical
   underlying Static Optimization run, and no Residual Reduction Algorithm was run first to shrink
   it (`docs/MECHANISM_STATIC_OPT.md` §4/§8.2, inherited unchanged by `MECHANISM_HIP_FORCE.md` §7.6
   and `MECHANISM_ANKLE_FORCE.md` §6). Joint-level (hip/knee/ankle) reserve-actuator usage stayed low
   (≤12.5% of a small optimal_force) throughout, so leg-muscle recruitment is not *obviously*
   contaminated by it — disclosed as non-blocking, not resolved.
2. **Corpus-SO is dynamically implausible.** On the *separate* 228-clip video corpus
   (`docs/MECHANISM_FORCE_SCENES_BATCH.md` — MediaPipe-derived kinematics, no real GRF, no
   per-subject scaling: a materially weaker-tier dataset than subject2's lab mocap), the identical
   Static-Optimization method converged numerically on all 3 representative clips (squat/jump/other)
   but was machine-diagnosed `dynamically_plausible=False` on all 3 — knee 614–685%BW, hip
   770–1044%BW, "far outside any published literature" — **both before AND after the sign-bug fix**
   (re-verified live: squat sign-fixed ankle/knee/hip = 239.06/684.78/1044.16%BW, matching the
   doc's table exactly). The fix made these numbers slightly **larger**, not smaller — proof the
   sign bug was not the (sole) cause; root-caused instead to SO's own light 4-6Hz coordinate lowpass
   amplifying this noisier corpus's joint angles into spurious torque demands. Explicit verdict in
   that doc: "HONEST NEGATIVE for the absolute muscle-driven contact-force number on this corpus... 
   not hidden, not silently reported as if trustworthy." The bug-immune pure-Newton reaction-force
   numbers on that same corpus remain trustworthy; only the muscle-driven absolute numbers are the
   negative. This is a *different dataset* than the rest of this scorecard (subject2 lab mocap +
   real GRF) — kept separate deliberately, not blended into the master table.
3. **Shoulder needs scapula + wraps.** The grafted shoulder model has no scapulothoracic joint (real
   scapulohumeral rhythm — ~1/3 to 1/2 of total arm elevation — is entirely collapsed away) and no
   ported wrapping surfaces, which directly and measurably causes pose-dependent moment-arm **sign
   flips**: only 47/70 (67%) of pre-registered prime-mover sign checks pass; the middle deltoid (the
   textbook prime abductor) has the *wrong* sign at the neutral pose and the *right* sign at 90°
   abducted — a concrete, machine-measured symptom of the missing wrap geometry, not a mystery.
   `docs/MECHANISM_ARM_MUSCLES.md` §8 calls this explicitly "a first-step graft, not a validated
   shoulder model."
4. **Ankle is a weaker cadaveric tier.** No in-vivo instrumented-implant program has ever existed
   for the human ankle (verified absent from OrthoLoad and elsewhere, both by the ankle cert itself
   and independently by `docs/MECHANISM_MSK_BUILD_PLAN.md`'s anchor registry). The anchor is 3
   cadaveric/quasi-static/finite-element studies (1977–2000); one of the three source numbers
   (Giddings' 540%BW) is technically the adjacent talocalcaneal (subtalar) joint, not the talocrural
   joint this cert actually cuts at — disclosed, narrowing the anchor to 445%BW if restricted to the
   two talocrural-only studies.
5. **Sign-bug history.** A real, machine-confirmed sign bug (`static_opt_knee.py`'s
   `knee_crossing_muscles_and_forces`: an index ternary backwards in *both* branches) silently
   inverted the muscle-crossing-force direction at knee/hip/ankle — masked because the headline
   number is a vector **norm** (`np.linalg.norm`, always ≥0, structurally cannot expose an internal
   sign flip). It was caught not at knee/hip (where it produced a plausible-*looking* undershoot,
   ratio 0.86–0.90) but at the **spine**, where a *signed* axial projection turned the identical bug
   into an impossible negative (tensile) force — forcing the discovery
   (`docs/MECHANISM_SIGN_BUG_AUDIT.md` §3 measured this mechanism directly: buggy/corrected
   peak-instant vectors point 157.7°/162.6° apart at knee/hip, the signature of `R+M` vs `R−M`).
   Fixing it flipped the headline finding at knee/hip from "lands just under in-vivo" to "over-predicts
   in-vivo by ~1.4–1.5×" (see the before/after table in §5). A **separately-defined copy** of the
   identical bug was independently found, live, in the shoulder-with-muscles script (not auto-fixed
   by the shared-function patch, since it never called the shared function) and fixed in a dedicated
   follow-up (`docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md`). A **third**, non-twin illustrative
   copy (the elbow cert's arm26 reference comparison) remains flagged, live, low-priority, and
   explicitly out of scope (`docs/MECHANISM_SIGN_BUG_REMEDIATION.md` §7.2) — **not fixed by this
   scorecard either** (isolation: read-only synthesis, no code edits).

## 4. Shoulder / elbow: the model-fidelity ceiling, stated plainly

Both native-model findings (0/80 muscles cross the shoulder OR the elbow, at any tested pose) are
**one whole-arm structural gap**, not two independent joint-specific findings
(`docs/MECHANISM_ELBOW_FORCE.md` §2: "this is a whole-arm gap, not elbow-specific... both would be
diagnosing the same root cause: this model family's arm was built for kinematic arm-swing during
gait IK, not for upper-limb load-bearing dynamics"). The shoulder cert's own graft (donor:
Seth/Dong/Matias/Delp 2019, DOI 10.3389/fnbot.2019.00090, PMID 31780916) demonstrates that adding
real muscles is *tractable* (3.59× increase over the zero-muscle baseline, matched to 1.54% by an
independent JointReaction referee) but explicitly does not, by itself, produce a "validated shoulder
model" — the two missing pieces (scapulothoracic rhythm, wrapping surfaces) are named, concrete,
scoped follow-ons in `docs/MECHANISM_ARM_MUSCLES.md` §9, not vague gaps. The elbow has no analogous
graft attempted in this cert family; its only muscle-driven number comes from a **non-twin**
reference model (arm26/Holzbaur 2005) offered purely to demonstrate the method and the geometric
mechanism (a short ~2–5 cm flexor moment arm vs. a ~30 cm forearm lever → 8–12× force
multiplication) — explicitly disclaimed as "not a twin result" in `docs/MECHANISM_ELBOW_FORCE.md` §6.

## 5. Sign-bug correction, before vs. after (why the "twin over-predicts" headline exists at all)

| joint | published (bug artifact) | corrected | independent JointReaction referee | agreement |
|---|---:|---:|---:|---:|
| Knee | 233.20 %BW (ratio 0.903) | **391.10 %BW (ratio 1.515)** | 391.11 %BW | 0.0039% |
| Hip | 235.30 %BW (ratio 0.859) | **386.77 %BW (ratio 1.412)** | 387.04 %BW | 0.070% |
| Ankle | 286.68 %BW (ratio 0.601) | **484.49 %BW (ratio 1.016)** | 484.48 %BW | 0.0014% |
| Shoulder (grafted) | 24.18 %BW (ratio 0.335) | **19.07 %BW (ratio 0.264)** | 19.36 %BW | 1.54% |

Note the shoulder correction moved **down** (further from its anchor) while knee/hip/ankle moved
**up** — measured, not assumed, to be a property of whether the local muscle-crossing force `M` and
the kinematic reaction `R_old` are geometrically aligned or anti-aligned at that specific joint/trial
(`docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md` §4: acute angles, 45–70°, at the shoulder vs. the
knee/hip's near-antagonist >150° separation) — the sign bug does not push in one direction
universally; each joint's correction direction was measured on its own terms, not assumed from the
others.

## 6. QC finding: one gate's on-disk state now diverges from its origin doc's prose (confirmed benign)

Machine cross-check against `sign_bug_audit_results.json` found `part2.identity.{knee,hip}.pass =
false`, with the algebraic-identity residual at **3737.5 N (knee) / 4338.0 N (hip)** — starkly
different from `docs/MECHANISM_SIGN_BUG_AUDIT.md` §2's own printed "4.5×10⁻¹³ N both — PASS." This is
**not a new, undiagnosed problem**: the JSON file's mtime (15:54:20) postdates the audit doc's
(15:43:55) because `docs/MECHANISM_SIGN_BUG_REMEDIATION.md` explicitly re-ran `audit_sign_bug.py`
after patching the source function, and its own §2 already predicts exactly this: once
`static_opt_knee.knee_crossing_muscles_and_forces` is fixed at the source, the diagnostic's
"published_buggy" code path is no longer actually buggy — so the identity check comparing "buggy"
against "fixed" (`F_fixed == 2·R_old − F_buggy`) necessarily stops holding, since both sides are now
the same corrected answer. This scorecard's contribution is confirming that prediction **with the
exact numbers** (the remediation doc named the *phenomenon*, not this magnitude) — it does not
undermine the corrected 391.10/386.77%BW headline figures, which remain independently confirmed by
the fidelity-reproduction match and the JointReaction cross-check (both unaffected by this specific,
now-stale-by-construction internal diagnostic).

## 7. Cross-cutting honest gaps (inherited by every joint in this family, not repeated per-row above)

- **Static Optimization is an effort-minimizing solution, not measured EMG** — real co-contraction
  can exceed SO's activation-minimizing solution at every joint in this family.
- **Different populations throughout**: every OrthoLoad anchor is an elderly instrumented-implant
  patient cohort; the twin (subject2) is a healthy young(ish) OpenCap participant (78.2 kg, 1.96 m)
  — a plausibility/ballpark comparison, not a per-subject validation, at every joint.
  hip_gen1's subject-code text parsing is additionally coarse (9 apparent codes vs. 2 expected,
  immaterial to the per-trial %BW comparison).
- **Single trial, right-side primary, subject2 `walking1`/`spine_stoop_lift`/`arm_abduction`/
  `elbow_static`** — no claim of generality across subjects, trials, gait speeds, or the left side
  (left-side numbers, where reported, are coarse order-of-magnitude symmetry checks only).
- **Real force-plate GRF** for the walking-based joints (knee/hip/ankle) — the "best case," not the
  degraded kinematics-only fallback; the spine has no captured lifting trial at all (static pose
  substituted, disclosed); the shoulder/elbow have no captured arm-motion trial at all anywhere in
  the twin's 10-subject corpus (synthetic/prescribed kinematics substituted, disclosed).
- **Savitzky-Golay differentiation** (window swept 70–170 ms depending on joint, shown stable in
  every cert) is a design choice, not the literal raw signal.

## 8. Files

- This doc: `docs/MECHANISM_JOINT_FORCE_SCORECARD.md` (new, this session).
- Source docs synthesized (all read in full, not summarized-from-memory):
  `docs/MECHANISM_JOINT_FORCE_VALIDATION.md`, `docs/MECHANISM_STATIC_OPT.md`,
  `docs/MECHANISM_HIP_FORCE.md`, `docs/MECHANISM_ANKLE_FORCE.md`, `docs/MECHANISM_SPINE_FORCE.md`,
  `docs/MECHANISM_SHOULDER_FORCE.md`, `docs/MECHANISM_ARM_MUSCLES.md`,
  `docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md`, `docs/MECHANISM_ELBOW_FORCE.md`,
  `docs/MECHANISM_SIGN_BUG_AUDIT.md`, `docs/MECHANISM_SIGN_BUG_REMEDIATION.md`,
  `docs/MECHANISM_FORCE_SCENES_BATCH.md` (the corpus-SO negative, §3 item 2), `docs/MECHANISM_MSK_BUILD_PLAN.md`
  (anchor-tier cross-reference only, §1 legend).
- Source JSONs machine-cross-checked this session (§1 lists paths); every headline number in the
  master table was pulled directly from these files, not transcribed from doc prose.
- Not touched, not re-run: no scripts under `scripts/msk/` were executed or edited this session;
  this is a read-only synthesis over already-committed results.
