# MECHANISM SIGN BUG REMEDIATION — fix at source, verify, banner (2026-07-21)

Closes out `docs/MECHANISM_SIGN_BUG_AUDIT.md`'s read-only diagnosis (which found, but by explicit
task scope did not fix, a sign bug in `scripts/msk/static_opt_knee.py`'s
`knee_crossing_muscles_and_forces`, live in the knee/hip/ankle certs' own self-computed numbers).
This document: applies the minimal fix at the source, re-verifies it two independent ways (a
known-answer toy test and reproduction of the audit's own corrected numbers), enumerates and checks
the status of every caller in the repo, reconciles the spine cert's already-applied local fix
against the now-fixed shared helper, banners the three affected docs (plus a pointer in a fourth),
and states the honest new headline with a verified citation. Every number below is machine-measured
this session, not recalled. Isolation respected: bodytwin only, no git commit/push (coordinator
commits).

## 1. The one-line fix

`scripts/msk/static_opt_knee.py`, `knee_crossing_muscles_and_forces`, the crossing-index ternary
(now at line 346, shifted by an explanatory comment block):

```python
# before (bug): backwards in BOTH branches
inside_idx, outside_idx = (k, k + 1) if in_k1 else (k + 1, k)

# after (fix): when in_k1 is True, body k+1 -- not k -- is the inside point
inside_idx, outside_idx = (k + 1, k) if in_k1 else (k, k + 1)
```

Geometric derivation (not a guess): `in_k1 = bodies[k+1] in inside_set`. The branch only executes
when `in_k != in_k1` (a crossing). If `in_k1` is `True`, body `k+1` is the one INSIDE the free-body
cut, so `inside_idx` must be `k+1`. The pre-fix code assigned `inside_idx = k` in that exact
branch — backwards. This is the identical branch-swap `validate_spine_force.py`'s own
`crossing_muscles_unit_dirs` (line 238) already applied independently for the same pattern, and
identical to the audit's own local `crossing_muscles_and_forces_FIXED`.

## 2. Toy-test verdict — re-ran `scripts/msk/audit_sign_bug.py` post-fix (not hand-rolled)

Task instruction was to re-run the audit's own toy test to confirm `+ground_truth`. Ran
`.venv-msk/bin/python3 scripts/msk/audit_sign_bug.py` unmodified, post-fix. Two things must be
disentangled in its output — the raw per-case DATA (the actual confirmation) vs. this script's own
internal PASS/FAIL booleans (several of which were pre-registered to detect "is the bug still
present," so now correctly flip to FAIL, which is the desired outcome, not a defect):

| toy case (4 diverse geometries + the audit's own reference fix) | verdict | matches +truth | matches −truth |
|---|---|:---:|:---:|
| A: inside point first in path order | **CORRECT** | True | False |
| B: inside point second in path order | **CORRECT** | True | False |
| C: crossing mid-path (3-point path) | **CORRECT** | True | False |
| D: arbitrary 3-D axis orientation | **CORRECT** | True | False |
| E: audit's independent reference-fixed copy (sanity) | **CORRECT** | True | False |

**Confirmed exactly as instructed: `np.allclose(computed, +ground_truth)` is True and
`np.allclose(computed, -ground_truth)` is False, for all 4 diverse real-function cases.** (Full
JSON: `data/msk_smoketest/subject2_walking1/sign_bug_audit/sign_bug_audit_results.json`.)

The script's own `part1_gate_pass`/`fidelity_gate`/`identity_gate` booleans now read False (exit
code 1) — expected and diagnosed, not a red flag: those gates were pre-registered to check "does the
real function still reproduce the OLD buggy behavior / the OLD hardcoded published constants
(233.20/235.29 %BW)." Post-fix it correctly does not — e.g. `FIDELITY [knee]: published=233.20106
audit_rerun=391.09952 rel_diff=67.7%` — the "FAIL" here is the fix working, read backwards. The
`audit_sign_bug.py` script itself was NOT edited (it is a historical diagnostic artifact whose own
booleans describe a specific past state); its output is reinterpreted here, not altered.

## 3. Corrected numbers — reproduces the audit's 391/387 %BW, plus a fresh ankle measurement

Part 2 of the audit (which reuses the cached SO/JointReaction `.sto` output on disk — no Ipopt
re-run) recomputes knee+hip through the now-fixed `sok.knee_crossing_muscles_and_forces` (labeled
`published_buggy` in that script — a misnomer post-fix, kept because the script itself is
unmodified) and cross-checks it against the audit's own independently-coded local fixed copy
(`audit_fixed`):

| joint | published (bug artifact) | **corrected** | independent JointReaction | corrected-vs-JR agreement | corrected-vs-audit's-own-fixed-copy |
|---|---:|---:|---:|---:|---:|
| knee | 233.20 %BW (ratio 0.903) @ t=0.51s | **391.0995 %BW (ratio 1.515)** @ t=0.51s | 391.1148 %BW | **0.0039%** | 0.000000 %BW diff (bit-identical, 0.0° vector angle) |
| hip | 235.30 %BW (ratio 0.859) @ t=0.67s | **386.7679 %BW (ratio 1.412)** @ t=0.55s | 387.0372 %BW | **0.070%** | 0.000000 %BW diff (bit-identical, 0.0° vector angle) |

Ankle is outside `audit_sign_bug.py`'s scope (knee+hip only), so I re-ran `validate_ankle_force.py`
directly instead (also cheap — reuses the same cached sibling SO/JR output, no Ipopt re-run), which
independently confirms the task brief's approximate "≈484/~1.02" precisely:

| joint | published (bug artifact) | **corrected** (direct re-run, exit reflects an unrelated stale-constant issue, §6) | independent JointReaction | agreement |
|---|---:|---:|---:|---:|
| ankle | 286.68 %BW (ratio 0.601) | **484.4865 %BW (ratio 1.016)** @ t=0.58s | 484.4798 %BW | **0.0014%** |

I also directly re-ran `validate_hip_force.py` (not just the audit's parallel reimplementation) as a
second, independent cross-check: exit 0, all its own internal gates PASS, self-computed =
386.76785 %BW — identical to the audit's number to 5 decimal places. **Three joints, two
independent code paths each (the live cert script itself + the audit's standalone reimplementation
for knee/hip), all converge on the same corrected numbers, all matching their respective
bug-immune JointReaction cross-check to <0.1%.** JointReaction never calls
`knee_crossing_muscles_and_forces` at any joint (`docs/MECHANISM_SIGN_BUG_AUDIT.md` §2) — this is an
external, non-tautological anchor, not the fix confirming itself.

Sanity/scope check on the bug's exact boundary (explains why some numbers move and others don't): a
uniform sign flip applied to every crossing muscle negates the entire summed muscle-force VECTOR,
but `np.linalg.norm(v) == np.linalg.norm(-v)` — so any number that is itself a norm of a
muscle-crossing-only sum (not the final `R − Σmuscle` subtraction) is structurally blind to this
bug. Confirmed by direct measurement, not assumed: the ankle doc's Achilles/triceps-surae numbers
(320.2649 %BW scalar sum, 318.83 %BW vector-magnitude parallelism check) are bit-identical before
and after the fix. Only the final `F_bone_contact = R − Σmuscle-crossing` subtraction (where the
sign of the subtracted term matters, since `R − M ≠ R − (−M)` in general) is affected.

## 4. Every caller — checked, not assumed

`grep -rl "knee_crossing_muscles_and_forces" scripts/msk/*.py` found 9 files. Status of each:

| file | how it uses the function | effect of the source fix | status |
|---|---|---|---|
| `static_opt_knee.py` | defines it | fixed at source (§1) | **FIXED** |
| `validate_hip_force.py` (line ~325) | calls `sok.knee_crossing_muscles_and_forces` directly, unmodified | auto-fixed | **FIXED, confirmed by direct re-run** (§3) |
| `validate_ankle_force.py` (lines ~249, ~251) | calls `so_knee.knee_crossing_muscles_and_forces` directly for BOTH the ankle-r cut and the knee-r reproducibility cut | auto-fixed | **FIXED, confirmed by direct re-run** (§3); see §6 for a stale-constant side effect in this file's own internal gate |
| `validate_spine_force.py` | defines its OWN independent local copy, `crossing_muscles_unit_dirs` (line 201-246) — does **not** call the shared function anywhere (grep-confirmed: only comments reference the name) | none — was never wired to the shared function | **UNAFFECTED, NOT double-fixed** (§5) |
| `audit_sign_bug.py` | calls `sok.knee_crossing_muscles_and_forces` (now-fixed) for its `published_buggy` variant + its own independent local `crossing_muscles_and_forces_FIXED` for comparison | diagnostic tool itself; re-run and reinterpreted, not edited | **RE-VERIFIED** (§2-3) |
| `validate_elbow_force.py` (line 157) | calls `sok.knee_crossing_muscles_and_forces` for the TWIN's own elbow — returns an EMPTY contribution dict (0 muscles cross the elbow structurally; `f_bone_contact_equals_pure_reaction_identically: True`) | **no numerical effect on the twin's own elbow number** — a sign bug on an empty dict changes nothing | **UNAFFECTED** (twin); see §7 for a separate, independently-defined local copy in the same file used for an illustrative reference model |
| `force_scenes_batch.py` (line 488) | calls `so_knee.knee_crossing_muscles_and_forces` directly for a 3-clip ankle/knee/hip movement-type panel | auto-fixed for any future run | **FIXED going forward; on-disk doc/JSON still stale** — see §7 (out of this task's 4-doc scope, flagged not fixed) |
| `reconcile_diff_scheme.py` (lines 328, 382) | calls `sok.knee_crossing_muscles_and_forces` identically in BOTH its "Pass A" and "Pass B" comparison variants | auto-fixed for any future run | **FIXED going forward; on-disk JSON still stale** — see §7 (largely common-mode in this script's own A-vs-B comparison, so less distorted than the headline certs were; no dedicated doc exists) |
| `validate_shoulder_force_with_muscles.py` | **does not call the shared function at all** — comments only; defines an independent local copy `arm_crossing_muscles_and_forces` (lines 244-279) that duplicates the IDENTICAL pre-fix buggy ternary | **NOT auto-fixed by the source patch** (separate code, separate bug instance) | **STILL LIVE — flagged, not fixed, out of scope** (§7, most significant additional finding) |

## 5. Spine reconciliation — confirmed no double-fix

`validate_spine_force.py`'s `crossing_muscles_unit_dirs` (lines 201-246) is a **completely
independent function**, not a call site of `static_opt_knee.knee_crossing_muscles_and_forces` — it
has its own copy of the crossing-detection loop, with its own already-correct branch order (line
238: `inside_idx, outside_idx = (k + 1, k) if in_k1 else (k, k + 1)` — already matching the fix, and
the function's own inline comment already explains this exact reasoning, dated to the spine agent's
original session). Grepped this session to be sure: zero calls to
`sok.knee_crossing_muscles_and_forces` or bare `knee_crossing_muscles_and_forces` anywhere in
`validate_spine_force.py` — the two hits are both comments referencing the name for context. **The
shared-helper fix in `static_opt_knee.py` has zero effect on the spine cert — there is no double-fix
scenario anywhere in this codebase.** The spine cert's own published numbers (296.59/389.17 %BW
Tier-2, `docs/MECHANISM_SPINE_FORCE.md`) are untouched by this remediation and were not re-derived
here (out of scope, consistent with the original audit's own scope boundary).

## 6. A predictable side effect: two hardcoded "prior published" constants are now stale

Two scripts compare their own live re-derivation against a **hardcoded historical constant**
representing the OLD buggy knee number, as an internal self-consistency/reproducibility gate. Both
now correctly report a mismatch when re-run — not a new bug, a direct, mechanical consequence of the
fix working:

- `audit_sign_bug.py`: `PUBLISHED_KNEE_SELF_COMPUTED_PCT_BW = 233.20106467077753` /
  `PUBLISHED_HIP_SELF_COMPUTED_PCT_BW = 235.2956788485241` (lines 64-65) — its own `fidelity` gate now
  reads FAIL (§2), by design (this script is a historical diagnostic; not edited).
- `validate_ankle_force.py`: `PRIOR_KNEE_SELF_COMPUTED_PCT_BW = 233.20` (line 156) — its own §8 "FREE
  reproducibility gate" (knee-r cut vs already-published knee cert) now reports the self-computed
  candidate MISMATCHING at 67.710% relative diff when re-run (measured directly, §3), which flips
  `repro_pass`/`pipeline_valid`/the script's overall exit code to FAIL even though the ankle numbers
  themselves are now more correct than before. **Flagged, not fixed here** (out of this task's
  explicit scope of "fix `static_opt_knee.py`, banner the 4 named docs"): a follow-up should update
  this constant to `391.10` so `validate_ankle_force.py` exits clean again.

## 7. Additional findings discovered during the caller audit — flagged, NOT fixed (out of scope)

Symmetric QC requires surfacing these with the same evidence standard as the confirmed fixes above,
even though acting on them would exceed this task's explicit 1-source-file + 4-doc scope (the same
scope discipline `docs/MECHANISM_SIGN_BUG_AUDIT.md` itself used: found, disclosed, not silently
fixed). None of these were edited.

1. **`validate_shoulder_force_with_muscles.py` — a SEPARATE, independently-defined instance of the
   identical bug, most significant of these findings.** Its own `arm_crossing_muscles_and_forces`
   (lines 244-279) was written standalone (the grafted-arm-muscle model's muscles carry no `_r`/`_l`
   suffix, so the shared function's side-filter would have silently matched nothing) and copy-pasted
   the pre-fix ternary verbatim: `inside_idx, outside_idx = (k, k + 1) if in_k1 else (k + 1, k)`
   (line 272) — NOT auto-fixed by this remediation's source patch, because it never calls the shared
   function. Confirmed LIVE, not vacuous: the on-disk
   `data/msk_smoketest/subject2_arm_abduction/shoulder_force_with_muscles/shoulder_force_with_muscles_results.json`
   shows `n_crossing_muscles: 22` (a real, non-empty crossing set) and a headline
   `f_new_peak_pctBW: 24.177 %BW` (vs pre-graft 5.306 %BW, OrthoLoad glenohumeral anchor 72.2 %BW) —
   this number is likely itself a bug artifact, by the same mechanism as the knee/hip/ankle certs.
   No dedicated `docs/MECHANISM_*.md` file exists for this result (grepped repo-wide — none found; it
   lives only in the module docstring + JSON) — **recommend a dedicated follow-up remediation**
   (fix the branch swap in this file's own local copy, re-run, and either write or update whatever
   doc references the 72.2 %BW shoulder-with-muscles comparison).
2. **`validate_elbow_force.py`'s illustrative `arm26` reference-model comparison** (`crossing_muscles_generic`,
   lines 374-407) carries the same copy-pasted pre-fix ternary (line 400) — but this is explicitly a
   non-twin, non-subject-scaled illustrative demonstration (Holzbaur/Murray/Delp 2005's 6-muscle
   `arm26.osim`), not the twin's own number (which is separately, structurally unaffected — 0 muscles
   cross the twin's own elbow). `docs/MECHANISM_ELBOW_FORCE.md` line 27 cites this arm26 number
   ("455.2 N... 59.4 %BW") as illustrative context — low priority, but flagged for a follow-up
   spot-check since it is technically wrong under the same mechanism.
3. **`force_scenes_batch.py`** (§4) auto-inherits the source fix for any future run, but its
   dedicated doc `docs/MECHANISM_FORCE_SCENES_BATCH.md` (confirmed via grep to reference
   `knee_crossing_muscles_and_forces` directly, line 178, and to publish concrete pre-fix %BW numbers
   at lines 102-105 and 183-195, e.g. "knee 590-617%BW... hip 770-1029%BW... far outside any
   published literature") is stale and outside this task's 4-named-doc scope — **recommend the same
   banner-and-recompute treatment** as this remediation applied to the knee/hip/ankle docs.
4. **`reconcile_diff_scheme.py`** (§4) also auto-inherits the fix, but its bug exposure is
   structurally different: both of its comparison passes ("Pass A" raw-q and "Pass B" JR-consistent
   filtered-q") call the identical pre-fix function identically, so the sign error was largely
   common-mode WITHIN this script's own A-vs-B comparison (its actual subject, the
   differentiation-scheme gap) rather than a directional distortion of one method vs. the other. Its
   own on-disk result (`data/msk_smoketest/subject2_walking1/diff_scheme_reconcile/
   reconcile_diff_scheme_results.json`) already reports `overall_verdict: MIXED/PARTIAL` and a
   NEGATIVE hip `closure_frac_accel_only` (−0.207) — not a clean "differentiation scheme explains the
   gap" result even before this remediation — consistent with, not contradicted by,
   `docs/MECHANISM_SIGN_BUG_AUDIT.md` §4's finding that the differentiation-scheme hypothesis was
   never the dominant mechanism. No dedicated doc exists for this script. Lowest priority of the four
   findings in this section.

## 8. Honest new headline, with a verified citation

The corrected knee/hip self-computed numbers now match their own bug-immune `opensim.JointReaction`
cross-check to <0.1% at both joints (§3) — this convergence, from two structurally independent code
paths, is the strongest evidence in this remediation that the CORRECTION, not the original, is
physically right. **The honest new headline is not "the twin's muscle-driven estimate lands just
under the in-vivo anchor" (the old, bug-artifact framing) but "the twin OVER-predicts in-vivo
knee/hip contact force by ~1.4-1.5× (knee 1.515×, hip 1.412×)."**

This is not an idiosyncratic failure of this twin's pipeline — it is a **documented tendency of
exactly this class of method**. Verified live this session (NCBI eutils `esearch`/`efetch`, since
`WebSearch` was unavailable — session-wide quota exhausted, the same pre-existing constraint already
flagged in `docs/MECHANISM_MUSCLE_AUDIT.md`/`docs/MECHANISM_ANKLE_FORCE.md`):

> Moissenet F, Chèze L, Dumas R. "A 3D lower limb musculoskeletal model for simultaneous estimation
> of musculo-tendon, joint contact, ligament and bone forces during gait." *J Biomech.*
> 2014;47(1):50-58. DOI [10.1016/j.jbiomech.2013.10.015](https://doi.org/10.1016/j.jbiomech.2013.10.015),
> PMID [24210475](https://pubmed.ncbi.nlm.nih.gov/24210475/). Verbatim (machine-fetched abstract,
> not paraphrased first, quoted precisely second): *"Musculo-tendon forces and joint reaction forces
> are typically estimated using a two-step method, computing first the musculo-tendon forces by a
> static optimization procedure and then deducing the joint reaction forces from the force
> equilibrium. However, this method does not allow studying the interactions between musculo-tendon
> forces and joint reaction forces in establishing this equilibrium and the joint reaction forces are
> usually overestimated."*

This twin's own method is exactly this "two-step" pattern: the Newton's-law reaction `R` is computed
first (from kinematics+GRF alone), Static Optimization's muscle tensions are computed separately,
and the two are combined post-hoc (`F_bone_contact = R − Σmuscle-crossing`) rather than solved
simultaneously — precisely the architecture Moissenet et al. identify as prone to overestimating
joint reaction/contact forces. **The ankle is the exception, not a counter-example**: its anchor
(§0 of `docs/MECHANISM_ANKLE_FORCE.md`) is a weaker-tier published cadaveric/model-literature band,
not a direct in-vivo instrumented measurement, and the corrected number lands almost exactly on it
(ratio 1.016) rather than over-predicting it — there is no over-prediction-vs-in-vivo story at the
ankle because there is no in-vivo anchor there to over-predict against.

## 9. Files

- `scripts/msk/static_opt_knee.py` — the one-line fix (§1), with an inline comment pointing to this
  doc and the audit.
- `docs/MECHANISM_STATIC_OPT.md`, `docs/MECHANISM_HIP_FORCE.md`, `docs/MECHANISM_ANKLE_FORCE.md` — top
  correction banners + light inline pointers at the specific superseded sections (§5/§6 knee, §5/§6
  hip, §2.6/§5/§6 ankle); original prose preserved verbatim below each banner, not silently rewritten.
- `docs/MECHANISM_JOINT_FORCE_VALIDATION.md` — a short pointer banner (its own 102.17 %BW number is
  unaffected; it never calls the crossing-muscle function).
- `data/msk_smoketest/subject2_walking1/sign_bug_audit/sign_bug_audit_results.json` — regenerated by
  the post-fix audit re-run (§2-3).
- `data/msk_smoketest/subject2_walking1/hip_force_validation/hip_force_validation_results.json` and
  `data/msk_smoketest/subject2_walking1/ankle_force_validation/ankle_force_validation_results.json` —
  regenerated by direct post-fix re-runs of the actual live cert scripts (§3), replacing the stale
  buggy numbers on disk with the corrected ones.
- **Not touched** (explicit scope boundary, §5, §7): `validate_spine_force.py`,
  `validate_shoulder_force_with_muscles.py`, `validate_elbow_force.py`, `force_scenes_batch.py`,
  `reconcile_diff_scheme.py`, `docs/MECHANISM_SPINE_FORCE.md`, `docs/MECHANISM_FORCE_SCENES_BATCH.md`,
  `docs/MECHANISM_ELBOW_FORCE.md`, `audit_sign_bug.py`, and
  `data/msk_smoketest/subject2_walking1/static_optimization/static_opt_knee_results.json` (still
  holds the pre-fix 233.20 %BW number on disk — re-generating it requires a fresh, expensive Ipopt
  Static-Optimization solve via `static_opt_knee.py`'s `main()`, which unconditionally re-runs SO
  [no cache-reuse guard, unlike the hip/ankle scripts] and would also regenerate the shared SO/JR
  cache every other script in this family reads; not done here — the audit's cached-reuse recompute
  plus the direct hip/ankle re-runs already give three independently cross-validated, <0.1%-from-JR
  corrected numbers without that cost or risk).

## 10. Honest gaps (full list)

1. **`static_opt_knee_results.json` and `sign_bug_audit_results.json`'s own historical copy were not
   regenerated via a fresh Ipopt run** — the corrected knee number (391.0995/391.10 %BW) is
   triple-confirmed (audit cached-reuse recompute, bit-identical cross-check against the audit's
   independent fixed copy, and exact match to the already-published JointReaction number) but was
   never re-derived from a brand-new Static Optimization solve in this session; only hip and ankle
   were end-to-end re-run through their actual live scripts.
2. **Four additional call sites (§7) carry a live or illustrative instance of this bug and were
   deliberately left unfixed**, per this task's explicit scope (fix `static_opt_knee.py` + banner 4
   named docs) — most significant is `validate_shoulder_force_with_muscles.py`, whose own headline
   number (24.18 %BW vs a 72.2 %BW OrthoLoad anchor) is likely a live bug artifact with no
   correcting doc anywhere in the repo.
3. **Two hardcoded "prior published" constants are now stale** (§6) — `audit_sign_bug.py` lines
   64-65 and `validate_ankle_force.py` line 156 — flagged, not updated (the first is a historical
   diagnostic script's own record of the state it audited; the second causes a cosmetic exit-code
   FAIL on any future re-run of `validate_ankle_force.py` despite its actual numbers being correct).
4. **This remediation does not re-verify OrthoLoad anchors, model files, or SO's underlying
   activation-minimization limitation** — all inherited unchanged from the original certs and audit;
   out of scope (this remediation is about the sign convention of one downstream geometry function,
   nothing upstream or downstream of it).
5. **Single trial, right side, subject2 `walking1`** — same scope caveat as every cert in this
   family; not re-verified across subjects/trials/left leg.
