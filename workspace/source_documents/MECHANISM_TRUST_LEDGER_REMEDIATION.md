# MECHANISM TRUST LEDGER REMEDIATION — disposition of the §10 high-scrutiny flags (2026-07-21)

Fixes, at source, the specific rows `docs/MECHANISM_TRUST_LEDGER.md` §10 flagged as needing the
highest scrutiny. Method per flag: re-derive directly from the pipeline's own JSON/code (never
trust the ledger's or the original doc's prose alone), then either correct the label, persist the
missing evidence, banner the doc honestly, or — where a flag itself turns out to be wrong — refute
it with evidence. Nothing was deleted; every original number/claim is preserved verbatim alongside
its correction. Isolation respected throughout: bodytwin only, no git commit, no git push (files
listed at the end are on disk, ready for the coordinator to stage).

## Disposition summary

| # | flag | disposition |
|---|---|---|
| 1 | TRUNK_FLEXORS.md "gravity alone" mislabel | **CORRECTED** |
| 2 | 8-point robustness sweep unpersisted | **PERSISTED** (fresh re-run, both docs bannered) |
| 3 | SCAPULA_CLAVICLE.md 2:1 rhythm circularity | **ADJUDICATED: ledger's flag is CORRECT** (largely circular), with a defended, narrower exception |
| 4 | STATIC_OPT.md stale 233%BW JSON | **CORRECTED (stale-marked in-file)** |
| 5a | CLIMBING_SCENE gate table unpersisted | **PERSISTED** |
| 5b | YOUTUBE_CASCADE_POC "23/23" unpersisted | **BANNERED** (doc-prose-only, disclosed) |
| 5c | MSK_ELASTIC_BAND neutral-pose sub-table unpersisted | **PERSISTED** |
| 6 | STATE.md "~92/20" unreconcilable | **BANNERED** (stale/superseded) |

---

## 1. TRUNK_FLEXORS.md — "−41.749 N·m required by gravity alone" — CORRECTED

**Verified against the raw JSON** (`data/msk_smoketest/band_static_opt/band_static_opt_results.json`
→ `gravity_decomposition_lumbar_extension`, re-read directly this session):

| quantity | value (N·m) |
|---|---:|
| gravity alone (no GRF, no band) | **−50.10008001348936** |
| gravity + GRF | −50.10008001316368 (GRF ≈ 0 here) |
| gravity + GRF + band ("net of band") | **−41.74935672599269** |
| band's marginal contribution | +8.350723287170986 |

The ledger's flag is correct: −41.749 N·m is gravity **net of** the band's own +8.351 N·m
contribution, not gravity alone — gravity alone is the larger-magnitude number, −50.100 N·m. The
underlying computation was always right; only the label was wrong. **Fix**: added a `⚠ CORRECTION`
banner at the top of `docs/MECHANISM_TRUNK_FLEXORS.md` stating both values in a table, and corrected
the inline Bottom-line sentence that originally called −41.749 "required by gravity alone." No other
finding in that doc (moment-arm signs, capacity ceiling, SO recruitment, the required-force
reproduction gate) is affected — they never depended on this label.

## 2. The "8-point robustness sweep" — PERSISTED (fresh re-run, not a copy of memory)

Confirmed directly by reading `scripts/msk/band_static_opt.py`: `gravity_only_decomposition()` is
called **exactly once**, at the single published pose — no array of 8 (or any N) points existed
anywhere on disk. Fixed by writing `scripts/msk/band_static_opt_lumbar_sweep.py`, which imports (does
not duplicate) `band_static_opt.build_with_band_model()`/`gravity_only_decomposition()` and
`band_posterior_chain.set_pose()`/`generalized_force_multi()`, sweeps `lumbar_extension` over 8
evenly-spaced points across the doc-cited [−20°,+30°] range holding the rest of the hip-hinge pose
fixed, and persists the result to
`data/msk_smoketest/band_static_opt/lumbar_extension_robustness_sweep.json`.

**Result — the qualitative claim holds, with one precise correction to which quantity it describes**:

| lumbar_extension (deg) | required tau (N·m, gravity+GRF+band) | gravity_only (N·m) |
|---:|---:|---:|
| −20.000 | −3.184 | −15.167 |
| −12.857 | −14.477 | −25.745 |
| −5.714 | −25.642 | −35.924 |
| +1.429 | −36.496 | −45.546 |
| +8.571 | −46.858 | −54.459 |
| +15.714 | −56.559 | −62.529 |
| +22.857 | −65.431 | −69.627 |
| +30.000 | −73.322 | −75.645 |

Both columns are negative throughout, monotonic, no sign flip — **confirms** the "robust across an
8-point sweep" qualifier both docs cite. But the doc-cited range **"−15.2 → −75.6 N·m" matches the
`gravity_only` column** (−15.167 → −75.645, matching to stated rounding), **not** the `required tau`
column (−3.184 → −73.322) that the prose attributes it to — a real, disclosed labeling imprecision
in the original (never-persisted) claim, not a numerical error. Over-determination check: re-evaluated
at exactly lumbar_extension=5° with an independently-rebuilt model (a fresh `build_with_band_model()`
call, not a reused in-memory number), required tau = −41.749149 N·m (0.0005% from the published
−41.749357) and gravity_only = −50.099872 N·m (0.0004% from −50.100080) — an external cross-check
against a number already on disk, not a tautology, since the model was rebuilt from scratch.

**Fix**: both `docs/MECHANISM_BAND_STATIC_OPT.md` (Sec.5) and `docs/MECHANISM_TRUNK_FLEXORS.md`
(top banner) now point to the persisted JSON and disclose the gravity_only-vs-required_tau
labeling correction.

## 3. SCAPULA_CLAVICLE.md — 2:1 rhythm circularity — ADJUDICATED: the ledger's flag is CORRECT

**This is the flag the task marked most important; adjudicated against the actual code
(`scripts/msk/add_scapula_clavicle.py`), not just doc prose, on both sides.**

Facts established by reading the code directly:
- `RHYTHM_RATIO_SCAPULAR_TO_GH = 0.5` (line 113) is hard-coded as the `CoordinateCouplerConstraint`'s
  `LinearFunction` slope (line 818) — `scapula_upward_rot_r` is **forced, exactly, by construction**
  to equal `slope * arm_add_r` at every assembled pose.
- In `verify_and_sweep()` (line 830), the other 3 scapular DOFs (`scapula_abduction_r`,
  `scapula_elevation_r`, `scapula_winging_r`) are set once to fixed defaults and never driven during
  the sweep — confirmed by reading the coordinate-setup code (`set_coord(...)` calls, lines 533-539)
  and the sweep loop itself (only `arm_add_r`/`arm_flex_r`/`arm_rot_r` are touched).
- Consequently `total_humerothoracic_elev_deg_FK` (the "independent" quantity, computed via a
  humerus-long-axis FK measurement, lines 868-870) is, in this sweep, a **fully deterministic
  function of `arm_add_r` alone** — parameterized entirely by (a) the same hard-coded 0.5 and (b)
  this model's own fixed, unvalidated joint-axis/ellipsoid geometry. No real scapular
  motion-capture data enters this computation anywhere (the doc's own §7.5/§8 already say none
  exists in this corpus).

**Verdict**: the headline claim ("scapulohumeral rhythm measured 0.489, matches target 0.500") is
**largely circular**, exactly as the ledger flagged. A decisive test of circularity: if
`RHYTHM_RATIO_SCAPULAR_TO_GH` had been set to 0.3 instead of 0.5, this identical sweep+FK
methodology would report "ratio measured ≈0.29" with a similarly small residual — the check cannot
discriminate 0.5 from any other assumed constant, so it supplies ~no independent evidence for
0.489-vs-0.500 specifically.

**Symmetric fairness — what the doc gets right and the ledger doesn't fully credit**: the ~1.06°
max discrepancy between the FK-derived total elevation and the coupler's own internal state is
**not** a pure tautology in the "compute the same number twice" sense. `total_elev_deg` is computed
via a genuinely different mathematical operation (full rigid-body FK through a real 2-joint chain)
than reading the coupler's state variable directly. If the GH-joint axis and the scapula's
upward-rotation axis were exactly parallel/collinear in this model's geometry, 3D-rotation
composition would force these two numbers to agree to floating-point precision — a pure identity,
zero information. They do **not** agree to floating-point precision (1.06° residual, from real
axis non-collinearity + the ellipsoid joint's translation effects) — proof this is a genuinely
separate code path that could have diverged far more, and, before this session's own sign-bug fix
(§6 of the doc), **did** diverge by ~150-160°. This is precisely the mechanism that caught that
real bug. So the FK cross-check is legitimate and has demonstrated teeth **as a model-wiring/
geometry-composition self-consistency device** (differential testing against a second code path) —
it is simply not evidence that 0.489/0.500 reflects true human anatomy, which is a different
question this corpus has no data to answer.

**Conclusion**: the ledger's flag stands — not too harsh. **Fix**: banner added to the top of
`docs/MECHANISM_SCAPULA_CLAVICLE.md` separating "Claim A" (circular, the ratio-validation framing)
from "Claim B" (real, the wiring-QC value), with the code line numbers and reasoning above. Original
headline table and §6 prose preserved verbatim below the banner.

## 4. STATIC_OPT.md — stale 233%BW JSON — CORRECTED (stale-marked in-file)

Confirmed: `data/msk_smoketest/subject2_walking1/static_optimization/static_opt_knee_results.json`
still holds `self_computed.bone_contact_residual_peak_pct_bw = 233.20106467077753` (the pre-sign-fix
value) — a naive automated check reading only this file gets the wrong, superseded number, exactly
as flagged. **Fix (does not silently overwrite history)**: added a new top-level
`_STALE_PRE_SIGN_FIX_NOTICE` key to that same JSON file, stating the stale value, the corrected
value (391.10%BW), the root cause, where it was fixed, and why the file wasn't regenerated (doing so
requires an expensive, unconditional fresh Ipopt Static-Optimization solve that would also
regenerate a shared SO/JR cache other scripts read — the corrected number is instead
triple-confirmed in `docs/MECHANISM_SIGN_BUG_REMEDIATION.md` §3 via cached-reuse recompute + an
independent fixed-copy implementation + <0.1% agreement with this same file's own never-buggy
`jr_summary.peak_pct_bw` = 391.11%BW). Every original key/value in the file is untouched. Also
updated `docs/MECHANISM_STATIC_OPT.md`'s own "Files" section to explicitly warn the JSON is stale and
point to the corrected number.

## 5a. CLIMBING_SCENE.md gate table — PERSISTED

`docs/MECHANISM_CLIMBING_SCENE.md` Sec.5's "~3% of world record" anchor-demo table had no persisted
array anywhere (confirmed absent twice in the prior session, both times reproduced only by
live-re-executing the script). Fixed by writing `scripts/msk/climbing_scene_persist_anchor_table.py`,
which imports (does not duplicate) `climbing_scene.anchor_demo()` and calls it for the 3
`gt_solo_run` clips, persisting to `data/msk_pose/climbing/anchor_demo_gate_table.json`.

Fresh re-run reproduces the doc's table exactly: climb-only velocities 3.322 / 3.387 / 3.375 m/s
(clip 1 tuned, clips 2-3 stock, matching the doc's own per-clip source selection), 2/3 clips full
gate pass, max deviation from the IFSC men's-WR-implied pace (15m/4.54s = 3.304 m/s) = **2.51%** —
within the doc's claimed "~3%." Banner added to the doc pointing to the new JSON.

## 5b. YOUTUBE_CASCADE_POC.md "23/23 fetches failed" — BANNERED (doc-prose-only, disclosed)

Confirmed: the qualitative claim ("zero clips downloaded") is independently, trivially re-verifiable
right now via `find data/youtube -type f`. The specific "23/23" count and the list of which 23 video
IDs were tried has no persisted array anywhere. **Not re-run this session** (a deliberate, disclosed
choice, not an oversight): reproducing it means firing ~23+ fresh network requests at YouTube whose
outcome depends on YouTube's current server-side gating, not on anything this repo controls — a
live re-run could legitimately return a different count without that changing anything about what
was true when the doc was written, and would not be a "verification" of a historical claim so much
as a new, separate measurement. Marking prose-only is the lean, honest choice; banner added stating
this explicitly and explaining why re-running wasn't attempted.

## 5c. MSK_ELASTIC_BAND.md neutral-pose ID-bug sub-table — PERSISTED

`docs/MECHANISM_MSK_ELASTIC_BAND.md` Sec.4's "neutral standing pose, gravity only" table
(ankle/knee/hip/lumbar verified-dPE/dq vs ID-reported, showing the ~1000-1800x knee/hip corruption)
had no persisted artifact — only the separate hip-hinge-pose propagation table (Sec.5) was
reproducible from a raw `.sto` file. Fixed by writing
`scripts/msk/elastic_band_neutral_pose_idbug_check.py`, which imports (does not duplicate)
`attach_band.py`'s `MODEL_PATH`/`generalized_force`/`total_effective_pe`/`read_sto_row`/`find_col`,
builds the model at its own **file-default pose** (confirmed live: all rotational coordinates = 0.0,
matching the doc's "all coordinates 0" description), with **no GRF, no band, no external loads at
all** — exactly the doc's own stated falsifier — and runs a real `InverseDynamicsTool` pass,
persisting to `data/msk_smoketest/elastic_band/neutral_pose_idbug_results.json`.

**Fresh reproduction matches the doc almost to the reported digit**:

| coordinate | verified dPE/dq | ID reported | ratio |
|---|---:|---:|---:|
| ankle_angle_r | 1.2229 | 1.17025 | matches |
| knee_angle_r | −0.6858 | **−1252.616** | 1826× |
| hip_flexion_r | 1.3133 | **1337.124** | 1018× |
| lumbar_extension | −7.5911 | −7.61169 | matches |

(doc's own table: −0.69/−1252.62/"1815×", 1.31/1337.12/"1021×" — matches to stated rounding.) As a
bonus over-determination the original table didn't report, the **left side** independently
reproduces the same corruption: knee_angle_l −0.8024 vs ID −1238.532 (1544×), hip_flexion_l 1.2973
vs ID 1340.171 (1033×). Banner added to the doc pointing to the new JSON.

## 6. STATE.md "~92/20" vs the live graph — BANNERED (stale/superseded)

Independently re-verified this session (not just trusted from the ledger): `data/MECHANISM_ANCHOR_GRAPH.json`
currently has **998 total nodes**, status `{OPEN: 923, ASSUMED: 66, REFUTED: 5, DEFERRED: 3, PROVEN: 1}`
— its own `_meta.note` states "status OPEN = designed, NOT measured." This does not reconcile with
`docs/MECHANISM_STATE.md`'s "~92 designed cells are EXECUTED / 20 re-verified," and `bt_memory/LEDGER.jsonl`
is an event stream (not a per-cell tally), so it cannot arbitrate either. Most likely explanation
(not chased further, per lean scope): `data/MECHANISM_ANCHOR_GRAPH.json` is a **different, broader**
whole-body design catalog than whatever STATE.md's "~92" was originally counting (confirmed
different in kind by the ledger's own §11) — i.e. probably two different tracking schemes that were
never the same graph, not a regression in either one. **Fix**: banner added to the top of
`docs/MECHANISM_STATE.md`, dated 2026-07-18 (the oldest doc in this ledger family), marking the
"~92/20" claim stale/likely-superseded rather than asserting it false, with the freshly re-verified
998-node breakdown quoted directly. Rest of the doc preserved verbatim.

---

## Files changed / added this session

**Docs corrected/bannered** (originals preserved, banners/inline fixes only, nothing deleted):
- `docs/MECHANISM_TRUNK_FLEXORS.md`
- `docs/MECHANISM_BAND_STATIC_OPT.md`
- `docs/MECHANISM_SCAPULA_CLAVICLE.md`
- `docs/MECHANISM_STATIC_OPT.md`
- `docs/MECHANISM_CLIMBING_SCENE.md`
- `docs/MECHANISM_YOUTUBE_CASCADE_POC.md`
- `docs/MECHANISM_MSK_ELASTIC_BAND.md`
- `docs/MECHANISM_STATE.md`

**Evidence JSON edited in place** (additive top-level key only, all original values untouched):
- `data/msk_smoketest/subject2_walking1/static_optimization/static_opt_knee_results.json`

**New scripts** (each imports/reuses existing machinery, no physics/analysis re-implemented):
- `scripts/msk/band_static_opt_lumbar_sweep.py`
- `scripts/msk/climbing_scene_persist_anchor_table.py`
- `scripts/msk/elastic_band_neutral_pose_idbug_check.py`

**New evidence JSON (persisted, previously-missing arrays)**:
- `data/msk_smoketest/band_static_opt/lumbar_extension_robustness_sweep.json`
- `data/msk_pose/climbing/anchor_demo_gate_table.json`
- `data/msk_smoketest/elastic_band/neutral_pose_idbug_results.json`

(all three `data/` JSON paths are gitignored, matching this repo's existing convention for every
sibling evidence file in this cert family — e.g. `band_static_opt_results.json` /
`static_opt_knee_results.json` are gitignored the same way — not a new gap.)

**Not touched, out of this task's scope**: `docs/MECHANISM_CORPUS_IK.md`, `scripts/msk/pose_to_opensim_ik.py`,
`opensim.log`, and untracked `scripts/msk/{add_trunk_flexors_evidence.json,bone_stress.py,emg_driven.py,
rra_task_gain_tune.py,wobbling_mass.py}` showed as modified/untracked in `git status` at the start of this
session — these are a **concurrent mechanism instance's own work** (this repo is explicitly shared between
mechanism1/mechanism2 per `COORDINATOR.md` §0), not this task's changes; left alone, not staged, not reverted.

No git commit, no git push performed (coordinator commits, per isolation directive).
