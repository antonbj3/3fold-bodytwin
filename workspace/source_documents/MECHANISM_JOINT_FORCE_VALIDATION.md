# MECHANISM JOINT FORCE VALIDATION — first predicted-vs-in-vivo cert (2026-07-21)

## ℹ POINTER (2026-07-21) — this doc's own number is NOT affected, but three downstream docs were

A sign bug was found and fixed in `scripts/msk/static_opt_knee.py`'s
`knee_crossing_muscles_and_forces` (full story: `docs/MECHANISM_SIGN_BUG_AUDIT.md`,
`docs/MECHANISM_SIGN_BUG_REMEDIATION.md`). **This doc's own headline number (102.17 %BW, PURE
kinematics+GRF reaction, no muscles) is UNAFFECTED** — `validate_joint_force.py` never calls that
function; the bug lives entirely in the downstream muscle-subtraction step this doc's own method
explicitly does not perform (§3: "never calls `InverseDynamicsTool`/`InverseDynamicsSolver`/
`JointReaction` — cannot be touched"; the crossing-muscle detector is a separate, later addition
in the docs listed below, not present here). The three docs that build on this one by ADDING
muscle forces did have the bug, now corrected — see their own top banners:
`docs/MECHANISM_STATIC_OPT.md` (knee: 233.20 → 391.10 %BW), `docs/MECHANISM_HIP_FORCE.md` (hip: 235.30
→ 386.77 %BW), `docs/MECHANISM_ANKLE_FORCE.md` (ankle: 286.68 → 484.49 %BW).

---

Executes the operator's task: does the twin's PREDICTED knee joint force match in-vivo OrthoLoad,
for a self-contained case needing no new data? Every number below is machine-measured this session
(`scripts/msk/validate_joint_force.py`, exit 0), not recalled. Isolation respected: `.venv-msk` only,
OrthoLoad + LabValidation data read in place (never written to), no git commit/push.

## Headline result

| | peak (%BW) |
|---|---:|
| **PREDICTED knee reaction force** (this twin, subject2 `walking1`) | **102.17** |
| **IN-VIVO OrthoLoad knee** (median, n=72 unassisted level-walking trials, 9 subjects) | **258.22** |
| discrepancy | 156.05 percentage points (ratio 0.396) |

**Both numbers were independently sanity-checked against the task's own stated ballpark (~2.5–3×BW
in-vivo) before comparing them to each other:** the OrthoLoad median (258.2%BW, range 186.7–353.6%)
lands exactly in that band — the external anchor itself is sane. The prediction (102.2%BW) is
**expected to be lower, not a pipeline defect** — see §3, the single most important section of this
document.

## 1. Inputs (all pre-existing, nothing new fetched)

- Model: `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` (subject2, build step 1; 78.2 kg
  total mass, exact match to `sessionMetadata.yaml`'s `mass_kg: 78.2` — confirms the scaling preserved
  measured mass exactly).
- Kinematics: `OpenSimData/Mocap/IK/walking1.mot` — 158 frames, 100 Hz, 0–1.57 s (the same trial
  anchored to <1.5° in the Step-1 IK smoke test, `docs/MECHANISM_MSK_ENV.md`).
- Ground reaction force: `ForceData/walking1_forces.mot` — **present** (both feet, 2000 Hz, 0–1.579 s,
  same time origin as the IK file). The task's "if GRF is missing" fallback does not apply here.
- OrthoLoad in-vivo knee anchor: `data/external/orthoload/knee/database_api/akf/*.akf` (609 files).

## 2. Pipeline

1. **Official Inverse Dynamics** (`opensim.InverseDynamicsTool`, patched from the pre-existing
   `walking1_setup_id.xml` / `..._setup_externalLoads.xml` — path-patch only, same discipline as
   `scripts/msk/smoke_test_ik.py`). Run for the record; **its `knee_angle_r/l` and `hip_flexion_r/l`
   generalized-force columns are never read** — see §4. Uncorrupted columns are reported as context:
   `ankle_angle_r_moment` peak 106.8 N·m, `lumbar_extension_moment` peak 28.6 N·m (physiologically
   sane order of magnitude for walking).
2. **Predicted joint reaction force** — Newton's second law for a system of rigid bodies (§3), computed
   from forward kinematics + body mass properties + the measured GRF only. Never calls
   `InverseDynamicsTool`/`InverseDynamicsSolver`/`JointReaction` — cannot be touched by the §4 bug.
3. **OrthoLoad anchor** — parse all knee AKF files, select "Level Walking" trials, exclude gait-aid-
   assisted trials, report the full distribution (§5).
4. **Comparison + honesty gates** (§6).

## 3. What "predicted joint reaction force" IS and IS NOT (read this before the numbers)

For the sub-chain of rigid bodies distal to a joint — found by BFS on the model's **own** joint tree,
not assumed (this model's `patellofemoral_r` joint hangs `patella_r` off the **femur**, not the tibia;
a hardcoded "shank = tibia+patella+foot" guess would have been wrong) — Newton's second law gives:

```
sum_i(m_i * a_i) = gravity + GRF + R        (R = the net force crossing the joint boundary)
```

`a_i` (each body's ground-frame COM acceleration) comes **purely from forward kinematics** — the
already-IK-validated joint angles, smoothed + double-differentiated (Savitzky-Golay, validated on a
synthetic known-derivative signal and via a whole-body residual check, §6) — no dynamics engine, no
actuator/muscle model, no assumption about how the net moment is produced.

**This is the classical biomechanics "link-segment joint reaction force" (Winter's textbook method) —
and it is a well-documented, systematically SMALLER quantity than the true bone-on-bone / instrumented-
implant "contact force" that OrthoLoad measures.** The reason: every muscle and ligament crossing the
joint is invisible to a kinematics+external-force-only calculation — Newton's law only "sees" their
combined NET effect on the segment's acceleration. Concretely: both the quadriceps (via the patellar
tendon) and the hamstrings pull the tibia **proximally** when active — their MOMENT contributions
oppose and partially cancel (that's how a net flexion/extension moment gets produced at all), but their
FORCE contributions do not cancel the same way, so real antagonist co-contraction compresses the joint
substantially more than the net reaction force implies. Resolving that requires knowing individual
muscle forces (EMG-informed or Static-Optimization/CMC muscle-driven simulation) — explicitly descoped
here (§8), exactly as the task anticipated ("HONEST caveats... joint-reaction vs contact force").

This also explains, without needing to be tuned to match, why the ratio (0.396) is physically
reasonable rather than suspicious: it is in the same regime as this repo's own COMAK/Grand-Challenge
flagship needing a full muscle-driven pipeline to close an analogous gap (`bt_memory/knee-comak-
flagship-close-and-fork-build.md` — uncalibrated stock-weights COMAK still landed within ~3% of
measured *because* it used muscle forces; a pure kinematic reaction force has no such term at all).

## 4. A real OpenSim 4.6 engine bug, independently re-confirmed this session (not just cited)

`docs/MECHANISM_MSK_ELASTIC_BAND.md` §4 documents a machine-verified `InverseDynamicsTool` defect on
this exact model: `knee_angle_r/l` and `hip_flexion_r/l` generalized-force output is wrong by
~1000–1800× at that document's hip-hinge test pose (independent `dPE/dq` virtual-work cross-check:
verified −0.69 N·m vs ID's reported −1252.62 N·m for knee; verified 1.31 N·m vs ID's 1337.12 N·m for
hip). Before trusting this as policy, I re-ran the diagnosis myself this session rather than accepting
it at face value:

- At the model's literal all-zero pose, gravity only, `InverseDynamicsSolver` gave SMALL, plausible
  knee/hip values (−1.27 / +1.29 N·m) — NOT reproducing the bug, which initially looked like a
  refutation.
- Tracing the discrepancy to `scripts/msk/attach_band.py` / `band_config.json` showed the documented
  bug pose is the **hip-hinge posture** (pelvis_tilt 25°, hip_flexion 45°, knee_angle 20°), not the
  zero pose. Re-running `InverseDynamicsTool` at that exact pose reproduced `ankle_angle_r` = 1.170
  N·m and `lumbar_extension` = −7.612 N·m — matching the doc's reported ID values (1.17, −7.61) almost
  exactly, confirming this reproduction's fidelity, though a residual GRF-construction difference left
  the knee/hip figures themselves not bit-identical to the doc's specific table.
- Net conclusion: the defect is pose-dependent (consistent with its likely root cause — a nonlinear
  `PolynomialFunction`-coupled "rolling knee" `CustomJoint` — being most wrong away from the
  fully-extended reference configuration), confirmed independently enough to justify inheriting the
  documented workaround as policy for this model family, exactly as that document's own §8 recommends
  for "any FUTURE script."

**Consequence for this script:** the predicted-force computation in §3 never calls
`InverseDynamicsTool`/`InverseDynamicsSolver`/`JointReaction` at all — it is immune to this defect by
construction, not by luck.

## 5. OrthoLoad anchor selection — a real bug caught and fixed in this script's own first draft

AKF format (verified live): `BodyWeight [N]:` header + a `Time Fx Fy Fz F Mx My Mz Marker` data block,
`Comment #1` naming the activity (e.g. `"Knee Joint; Walking: Level walking"` or `"Knee Joint;
Gaitanalysis; Level Walking;"`).

- **First-draft bug**: filtering `Comment #1` on the literal substring `"Walking: Level walking"`
  found only 22 files. Cross-checking against a differently-coded OrthoLoad parser built this session
  by a sibling MSK agent (`data/external/orthoload/_index/orthoload_akf_trials.jsonl`, 3942/3942 files
  parsed) surfaced trials from subjects (`k6l`, `k7l`, `k9l`) that never appeared in my 22 at all.
  Root cause: some files phrase the activity as `"Gaitanalysis; Level Walking"` with no `"Walking:"`
  prefix — the narrow pattern silently undercounted. Fixed to a broad case-insensitive `"Level
  Walking"` substring match → **72 files across 9 subjects**.
- **Second bug caught by the same fix pass**: `BodyWeight [N]:\t902,5` — 28/609 knee files use a
  **comma** decimal separator; a naive `float()` call crashes. Fixed (`.replace(",", ".")`).
- **Gait-aid exclusion**: 8 "Level Walking" trials are annotated `"3 Point"` / `"4 Point"` — clinical
  gait-pattern terms for **walker/cane-assisted** ambulation, a different, lower-load activity (they
  read 136–238%BW, consistently lower — physically sensible, not cherry-picked to fit). Excluded from
  the primary 72-trial comparison set; reported separately.
- **Cross-validation, not just a citation**: this script's own AKF parser was re-run independently
  against all 72 selected files and compared to the sibling index's `peak_force_pctBW` field —
  **0 mismatches >0.5%BW**. Internally, `sqrt(Fx²+Fy²+Fz²)` recomputed from the raw force components
  matches each file's own printed `F` (resultant) column to within 0.0109 N max across all 72 trials.
- This script is **self-contained**: it re-parses the raw AKF files itself and does not depend on the
  sibling index file continuing to exist; the cross-validation is documented here as supporting
  evidence, not a runtime dependency.

Representative (closest-to-distribution-median, not cherry-picked for fit) trial: `k3r_291008_1_27p`,
259.25 %BW, body weight 960 N, comment "Knee Joint; Gaitanalysis; Level Walking;".

## 6. Machine-checked gates (PASS/FAIL, not eyeballed)

| gate | pre-registered threshold | measured | verdict |
|---|---|---:|---|
| GRF sanity (peak vertical R/L GRF) | 85–140 %BW | 108.94 / 107.64 %BW | PASS |
| Whole-body Newton residual (RMS, x/y/z) | < 8 %BW | 2.75 / 4.10 / 1.17 %BW | PASS |
| Sensitivity to smoothing window (70–130 ms) | knee peak range < 10 pct-pts | 0.13 pct-pts | PASS |
| Regime sanity: predicted/in-vivo ratio | in (0.10, 1.05) — positive, below (not above/equal to) the contact-force anchor | 0.396 | PASS |

The **whole-body** residual check is independent of the knee/hip sub-chain choice: total mass ×
whole-body-COM-acceleration must equal total GRF (both feet) + gravity, using the same
kinematics-differentiation pipeline. It passing at 1–4%BW RMS (not the knee/hip-specific calculation)
is what gives confidence the differentiation + GRF-time-alignment machinery itself is sound, before
trusting it for the harder-to-externally-verify sub-chain numbers.

Bonus consistency checks (reported, not gated): ankle_r reaction peak 107.50%BW ≥ knee_r 102.17%BW ≥
hip_r 88.59%BW — each more-proximal cut sheds roughly that segment's own weight, the expected signature
for modest-acceleration stance-phase walking. Left knee peak 100.29%BW vs right 102.17%BW — a
same-order-of-magnitude L/R symmetry check (not exact, real gait/anthropometric asymmetry expected).

## 7. Honest caveats (full list)

1. **Reaction force ≠ contact force (§3)** — the primary, load-bearing caveat. The 0.396 ratio is the
   expected direction and rough order of magnitude for "muscle co-contraction missing," not a precision
   estimate of the gap; closing it needs Static Optimization or EMG-informed muscle force estimation
   (explicitly descoped here, see §8).
2. **Different subjects.** OrthoLoad's 9 subjects are instrumented-knee-implant patients (elderly,
   post-arthroplasty anthropometry/gait pattern); subject2 is a healthy young(ish) OpenCap participant
   (78.2 kg, 1.96 m). Neither population nor anthropometry is matched — this is a plausibility/ballpark
   comparison, not a per-subject validation.
3. **Different activities within "walking."** OrthoLoad trials are overground level walking at each
   patient's own comfortable pace on their own recording day; subject2's `walking1` is one OpenCap
   capture session trial. Speed/cadence were not matched.
4. **GRF source**: real force-plate data for subject2 (not estimated) — this is the "best case" the
   task anticipated, not the degraded "kinematics + segment mass only" fallback.
5. **The OpenSim ID Tool bug (§4)** is a real, independently-reconfirmed defect in this exact
   model/joint style; the FORCE result here is architecturally immune (never calls the affected code
   path), but this script's own moment context (ankle/lumbar) still leans on the tool for those two
   uncorrupted columns.
6. **Numerical differentiation**: Savitzky-Golay smoothing (11-sample/110 ms window, 4th-order local
   polynomial) is a design choice; shown stable across a 70–130 ms sweep (§6) but is not literally the
   raw signal.
7. **Single trial, right leg primary** (left leg reported only as a coarse symmetry check, not a
   full independent replicate).

## 8. Next step

Static Optimization (or full muscle-driven / COMAK-style analysis) using this document's verified
reaction-force time series as a lower-bound target, to estimate individual muscle forces and recover
the co-contraction compression term missing from §3 — the documented, understood path to close the
0.396 ratio toward 1.0, per `bt_memory/knee-comak-flagship-close-and-fork-build.md`'s own account of
what that step costs (SO is known-finicky; a naive linear-EMG-informed attempt there made things
*worse*, not better — a real, bounded, non-trivial follow-on task, not a quick fix).

## Files

- `scripts/msk/validate_joint_force.py` — the full pipeline (self-contained, re-runnable).
- `data/msk_smoketest/subject2_walking1/joint_force_validation/` — patched ID setup XMLs, ID `.sto`
  output, `joint_force_validation_results.json` (every number in this document, machine-written).
