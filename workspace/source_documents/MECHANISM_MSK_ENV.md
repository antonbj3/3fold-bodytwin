# MECHANISM MSK ENV — verified working OpenSim environment (2026-07-21)

Executes `docs/MECHANISM_MSK_BUILD_PLAN.md` §5 step 1. Every number below is machine-measured this
session (real command output), not recalled. Isolation respected: nothing written outside
`~/projects/bodytwin`; the two NTFS drives were read-only inputs, mtimes verified unchanged.

## Environment

- **Path:** `source_repository/.venv-msk` (uv-managed, ext4 — NOT the broken NTFS venvs).
- **Python:** 3.12.8 (`uv`'s own managed interpreter — system has no `python3.12` in PATH; uv already
  had `cpython-3.12.8-linux-x86_64-gnu` provisioned locally, so venv creation needed **zero network**).
- **opensim:** `4.6-2026-06-22-85aaf64` (verified via `opensim.GetVersion()`), installed **fully offline**
  from the pre-cached wheel (`uv pip install --offline` resolved+installed in <20ms — confirms the build
  plan's "should take seconds" claim empirically, not just by assertion).
- **numpy:** 2.5.1 (satisfies opensim's `numpy>=2.1` pin; also cached, also offline).
- Repro:
  ```
  ~/.local/bin/uv venv --python 3.12 source_repository/.venv-msk
  PYTHONNOUSERSITE=1 ~/.local/bin/uv pip install --python source_repository/.venv-msk/bin/python --offline opensim numpy
  PYTHONNOUSERSITE=1 source_repository/.venv-msk/bin/python scripts/msk/verify_env.py
  ```
- Verified present in the stock pip wheel (closes Honest Gap #8 in the build plan — previously only
  confirmed in the JAM C++ fork): `PathSpring`, `Blankevoort1991Ligament`, `InverseKinematicsTool`,
  `InverseDynamicsTool`.
- **Not covered by this venv:** JAM-plugin-only tools (`COMAKTool`, `ForsimTool`, `JointMechanicsTool`) —
  those remain C++-only via the separate JAM build at
  `/media/anton/8838D60F38D5FBDE/mechanism_data/opensim_jam_build/opensim-build/opensim-cmd`
  (re-verified working this session: `OpenSim version 4.5.2-2026-01-19-a3c872a2b`, `ldd` clean, exit 0).
  Two independent, both-working OpenSim installs now exist for two different purposes: `.venv-msk`
  (Python API, standard IK/ID/PathSpring/Ligament) vs. JAM binary (COMAK contact mechanics).

## Model

`LaiArnoldModified2017_poly_withArms_weldHand` — scaled instance used:
`/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/OpenSimData/Mocap/Model/LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`
(read-only; 10 subjects total exist under `LabValidation_withVideos/`, same filename pattern).

Detailed-foot joints confirmed by direct grep on the XML (not recalled):
```
3845:  <PinJoint name="subtalar_r">
3903:  <PinJoint name="mtp_r">
4473:  <PinJoint name="subtalar_l">
4531:  <PinJoint name="mtp_l">
```
Wrist gap confirmed (0 DOF, matches build plan §2/#9):
```
5045:  <WeldJoint name="radius_hand_r">
5370:  <WeldJoint name="radius_hand_l">
```
Model load report (35 coordinates, 22 bodies, 22 joints, 80 muscles, 93 actuators, 49 markers) — `.vtp`
mesh-geometry files are absent alongside this copy (non-fatal `Couldn't find file '*.vtp'` warnings at
load; IK/ID don't need visual meshes, only the kinematic tree + markers).

## Smoke test (Task 3) — PASS

Scripts: `scripts/msk/verify_env.py` (env check), `scripts/msk/smoke_test_ik.py` (IK end-to-end +
external-anchor regression check). Output data: `data/msk_smoketest/subject2_walking1/`.

**Method:** patched only the 4 path-bearing tags (`results_directory`, `model_file`,
`output_motion_file`, `marker_file`) of the *original* `walking1_setup_ik.xml` (Windows paths from the
2021 mobilecap/OpenCap run) to real Linux paths — left all 51 `IKMarkerTask` weights byte-identical, so
the comparison isolates "does the rebuilt toolchain reproduce" rather than confounding it with different
task weights. Ran `opensim.InverseKinematicsTool` on subject2's real `walking1.trc` (158 frames, 51
markers, 100 Hz).

**Result:** `InverseKinematicsTool completed 158 frames in 3 second(s)`, output
`walking1_smoketest.mot` (97,420 bytes) — independently confirmed on disk (`ls`, not just the script's
self-report), nRows=158/nColumns=36 matching the reference file's own header exactly. Raw `sed`/`tail`
inspection (bypassing this session's own parser) shows genuine smoothly-varying joint-angle series, not
placeholder/degenerate values.

**External anchor:** compared against the pre-existing `walking1.mot` in the same IK folder — an
independently-computed reference from the original 2021 Windows/mobilecap pipeline, never touched or
regenerated here (mtime verified unchanged: `2021-12-06`). Pre-registered thresholds (set before running,
not fit after): PASS iff median per-column RMSE < 0.5° AND global max abs diff < 5.0°.

| joint | RMSE (deg) | max abs diff (deg) |
|---|---|---|
| knee_angle_r | 0.049 | 0.190 |
| ankle_angle_r | 0.169 | 0.808 |
| subtalar_angle_r | 0.271 | 1.328 |
| mtp_angle_r | 0.044 | 0.111 |
| knee_angle_l | 0.048 | 0.279 |
| ankle_angle_l | 0.137 | 0.564 |
| subtalar_angle_l | 0.203 | 0.879 |
| mtp_angle_l | 0.044 | 0.111 |
| hip_flexion_r/l | 0.110 / 0.109 | 0.439 / 0.376 |

**Across all 35 shared columns: median RMSE = 0.055°, global max abs diff = 1.328° (subtalar_r, the
least marker-observable DOF — physically sensible). Both well inside the pre-registered ceiling — not a
knife-edge pass.** This also proves the detailed-foot DOFs (subtalar/mtp) are not just structurally
present (§ above) but actually **solved for** by IK, tracking an independent reference to sub-1.5°.
Sanity: 0 NaN, only 1/5688 exactly-zero value, strictly monotonic time column.

**Scope honesty:** single trial (subject2/walking1), Mocap markers (not yet the video/Pose2Sim chain).
This is a smoke test per the task definition, not a claim of generality across subjects/trial types.

## Next build step

Per `docs/MECHANISM_MSK_BUILD_PLAN.md` §5 step 2: attach one `PathSpring` (Tier-1 elastic band) to this
same scaled model, anchored to a fixed ground body and an existing foot-segment path point, and smoke-test
in isolation (tension==0 below resting length, linear above, reaction visible in ID moments) — before
touching any real video. `.venv-msk` already has `PathSpring`/`Blankevoort1991Ligament` confirmed present,
so this is a zero-new-dependency next step. In parallel/subsequently: re-provision Pose2Sim on ext4 (same
pattern, separate venv) to re-establish the video→kinematics chain and extend this same regression-check
method to a video-derived trial.

## Blockers / honest notes

- None blocking. System has no `python3.12` in `PATH` — irrelevant, the venv's own `bin/python`/`python3`/
  `python3.12` all symlink directly to uv's managed 3.12.8 interpreter, no PATH dependency for callers that
  use the venv's own binary (as both scripts here do, via absolute path).
- Root partition (`/`, hosts `/home`) is at 95% used, 35G free — this venv cost 339MB, negligible, but
  flag before any large bulk data operation.
- This pip wheel does not include JAM's COMAK/Forsim/JointMechanics tools (plugin-only, C++) — use the
  separate JAM `opensim-cmd` build for that class of work, unaffected by anything in this doc.
