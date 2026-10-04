# MECHANISM MUSCLE AUDIT — LaiArnoldModified2017 muscle-tendon actuators (2026-07-21)

Answers the operator's explicit ask: are this model's muscles "correctly modeled"? Executes via
`scripts/msk/audit_muscles.py`, run for real this session with `.venv-msk`'s OpenSim 4.6 Python API
(`opensim.GetVersion() == '4.6-2026-06-22-85aaf64'`, per `docs/MECHANISM_MSK_ENV.md`). Every number
below is machine-printed console output from that one run (`exit 0`), not recalled or assumed.
Isolation respected: read-only against the two external NTFS data drives, nothing written outside
`~/projects/bodytwin` except the two deliverables (`scripts/msk/audit_muscles.py`,
`scripts/msk/audit_muscles_evidence.json`, this doc) — no git commit, no git push.

## Method

Loaded 4 models, all via `opensim.Model(path); model.initSystem()`:

| role | internal model name (printed) | path |
|---|---|---|
| **SCALED** (primary target, per MECHANISM_MSK_ENV.md) | `LaiArnoldModified2017_poly_withArms_weldHand_scaled` | `.../LabValidation_withVideos/subject2/OpenSimData/Mocap/Model/..._scaled.osim` |
| GENERIC (unscaled nominal template, same lineage) | `LaiArnoldModified2017_poly_withArms_weldHand` | `.../subject2/OpenSimData/Mocap/Model/..._generic.osim` |
| ANCHOR-A (independent artifact) | `Pilot2_Scale3` | `opensim_jam_build/opensim-core-jam/Applications/opensense/test/model_Rajagopal2015_posed.osim` |
| ANCHOR-B (independent artifact) | `subject_walk_scaled` | `opensim_jam_build/opensim-core-jam/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled.osim` |

ANCHOR-A/B are bundled with the separate JAM C++ build tree (opensense test fixture, official Moco
tutorial example) — genuinely different provenance from the LabValidation/OpenCap dataset the scaled
model comes from. **Honesty note on the external anchor:** a live internet fetch of the canonical
`opensim-models` GitHub repo's `Rajagopal_2015.osim` was attempted this session (`curl` reachable,
guessed paths 404'd; one `WebSearch` call to find the right path returned "quota exhausted, 2000/2000
used" — a pre-existing constraint already flagged in `docs/MECHANISM_MSK_BUILD_PLAN.md`). So the anchor
class here is **"two independent on-disk artifacts agreeing,"** not "internet-independently verified" —
disclosed, not hidden. Their internal names (`Pilot2_Scale3`, `subject_walk_scaled`) turned out to be
important: both are **themselves specific scaled subjects**, not nominal references (see Anomaly-check
§3 below — this changed the interpretation of one check mid-audit).

## 1. Top-line counts (real `getSize()` calls, scaled model)

| quantity | value |
|---|---|
| **Muscle-tendon actuators** (`getMuscles().getSize()`, all `Millard2012EquilibriumMuscle`) | **80** |
| **Total actuators** (`getActuators().getSize()`) | **93** |
| Non-muscle actuators (`93 - 80`) | **13**, all `CoordinateActuator` |
| Coordinates / Bodies / Joints | 35 / 22 / 22 |

The 13 non-muscle actuators are ideal torque motors on DOFs with no muscle: `lumbar_ext`,
`lumbar_bend`, `lumbar_rot` (trunk, 3), `shoulder_flex/add/rot_{r,l}`, `elbow_flex_{r,l}`,
`pro_sup_{r,l}` (arms, 10) — all `optimal_force = 10.0` (N·m for these rotational coordinates). This
exactly matches `docs/MECHANISM_MSK_ENV.md`'s previously-logged "80 muscles, 93 actuators" — cross-check
passes.

80 muscles = 40 unilateral × 2 sides (confirmed by name: every muscle has a `_r`/`_l` pair, verified
programmatically). Muscle set (Rajagopal 2016 lower-limb lineage, extended per Lai/Arnold 2017 — model
filename is the standard citation pointer; paper details recalled, not re-fetched this session):
adductors (6: brevis/longus/magnus×3/one more), biceps femoris (lh/sh), toe/ankle group (edl, ehl, fdl,
fhl, tibant, tibpost, perbrev, perlong), gastrocnemius (med/lat), gluteus max/med/min (×3 compartments
each = 9), gracilis, iliacus, piriformis, psoas, rectus femoris, sartorius, semimembranosus,
semitendinosus, soleus, tensor fasciae latae, vasti (int/lat/med).

## 2. Posterior-chain muscles — presence + params table

**16/16 present** (both sides; right side shown, left confirmed symmetric within 1.48% — §4). All are
`Millard2012EquilibriumMuscle`.

| Muscle | Fmax (N) — scaled | Lopt (m) — scaled | Lts (m) — scaled | Pennation (deg) | Fmax (N) — generic | scaled/generic Lopt,Lts ratio |
|---|---:|---:|---:|---:|---:|---:|
| Gastrocnemius medialis (`gasmed_r`) | 3115.51 | 0.0740 | 0.4852 | 9.49 | 3115.51 | 1.254 |
| Gastrocnemius lateralis (`gaslat_r`) | 1575.06 | 0.0865 | 0.4690 | 12.05 | 1575.06 | 1.254 |
| Soleus (`soleus_r`) | 6194.84 | 0.0552 | 0.3527 | 21.85 | 6194.84 | 1.254 |
| Biceps femoris long head (`bflh_r`) | 1313.18 | 0.1192 | 0.4060 | 10.08 | 1313.18 | 1.221 |
| Biceps femoris short head (`bfsh_r`) | 557.11 | 0.1332 | 0.1278 | 15.14 | 557.11 | 1.208 |
| Semimembranosus (`semimem_r`) | 2200.99 | 0.1047 | 0.4078 | 14.59 | 2200.99 | 1.217 |
| Semitendinosus (`semiten_r`) | 591.30 | 0.2367 | 0.3032 | 13.83 | 591.30 | 1.227 |
| Gluteus maximus 1 (`glmax1_r`) | 983.78 | 0.1624 | 0.0964 | 20.28 | 983.78 | 1.105 |
| Gluteus maximus 2 (`glmax2_r`) | 1406.05 | 0.1776 | 0.1233 | 21.05 | 1406.05 | 1.131 |
| Gluteus maximus 3 (`glmax3_r`) | 947.75 | 0.1951 | 0.1203 | 21.91 | 947.75 | 1.168 |
| Gluteus medius 1 (`glmed1_r`) | 1093.47 | 0.0819 | 0.0626 | 18.14 | 1093.47 | 1.069 |
| Gluteus medius 2 (`glmed2_r`) | 765.09 | 0.0902 | 0.0807 | 18.14 | 765.09 | 1.069 |
| Gluteus medius 3 (`glmed3_r`) | 871.20 | 0.0857 | 0.0531 | 18.14 | 871.20 | 1.099 |
| Gluteus minimus 1 (`glmin1_r`) | 374.05 | 0.0866 | 0.0206 | 10.00 | 374.05 | 1.064 |
| Gluteus minimus 2 (`glmin2_r`) | 394.82 | 0.0731 | 0.0341 | 0.00 | 394.82 | 1.064 |
| Gluteus minimus 3 (`glmin3_r`) | 446.77 | 0.0389 | 0.0520 | 1.00 | 446.77 | 1.097 |

(`vmax` = 10.0 fiber-lengths/s, path points ≥ 2 for all 16 — no degenerate paths.) Full 80-muscle dump
is in the script's console output and `scripts/msk/audit_muscles_evidence.json`; omitted here for
length (LEAN — the operator's explicit table ask is the posterior chain).

## 3. Physiological plausibility check — 0/80 flagged

Pre-registered bounds (broad, literature-level, set before inspecting any of this model's numbers):
Fmax ∈ [1, 10000] N, Lopt ∈ [5, 600] mm, Lts ∈ [0, 550] mm, pennation ∈ [0°, 45°], vmax ∈ [1, 20]
fiber-lengths/s, path points ≥ 2. **Result: 0/80 muscles flagged** — no NaN, no non-positive/degenerate
value, no muscle outside these bounds anywhere in the model. Machine-checked (`bounds_flags()` in the
script), not eyeballed.

**Forcing the adversary on this negative** (per the "honest-negative is not a free pass" rule — a
one-shot "nothing flagged" is not enough): three follow-up checks, since loose bounds could be masking
a subtler problem:

- **Scaled-vs-its-own-generic ratio, all 80 muscles** (the clean, confound-free test — same lineage,
  same template, only subject2's own scale factors differ; pre-registered bound [0.5, 2.0]×): **Fmax
  ratio = exactly 1.000000 for all 80** muscles, **pennation angle diff = exactly 0.0° for all 80** —
  i.e. this scaling pipeline **never rescales strength or pennation**, only geometry. Lopt/Lts scale
  together by an identical per-muscle factor in **[1.06, 1.25]×** — sane for one adult subject vs a
  nominal template, **zero muscles outside the pre-registered bound**.
- **Cross-anchor divergence** (generic vs ANCHOR-A vs ANCHOR-B, threshold 15%): raw result = 4
  (muscle, parameter) pairs flagged, all `fhl` (flexor hallucis longus) Lopt/Lts, bilateral, ~15-16%.
  **Adversary forced and it survives as an explanation, not the "bug" hypothesis**: both anchors are
  themselves independently-scaled subjects (`Pilot2_Scale3`, `subject_walk_scaled`, not nominal
  references — see §Method), and the divergence is not isolated to `fhl` — it's a graded pattern across
  the *entire* distal shank/foot muscle group (`edl` 13.5%, `ehl` 14.2%, `fdl` 14.7%, `fhl` 16.3%,
  `tibant` 10.1%, `tibpost` 8.8%, `perbrev` 8.3%, `perlong` 9.7% vs ANCHOR-A alone), consistent with
  ordinary inter-subject anthropometric variation concentrated in the shank/foot segment, not a
  parameter-entry defect. Meanwhile Fmax and pennation showed **zero** divergence across all 80 muscles
  in this same 3-way comparison — a second, independent confirmation (beyond §3's own-lineage test) that
  those two parameters are universal constants in this model family, never subject-scaled by any of the
  3 independent scaling runs checked. **Conclusion: no anomaly here, retracted after forcing the
  confound** — reported for transparency (kills/near-misses are auditable), not as a live flag.
- **Bilateral (L/R) symmetry**: GENERIC max diff = 0.55% (a template should be ~symmetric — passes);
  SCALED (subject2) max diff = 1.48% (worst: `soleus_r` Lopt vs `soleus_l`) — small, physiologically
  plausible subject-specific asymmetry from marker-based scaling, not a red flag.

## 4. Anomalies / gaps

1. **Erector spinae / multifidus / longissimus / iliocostalis: ABSENT.** Keyword search across all 80
   muscle names in all 4 loaded models (`scaled`, `generic`, `anchor_a`, `anchor_b`) for
   erector/longissimus/multifid/iliocost/spinae/semispinalis/obliq/rectusabd/quadratuslumb/transvers
   returns **zero hits in every one of the 4 models**. Trunk extension/bending/rotation is driven
   entirely by 3 ideal `CoordinateActuator`s (`lumbar_ext`, `lumbar_bend`, `lumbar_rot`,
   `optimal_force = 10 N·m` each) — **no discrete muscle-tendon actuator represents the lumbar
   extensors anywhere in this model lineage.** This is a real architecture gap, not a parameter bug,
   and it directly answers the operator's explicit question about erector spinae/lumbar extensors: they
   are not modeled as muscles at all. (Caveat on the 10 N·m figure: this is far below published maximum
   voluntary trunk-extension moments, ~150-300 N·m per Marras/McGill-era literature — recalled from
   training-data knowledge, not independently re-fetched this session. The low value is plausibly a
   deliberate placeholder that biases muscle-driven optimization tools like Static Optimization/CMC/Moco
   away from relying on this ideal actuator, rather than a claimed physiological torque ceiling — this
   interpretation is inference, not confirmed from a source document, and is flagged as such.)
2. **No subject-specific muscle strength scaling.** §3's ratio check shows `Fmax` is bit-identical
   between subject2's scaled model and the nominal generic template for all 80 muscles. If the operator
   expects subject2's muscles to reflect subject2's actual strength/size (vs. population-average
   cadaver-derived Fmax), that is **not what this scaled model contains** — this is standard/common
   OpenSim practice (many published gait studies use exactly this convention) but is a real modeling
   limitation worth being explicit about given the "correctly modeled" framing of the ask.
3. **Cross-anchor `fhl`/distal-foot-muscle divergence** — investigated in §3, explained by inter-subject
   variation once the anchors' own scale-provenance was checked, not a live anomaly. Listed here only
   for the record (symmetric QC: a forced-and-killed flag is disclosed, not silently dropped).
4. No anomalies found in: posterior-chain presence (16/16), physiological bounds (0/80 flagged), path
   degeneracy (0 muscles with <2 path points), NaN/degenerate values (0/80), bilateral symmetry (both
   models pass), or scaling sanity (0/80 muscles outside the pre-registered ratio bound).

## Next step

The muscle-tendon actuator layer of this model checks out against every test run this session: full
posterior-chain coverage, zero physiologically-implausible parameters, and a scaling pipeline that
behaves exactly as it should (geometry scaled sanely per-subject, strength/pennation held at the
literature constant). The one real, actionable gap is **#1 above** — if the operator's downstream work
needs a lumbar-extensor mechanism with realistic force capacity (e.g. anything coupling to the
posterior-chain/erector-spinae interest driving this audit), the current 10 N·m ideal actuator will not
supply a realistic trunk-extension moment and either needs (a) a literature-sourced
`Millard2012EquilibriumMuscle`/`Blankevoort1991Ligament` erector-spinae addition (same zero-new-
dependency pattern already used for the Tier-1/2 elastic-band work in
`docs/MECHANISM_MSK_BUILD_PLAN.md` §3 — `scripts/msk/attach_band.py` is the sibling in-progress example
of adding a new path-actuator to this exact model), or (b) an explicit note in any spine-loading result
that the ideal actuator's 10 N·m ceiling is a known non-physiological placeholder. Separately, closing
Honest Gap in this audit itself: the internet-independent literature anchor (canonical
`opensim-models` GitHub `Rajagopal_2015.osim` or the published Rajagopal 2016 parameter table) was not
reachable this session (WebSearch quota exhausted) — worth a fresh attempt in a future session with
quota available, though the two independent on-disk anchors used here already gave an over-determined
check (2 independent artifacts agreeing) for the claims made.
