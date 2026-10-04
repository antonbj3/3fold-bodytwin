# MECHANISM — Grand Challenge DM: same-subject in-vivo knee-force validation (FLAGSHIP)

> ## ⚠️ CORRECTED 2026-07-22 — UNITS BUG found; this flagship's headline is RETRACTED
> **The measured anchor below (58–68 %BW) is wrong.** The raw eTibia CSV
> (`data/external/simtk_kneeloads_DM_measured_force/DM_ngait_og1_knee_forces.csv`, columns `Fx,Fy,Fz`
> **unlabeled**) is in **pounds-force (lbf)**; the resultant `sqrt(Fx²+Fy²+Fz²)` was treated as Newtons — the
> ×4.4482 lbf→N conversion was never applied. My "independent re-parse" (63.8 %BW) reproduced the bug because
> it shared that same hidden units assumption — so it was **not** a decorrelated check.
>
> **Decisive tell I missed:** a peak tibiofemoral force of **0.67×BW is physically impossible** (it is *less
> than the ground reaction force*; every instrumented knee — D'Lima / Bergmann / Kutzner — shows ~2–3×BW at
> walking peak). Corrected: **437.9 lbf → 1948 N ≈ 296 %BW** (og1 peak), matching the cohort (JW 271, SC
> 262 %BW; `MECHANISM_GRAND_CHALLENGE_MULTISUBJECT.md`) and the instrumented-knee literature.
>
> **Corrected conclusion:** the twin's muscle-free **Tier-1 (104 %BW) *UNDER*-predicts** DM's ≈296 %BW — exactly
> as a muscle-free lower bound should. The "twin over-predicts DM by 1.5–1.8× at the muscle-free floor" claim is
> a units artifact and is **RETRACTED.** **Crucially, the twin's full deformable-contact model (OpenSim-JAM COMAK)
> already ran on this exact DM gait and predicted 2.503 BW ≈ a ~3% MATCH** to the 2.584 BW measured (`KNEE-CELL`,
> `MECHANISM_STATE.md`, r=0.802) — so on same-subject in-vivo ground truth the twin's *best* knee model MATCHES;
> the muscle-free Tier-1 (104 %BW) here is only the floor. The separate cross-cohort over-prediction (LaiArnold
> **static-optimization free-body-cut** on *subject2* vs OrthoLoad, 391 vs 258 %BW) is a different model /
> subject / quantity (free-body-cut reaction, not deformable contact); it stands on its own but is now suspect
> as an **SO-method** property rather than a twin-physics error, since COMAK contact matches same-subject in-vivo. Everything below is preserved for provenance
> but **SUPERSEDED.**
>
> **Independently confirmed 2026-07-22** (`MECHANISM_DM_UNITS_ADJUDICATION.md`): a neutral adjudication agent —
> instructed to derive its own verdict, not find this one — independently ruled **lbf**. Its decisive,
> mass-and-literature-independent check: the per-trial peak JCF/GRF ratio is **0.51–0.64 under a Newton reading
> (impossible — knee force below the ground reaction) but 2.27–2.83 under lbf (textbook 2–3×) in every one of 7
> trials.** DM mass read from `DM.osim` = 69.98 kg; D'Lima 2012 + Kutzner 2010 literature = 2–3.3×BW. One honest
> residual: the `1legstand` quasi-static plateau fits a Newton reading better — disclosed, outweighed by the three
> gait checks.

**The decisive test the whole force investigation was built toward.** Until now every force comparison was
*cross-cohort* (twin on subject2 vs OrthoLoad medians from other people), leaving one escape hatch open: maybe
the ~1.5× over-prediction was a cohort mismatch, not a real model error. This closes that hatch — for the first
time the model's kinematics AND the measured in-vivo knee contact force come from the **same subject** (Grand
Challenge subject **DM**, instrumented eTibia tibial implant, Fregly/D'Lima).

> **Provenance / salvage note:** the validation agent computed the result and wrote its evidence JSON but stalled
> before writing this doc and was stopped by the coordinator. This report is authored from the salvaged
> `data/msk_smoketest/DM_ngait_og1/final/dm_grand_challenge_results.json`, and the **measured anchor was
> independently re-parsed from the raw eTibia CSV by the coordinator** (63.81 %BW global / 58.18 %BW window —
> exact match). The trial identity is cross-checked (measured GRFz vs independently parsed GRF r=0.996).

## Result: CONFIRMS the over-prediction against same-subject ground truth
| framing | twin (Tier-1, muscle-free) | DM measured in-vivo | ratio |
|---|---:|---:|---:|
| same-stance-phase peak (t=[2.30,3.10]s) — **PRIMARY** | 104.4 %BW | 58.2 %BW | **1.79×** |
| trial-wide peak-vs-peak | 104.4 %BW | 63.8 %BW | 1.64× |
| twin vs cross-trial median (17 trials) | 104.4 %BW | 68.1 %BW | 1.53× |
| exact-same-instant (strictest, phase-sensitive) | 104.4 %BW | 35.6 %BW | 2.93× |

DM's measured in-vivo knee contact peak is **56–83 %BW across 17 trials (median 68 %BW)** — my independent CSV
re-parse confirms 63.8 %BW (global) / 58.2 %BW (validated low-residual window) for the primary trial. The twin
over-predicts by **1.5–1.8×** on same-subject in-vivo ground truth. The pre-registered "matches within 20%"
alternative is **refuted under every framing.**

## Why this is a *conservative* (lower-bound) confirmation
The twin number here is **Tier-1 — the pure external-reaction contact force with NO muscle contribution.** The
full muscle+ligament Tier-2 cert is **BLOCKED on DM** (diagnosed-gap, below). But muscle and ligament crossing
forces only ever **add** to the joint contact magnitude beyond the pure reaction (at subject2's knee/hip, Tier-2
was 2.3–3.8× *larger* than Tier-1). So the true over-prediction, if Tier-2 could run, would most plausibly be
**larger**, not smaller — 1.5–1.8× is a floor.

## Honest caveats (do not over-read)
- **Tier-2 (full SO cert) is blocked on DM**: all 241 frames failed Ipopt ("model too weak for static
  optimization"). Root cause is real and disclosed: DM's TKA knee+patellofemoral carry 12 *free* secondary DOF
  (vs LaiArnold's spline-coupled knee), a 149-bundle Blankevoort ligament apparatus (aggregate passive force up
  to ~18,700 N) the muscles+reserves must balance every instant, and secondary-DOF reserves shipped at
  optimal_force=1 N. Also: DM's quadriceps insert on the patella and reach the knee cut via the 30 patellar-tendon
  ligament bundles, not as direct muscle crossers — a genuine structural difference from LaiArnold's free-body cut.
- **DM is a TKA (prosthetic) knee** — the in-vivo force is through an implant; geometry differs from a native
  knee. It is the best available same-subject in-vivo data, not a native-knee analog.
- IK is 1.4 cm RMS with 23.5% of marker samples occluded (cleaned); the whole-body Y-residual fails *overall* due
  to overground force-plate-coverage gaps (all plates read zero between footfalls) but is **clean (3.46 %BW) at
  the twin's own peak instant**.

## What it changes
The over-prediction is no longer explainable as a cross-cohort artifact **at the knee**: on the one subject where
we have same-subject measured in-vivo force, the twin over-predicts by ≥1.5× even at its muscle-free floor. This
is the strongest structural evidence in the investigation that the twin genuinely over-predicts real in-vivo
knee contact force — consistent with, and now anchoring, the whole `MECHANISM_FORCE_INVESTIGATION_SYNTHESIS.md`
picture. Getting the full Tier-2 number on DM (adapting the SO machinery to the free-DOF TKA knee) is the top
open follow-up; the other 5 Grand Challenge subjects need the operator's free SimTK login.

**Confidence tier:** in-vivo-anchored (same-subject eTibia implant, coordinator-re-verified anchor) for the
Tier-1 ratio; Tier-2 = diagnosed-gap. Files: `data/msk_smoketest/DM_ngait_og1/final/dm_grand_challenge_results.json`,
`data/external/simtk_kneeloads_DM_measured_force/` (measured), `data/external/opensim_jam_grand_challenge_DM/` (model+motion).
