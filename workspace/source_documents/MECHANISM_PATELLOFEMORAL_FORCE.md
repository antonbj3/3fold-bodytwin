# MECHANISM PATELLOFEMORAL JOINT (PFJ) CONTACT FORCE — a genuinely new joint (2026-07-21)

**Question.** The mission's force certs so far cover tibiofemoral ("knee"), hip, ankle, spine —
not the **patellofemoral joint** (patella against the femoral trochlea, loaded by the quadriceps
via the patellar tendon), a distinct articulation from the tibiofemoral joint already certified.
This adds it: peak PFJ contact force over subject2/walking1, via two independent, decorrelated
pathways, tested against the literature and against the task's own named symmetric-QC concern —
does the patellar-tendon:quadriceps-tendon force ratio (van Eijden 1986; Buff/Jones/Hungerford
1988) actually emerge from each model's own mechanism, or is it silently assumed 1:1?

**Which model has a patella — the task's own framing was WRONG, verified live.** The task assumed
the base LaiArnold subject2 model "likely lacks" a patella/PFJ. Direct XML inspection of
`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` shows it has a `patella_r` body (line
954) and a `patellofemoral_r` `CustomJoint` (line 3646) — it DOES have a PFJ. **Both** models used
here have one, but built completely differently:

| | LaiArnold subject2 (Pathway A) | lenhart2015.osim / JAM (Pathway B) |
|---|---|---|
| Patella DOF | `CoordinateCouplerConstraint` — a function of `knee_angle_r`, not free | 6 FREE coordinates (`pf_flex/rot/tilt/tx/ty/tz_r`) — settle via real forward dynamics |
| Quad→patella→tibia | `recfem_r`/`vasmed_r`/`vaslat_r`/`vasint_r` each thread THROUGH patella via-points and continue, as **one** `Millard2012EquilibriumMuscle` `PathActuator`, straight to `tibia_r` — **zero** separate patellar-tendon element (grep for Ligament/PathSpring: 0 hits) | Same 4 muscles **terminate at** patella_r; a **separate** 6-bundle `Blankevoort1991Ligament` group `PT1`..`PT6` (patellar tendon) carries force from patella to `tibia_r` |
| Consequence for the quad:PT ratio | **Structurally forced to 1.000, always** (one PathActuator = one tension along its whole length — a mathematical certainty, not a simulation result) | Genuinely free to differ (patella is a rigid body in force balance among quad pull + PT tension + trochlear contact + retinacular ligaments) — the same conceptual architecture van Eijden's own 1986 model uses |

Both are used, as two decorrelated pathways; **Pathway A is the primary/higher-confidence
number** (verified below to reproduce the established tibiofemoral cert bit-exactly and to use
real per-frame Static-Optimization tension); **Pathway B is architecturally the "correct" model
for testing the ratio question but its absolute magnitude is a likely floor** (a newly-discovered,
disclosed artifact in the reused data — Sec. 5).

---

## Headline

| | Pathway A (LaiArnold, JointReaction) | Pathway B (JAM/lenhart2015, deformable contact) |
|---|---:|---:|
| Peak PFJ force | **591.79 N = 77.17 %BW** | **331.59 N = 54.54 %(lenhart's own BW)** |
| Time / knee angle at peak | t=1.460s, knee_angle_r=**21.53°** | t=0.760s, knee_flex_r=**40.60°** (swing) |
| Survives forced adversary? | YES — edge-margin-exclusion applied (Sec. 3), peak unaffected | YES — settle-window corrected from area-based (t=0.06s) to force-based (t=0.09s); peak at t=0.76s is far from either boundary |
| vs pre-registered walking band (50–150 %BW) | **inside** (mid-upper) | **inside** (mid) — but likely an under-estimate (Sec. 5) |
| Quad:PT ratio (van Eijden/Buff symmetric QC) | **≡ 1.000 always** (structural, by construction) — FALSIFIES "reproduces the literature ratio" for this architecture | **varies genuinely**, load-bearing-regime median **1.152** (n=110 frames), low-flexion (<20°) median **1.152** vs Buff 1988's 0.86 @ 10° |

**Neither pathway over-predicts** the task's own stated walking-gait PFJ band — a genuine,
disclosed **difference** from the tibiofemoral joint (which over-predicted OrthoLoad by
1.4–2.0×, `docs/MECHANISM_CROSS_SUBJECT.md`). PFJ force's own peak IS, however, concentrated
almost exactly where quadriceps load is highest (Sec. 4) — the same qualitative mechanism as the
tibiofemoral finding (`docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`: recfem_r dominant quad
contributor), just without an accompanying magnitude excess here.

**Confidence tier: cadaveric/published-plausibility** (not in-vivo-anchored) — verified live this
session that **no genuine in-vivo instrumented PFJ force measurement exists in the literature**
(Sec. 6); the task's own given PMID for "Kutzner 2013 PFJ" is wrong.

Script: `scripts/msk/patellofemoral_force.py`. Data:
`data/msk_smoketest/patellofemoral_force/patellofemoral_force_results.json` (+
`patellofemoral_force_series.json`, full per-frame series for both pathways). **No new
simulation was run** — both pathways extract already-existing, previously-computed data
(`walking1_JointReaction_ReactionLoads.sto`, already used by
`docs/MECHANISM_BONE_STRESS_CONSEQUENCE.md`; `gait_driven.h5`, already used by
`docs/MECHANISM_JAM_CONTACT_DECORR.md`) — this document is a **new, deeper extraction** from data
that already existed, not a re-run.

---

## 1. Pre-registration (thresholds fixed before the final numbers were read)

- **Walking-gait PFJ band**: 50–150 %BW — the task's own stated classic-literature range
  (Reilly & Martens 1972 family). Its exact primary-source number could not be independently
  re-verified this session (Sec. 6) — used as the task's own pre-registered band, flagged, not
  fabricated.
- **"Materially over-predicts"**: > 2.0× the band's upper edge (300 %BW) — mirrors the
  tibiofemoral cert's own over-prediction framing.
- **Absolute MVC ceiling** (sanity, not walking): 8000 N — van Eijden 1987's own maximal-
  voluntary-contraction quadriceps force at ~75° flexion (PMID 3630619). Walking-gait PFJ force
  must sit well below this or something is badly wrong (both pathways do: 591.8 N and 331.6 N).
- **Buff 1988 ratio anchor**: FQ/FP (quadriceps tension / patellar-tendon tension) = 0.86 at 10°,
  1.55 at 70° (PMID 3339023, fetched verbatim, Sec. 6). Pass criterion for Pathway B: the
  load-bearing-regime ratio at low flexion (the only part of Buff's curve level-walking stance can
  test — Sec. 5) falls within a **factor of 2** of Buff's 10° value (band: 0.43–1.72).
- **Falsifier for "the model's patellar mechanism reproduces the literature ratio"**: ratio ≡ 1.0
  identically (no angle-dependence at all) — checked structurally (not by eyeballing a plot) via
  (a) exhaustive count of Ligament/PathSpring-class `Force` elements in the whole model (expect >0
  for a model that CAN reproduce it) and (b) live inspection of the quad muscles' own
  `GeometryPath` body sequence.

---

## 2. Method — geometric, reusing already-validated infrastructure

**Pathway A** parses the SAME `walking1_JointReaction_ReactionLoads.sto` (OpenSim's own
`JointReaction` analysis, subject2/walking1, `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`)
that `docs/MECHANISM_BONE_STRESS_CONSEQUENCE.md`'s free body already reads one row of
(`patellofemoral_r_on_patella_r_in_patella_r_f{x,y,z}`) — this column has existed all along but its
own PEAK, %BW, and flexion-angle-at-peak were never previously searched for or reported as a
headline (only its value at the KNEE's own peak instant was used, inside a different free body).
`vjf.parse_mot`/`MODEL_FILE`/`G` reused unedited from `validate_joint_force.py`, same convention as
every other joint-force cert in this mission.

**Pathway B** reads the SAME `gait_driven.h5` (JAM's `ForsimTool`+`JointMechanicsTool` output,
`docs/MECHANISM_JAM_CONTACT_DECORR.md`) — its `pf_contact` group (`Smith2018ArticularContactForce`,
patella_cartilage vs femur_cartilage mesh-side conventions) already existed but only its mean
PRESSURE was previously reported; its own peak FORCE, and the `Blankevoort1991Ligament` PT-bundle
tensions alongside the `Muscle` group's per-muscle actuation, are new extractions from this
already-existing, unedited file. `jam_contact_decorr.py`/`jam_contact_decorr_analyze.py`'s own
`load_series`/`detect_settle_index`/path constants reused unedited.

Both extractions are read-only against externally-fetched or previously-computed data; this
script performs zero new heavy simulation.

---

## 3. Forced adversary #1 — is the "peak" an edge/settling artifact? (both pathways checked)

A first pass found suspiciously convenient-looking peaks that were forced to their strongest fair
form before being trusted, per the discipline (a false positive is a premature negative in
reverse — an unforced adversary):

**Pathway A**: raw argmax landed at t=1.460s (93% through the 1.57s trial) — an interior,
non-edge frame. But a companion check (quadriceps SO-tension peak) landed at the trial's literal
**last frame** (t=1.570s) — exactly the Savitzky-Golay-differentiator edge-artifact class
`bone_stress.py` already discloses and guards against via `EDGE_MARGIN_FRAMES`
(`=SAVGOL_WINDOW_PRIMARY//2` frames excluded each end). Reused that SAME constant here (not a new
heuristic): the PF-force peak (t=1.46s, index 146 of 158) sits safely inside the interior
`[5, 153)`, survives exclusion unchanged; the quad-tension "peak" at the literal last frame does
NOT survive (excluded) — its properly-excluded interior peak is at t=0.060s (975.3 N, vs the
excluded raw value 1010.4 N — the tail-end bump is disclosed as the artifact it looks like: quad
tension plateaus at 938–970 N from t=1.44–1.56s then ticks up only at the very last sample).

**Pathway B**: raw inspection of the FIRST 20 frames (before any exclusion) shows a violent
disengage/re-engage transient BOTH `tf_contact` and `pf_contact` share: at t=0.01s, contact area
collapses to near-zero for both (tf: 2.8 mm², pf: 0.0 mm² — total disengagement), then at t=0.03–
0.04s `tf_contact` force **overshoots to 2007–2081 N** (within 1% of the run's own eventual
"real" 2099.9 N peak, at completely the wrong time) before settling. The REUSED
`detect_settle_index` (from `jam_contact_decorr_analyze.py`, which probes contact AREA) declared
this settled by frame 6 (t=0.06s) — too early: force keeps ringing non-monotonically through frame
8–9 (tf: 1246→1474→1687→1648 N). **Fixed at the source**: reran `detect_settle_index` on FORCE
MAGNITUDE (not area) for both contacts and took the max of all four candidate indices — moved the
cutoff from t=0.06s to **t=0.09s**. The reported peak (t=0.76s, mid-swing) sits far from either
trial boundary and is unaffected by this correction — but the earlier (pre-fix) "peak" at t=0.08s
(367.66 N) WAS a tail-of-ringing artifact and has been discarded, not reported.

Both checks are machine-recorded (`laiarnold_pf_peak_survives_edge_margin_exclusion`,
`jam_settle_window_force_based_ge_area_based` — both PASS in the evidence JSON).

---

## 4. Timing — PFJ's peak does NOT coincide with the tibiofemoral joint's own peak, and this is mechanistically explained, not just observed

At t=0.51s — the already-established tibiofemoral (`knee_joint_contact`) peak instant (391.11 %BW,
exactly reproduced here, Sec. 5) — the SAME trial's PFJ force is only **236.45 N (30.8 %BW)**,
well below PFJ's OWN peak (591.79 N at t=1.46s). Reading `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`'s
own per-muscle breakdown explains why: the tibiofemoral peak at t=0.51s is dominated by
**gastrocnemius** (48.7% of the muscle-driven contribution; a muscle that does NOT cross the
patella at all) with quadriceps a secondary contributor (22.2%, `knee_angle_r`≈6° at that instant
— near full extension). The patellofemoral joint, by contrast, is loaded almost exclusively by the
**quadriceps** (it has no other tendon route) — so its own force peak tracks quadriceps tension's
own timing, not the tibiofemoral joint's gastrocnemius-dominated one.

Confirmed directly: quadriceps SO tension (`recfem_r`+`vasint_r`+`vaslat_r`+`vasmed_r`) at the
PFJ-force peak instant (t=1.46s) is **968.4 N** — within **0.7%** of quadriceps tension's OWN
interior-excluded maximum anywhere in the trial (975.3 N, at t=0.06s, the analogous
loading-response event one gait cycle earlier). **PFJ force IS concentrated almost exactly where
quadriceps load is highest** — directly answering the task's own question, and the same
qualitative mechanism (recfem/vasti dominance) already found for the tibiofemoral joint in
`MECHANISM_CONTACT_MUSCLE_DECOMP.md`, here without an accompanying over-prediction (Sec. 6).

---

## 5. Symmetric QC — the quad:PT ratio, forced (not assumed 1:1)

**LaiArnold (Pathway A)**: verified live via the OpenSim API — `recfem_r`'s `GeometryPath` touches
`['pelvis', 'patella_r', 'patella_r', 'patella_r', 'tibia_r']`, i.e. one continuous
`Millard2012EquilibriumMuscle`. An exhaustive scan of the model's entire `ForceSet` (93 forces)
found **zero** `Ligament`/`PathSpring`-class elements anywhere. Consequence, a mathematical
certainty (not a simulated result): the "quadriceps-tendon" and "patellar-tendon" segments of this
single actuator carry **identical** tension, always — ratio ≡ 1.000. **This is the model's own
falsifier for "reproduces the angle-dependent literature ratio" — it cannot, by construction.**
Disclosed as a genuine, clean limitation of this specific (very common, e.g. gait2392-style)
modeling convention, not a bug.

**JAM/lenhart2015 (Pathway B)**: the physically-correct architecture (Sec. 0) DOES let the ratio
vary. Measured directly: FQ/FP ranges from ~0.85 to >2 across the gait cycle. The naive
whole-trial correlation with flexion angle is weak and not cleanly monotonic (Pearson r=0.123)
because it mixes stance and swing frames with **opposite-sign** sub-trends (stance1 r=−0.60, swing
r=+0.20, stance2 r=−0.53 — genuine hysteresis, not noise) — reported honestly as inconclusive on
its own, not oversold. The **fair, primary comparison** restricts to frames where BOTH tendons are
substantially engaged (>20 N each, n=110/150 frames, angle range 5.6–65.6°): median FQ/FP =
**1.152**. Binned by angle:

| angle (deg) | n | median FQ/FP |
|---:|---:|---:|
| 0–10 | 29 | 1.213 |
| 10–20 | 39 | 1.129 |
| 20–30 | 16 | 0.857 |
| 30–40 | 3 | 0.850 |
| 40–50 | 4 | 1.142 |
| 50–60 | 6 | 1.667 |
| 60–70 | 13 | 2.023 |

The **20–40° bin (0.850–0.857) is strikingly close to Buff 1988's own 10° value (0.86)** — and the
low-flexion (<20°) load-bearing median (1.152, the only regime level-walking's actual stance phase
populates) sits within a factor of 1.34× of Buff's 10° anchor, comfortably inside the pre-registered
2× band. **The naive 1:1 assumption is falsified** (ratio is never pinned at 1.0 in this model) and
the ratio's rough magnitude in the load-bearing regime is in genuine, non-tuned agreement with a
direct cadaveric measurement (Buff 1988) — but see Sec. 6 for why "quad force" here is a constant
baseline, not real dynamic activation, so this should be read as testing the model's own
**geometric** patella force-balance in isolation, not a full dynamic validation. **Honest scope
limit**: this gait trial's load-bearing regime never reaches Buff's 70° point at all — level
walking's weight-bearing stance phase stays under ~20° of flexion; only unloaded swing reaches
40–65°, where PT tension is small and the ratio metric is comparatively noisy (a near-zero-
denominator regime, the same "ratio of two small numbers is unstable" pattern already flagged
elsewhere in this repo, e.g. `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`'s CCI floor-ratio warning).

---

## 6. Forced discovery — the reused `gait_driven.h5`'s muscle activation was NEVER actually time-varying

While extracting the quad:PT ratio, `quad_force(t)` (sum of the 4 quadriceps muscles' `actuation`
output) came back **bit-identical at every one of 150 frames** (278.605 N, `np.ptp`< 1e-6). Forced
to Orient rather than accepted as a quirk: checked whether this equals
`0.05 × max_isometric_force` for each muscle, using lenhart2015.osim's own live XML values —

| muscle | Fmax (N) | 0.05×Fmax | h5-reported actuation (N) |
|---|---:|---:|---:|
| recfem_r | 848.8 | 42.440 | 42.440 |
| vasint_r | 1024.2 | 51.210 | 51.210 |
| vaslat_r | 2255.4 | 112.770 | 112.770 |
| vasmed_r | 1443.7 | 72.185 | 72.185 |

**Exact match, max abs error 1.4e-14 N** — and the same identity holds for every OTHER muscle
checked (both the 13 "activation-driven" ones and a spot-check of 5 "constant-control" ones:
`gasmed_r`, `tfl_r`, `soleus_r`, `addbrev_r`, `glmax1_r`, `iliacus_r`, `psoas_r` — ALL flat at
0.0500 activation, ALL exactly `0.05×Fmax`). **0.05 is OpenSim's own `Muscle`-class
`default_activation`.** Conclusion: in this pre-existing `gait_driven.h5` (built by a prior
session for `docs/MECHANISM_JAM_CONTACT_DECORR.md`, reused here read-only), **neither the real
per-frame SO activation in `actuator_input_file.sto` nor `forsim_settings.xml`'s own
`constant_muscle_control=0.02` fallback ever took effect** — every one of the 40 muscles ran at
its own static default the entire 1.57s trial. This appears specific to
`use_muscle_physiology=false` mode (Force = activation×Fmax, no equilibrium solve — the setting
`jam_contact_decorr.py` was forced into after two other numerical blockers, per that document's
own Sec. 4) silently not consuming any external control input at all.

**This is a genuine, previously-undisclosed limitation of already-existing, already-published
infrastructure, found as a side effect of this investigation — disclosed here, not fixed** (a fix
would require re-running the ~90s ForsimTool+JointMechanicsTool pipeline and is out of this
task's scope; flagged for whoever next touches `jam_contact_decorr.py`). Consequence for THIS
document: **Pathway B's absolute force magnitude (331.59 N) is very likely a floor/under-estimate**
— real gait quadriceps activation peaks well above the 0.05-activation baseline at key phases
(LaiArnold's own SO gives `recfem_r` tension up to 586.7 N alone at ITS gait peak,
`docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`, vs this run's flat 42.4 N for the same muscle) — and
the ratio finding (Sec. 5) is best read as testing the model's geometric patella mechanism in
isolation against a constant quad baseline, not under realistic dynamic muscle loading. This is
exactly why **Pathway A is treated as this document's primary number**: it is independently
verified (below) to use real, per-frame, time-varying SO tension with no analogous defect.

---

## 7. Machine cross-checks (13/13 PASS)

| check | result |
|---|---|
| Subject2 BW: fresh recompute from `MODEL_FILE` vs literal `SUBJECT2_BW_N` constant | PASS, rel err 0.00e+00 |
| Pathway A reproduces the established `knee_joint_contact` peak (391.11478 %BW) bit-exactly | PASS, err 0.00e+00 |
| Pathway A PF-peak survives edge-margin exclusion (not a differentiator artifact) | PASS |
| Pathway B settle window: force-based ≥ area-based (the forced-adversary fix, Sec. 3) | PASS (used sidx=9 vs area-only's 6) |
| LaiArnold structural ratio-unity (0 Ligament/PathSpring elements + path geometry) | PASS |
| Pathway B mesh-pair force agreement (patella-side vs femur-side, over-determination) | PASS, rel diff 1.28% |
| quad-force-flat diagnostic (0.05×Fmax identity, Sec. 6) | PASS (confirmed, not a pass/fail on quality — a disclosed finding) |
| Ratio falsifies naive 1:1 (load-bearing median 1.152 ≠ 1.0) | PASS |
| Ratio low-flexion magnitude within 2× of Buff 1988 @10° | PASS (1.152 vs 0.86, ratio 1.34×) |
| Both peaks finite (no NaN/Inf) | PASS (both) |
| Both peaks well below the 8000 N MVC ceiling (non-walking sanity bound) | PASS (both) |

Full detail: `data/msk_smoketest/patellofemoral_force/patellofemoral_force_results.json`
(`checks`, `overall_checks_pass: true`).

---

## 8. Anchors — live-verified this session (NCBI eutils + Europe PMC), not recalled

- **Task-given PMID 23267768 ("Kutzner et al. 2013 PFJ") is WRONG** — verified live via NCBI
  esummary: it is "Chromogenic substrate from 4-nitro-1-naphthol for hydrolytic enzyme..." (Dang J
  et al., *Bioorg Med Chem Lett*) — an unrelated chemistry paper. Same "wrong task-given PMID"
  pattern `docs/MECHANISM_KNEE_LIGAMENTS.md` already found and disclosed for its own two citations.
- **No genuine in-vivo instrumented PFJ force measurement exists in the literature** — searched
  live: `Kutzner[au] AND patella` → 1 hit, an unrelated case report by a different "Kutzner H"
  (not "Kutzner I", the OrthoLoad researcher); `instrumented AND patellofemoral AND "contact
  force"` → 17 hits, none an instrumented-PFJ-implant telemetry paper. Confirmed explicitly,
  independently, by **Trepczynski A et al. 2025** (PMID 39881960, *Front Bioeng Biotechnol*,
  Julius Wolff Institute — the SAME OrthoLoad-lineage lab as Kutzner's own tibiofemoral telemetry,
  fetched verbatim): *"knowledge about the in-vivo loading conditions at the PF joint remains
  limited, as no direct measurements are available."* That paper's own hybrid method (in-vivo
  instrumented-TKA tibiofemoral force + fluoroscopy kinematics → patient-specific-model PF force)
  gives peak PF force **1.75–3.29× BW** for sit-stand-sit/squat (high-flexion ADLs, NOT walking) —
  used only as directional/order-of-magnitude context, not a walking anchor.
- **Reilly DT, Martens M.** "Experimental analysis of the quadriceps muscle force and
  patello-femoral joint reaction force for various activities." *Acta Orthop Scand.*
  1972;43(2):126-37. **PMID 5079747** — verified to exist and be correctly attributed (title/
  journal/year live-confirmed via NCBI efetch + Europe PMC) but its abstract is **not indexed** in
  either database (pre-1975, common for this era) — the task's own "~0.5–1.5× BW walking" figure
  attributed to this paper could not be independently re-verified against primary text this
  session; used as the task's own pre-registered band, flagged as unconfirmed-primary-text.
- **Buff HU, Jones LC, Hungerford DS.** "Experimental determination of forces transmitted through
  the patello-femoral joint." *J Biomech.* 1988;21(1):17-23. **PMID 3339023** — fetched verbatim:
  *"The ratio between the tensions in the quadriceps tendon and the patellar tendon (FQ/FP) ranged
  from 1.55 at 70 degrees of flexion to 0.86 at 10 degrees of flexion... This study demonstrates
  that FQ does not equal FP as several authors have reported."* THE primary, quantitative ratio
  anchor (Sec. 5).
- **van Eijden TM, Kouwenhoven E, Verburg J, Weijs WA.** "A mathematical model of the
  patellofemoral joint." *J Biomech.* 1986;19(3):219-29. **PMID 3700434** — fetched verbatim,
  confirms this is exactly the task-named 1986 model (autopsy-knee-derived, computes PF
  compression + patellar-ligament force as a function of flexion angle) — the same conceptual
  architecture lenhart2015/JAM has (Sec. 0).
- **van Eijden TM, Weijs WA, Kouwenhoven E, Verburg J.** "Forces acting on the patella during
  maximal voluntary contraction of the quadriceps femoris muscle at different knee flexion/
  extension angles." *Acta Anat.* 1987;129(4):310-4. **PMID 3630619** — fetched verbatim:
  quadriceps force 2000 N (extension) to 8000 N (~75°), patellar ligament force max 5000 N (~60°),
  PF reaction force 1000 N (extension) rising to ≈ quad force at 75–90°. MVC-effort regime (NOT
  walking) — used only as the absolute-ceiling sanity bound (Sec. 1).

---

## 9. Honest gaps

1. **No in-vivo PFJ anchor exists at all** (Sec. 6, a verified negative, not a search failure) —
   confidence tier is cadaveric/published-plausibility, and the walking-band comparison rests on
   a literature figure (Reilly & Martens 1972) whose primary text could not be independently
   re-verified this session (pre-1975 abstract-indexing gap).
2. **Pathway B's absolute magnitude is likely a floor** (Sec. 6) — the discovered flat-muscle-
   activation artifact in the reused `gait_driven.h5` means quadriceps loading in that run never
   rises above a constant baseline; a corrected re-run (out of scope here) would likely give a
   higher peak PF force and a ratio curve tested under genuinely dynamic loading rather than a
   constant quad baseline.
3. **This trial spans slightly more than one gait cycle** (~1.36, per the knee-flexion trace) —
   Pathway A's own global PFJ peak (t=1.46s) falls in the SECOND captured loading-response window,
   not the first (which peaks only ~4% lower, 567.8 N at t=0.06s) — both are genuine, comparable
   events, not a single-cycle artifact, but this is a within-subject repeat, not two independent
   trials.
4. **Buff 1988's ratio curve is only tested at low flexion** (Sec. 5) — level walking's
   load-bearing stance phase never reaches the 70° point; the 40–65° comparison is confined to
   unloaded swing, where the ratio metric is intrinsically noisier (small-denominator regime).
5. **Single subject, single trial** — same scope limit as every other joint-force cert in this
   mission; no cross-subject PFJ sweep attempted here.
6. **LaiArnold's structural ratio-unity finding is a property of this one (common) modeling
   convention** — other continuous-quadriceps-path models (e.g. gait2392-family) likely share it;
   not verified against a third model in this session.

---

## 10. Files

- `scripts/msk/patellofemoral_force.py` — extraction + analysis + machine cross-checks for both
  pathways (self-contained given already-existing data; no new simulation).
- `data/msk_smoketest/patellofemoral_force/patellofemoral_force_results.json` — full machine-
  readable evidence (both pathways' peaks, checks, thresholds, anchors, verdicts).
- `data/msk_smoketest/patellofemoral_force/patellofemoral_force_series.json` — full per-frame
  series (PF force, knee angle, quad/PT tension, FQ/FP ratio) for both pathways.
- Reused, unedited: `scripts/msk/validate_joint_force.py`, `scripts/msk/bone_stress.py`,
  `scripts/msk/jam_contact_decorr.py`, `scripts/msk/jam_contact_decorr_analyze.py`.
- Reused, unedited data: `data/msk_smoketest/subject2_walking1/static_optimization/jr/
  walking1_JointReaction_ReactionLoads.sto`, `data/msk_smoketest/subject2_walking1/
  static_optimization/so/walking1_StaticOptimization_force.sto`,
  `data/msk_smoketest/jam_contact_decorr/results/joint-mechanics/gait_driven.h5`.
