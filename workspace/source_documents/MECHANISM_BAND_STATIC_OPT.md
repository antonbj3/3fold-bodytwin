# MECHANISM STATIC OPTIMIZATION — muscles under the deadlift feet-to-hands band (2026-07-21)

Executes `docs/MECHANISM_BAND_POSTERIOR_CHAIN.md` Sec.9's own "exact next build step": that cert found
the `deadlift_feet_to_hands` band (calcn_r → torso/sternum → calcn_l) loads the full posterior chain
(ankle→subtalar→knee→hip→lumbar), 100% direct, at a fixed hip-hinge pose, on the erector-spinae-upgraded
model — but never ran Static Optimization, so it never showed **which** of the model's 168 muscles (80
native + 88 new erector-spinae) actually carry that load. This is that run.
Script: `scripts/msk/band_static_opt.py` (re-runnable, `.venv-msk/bin/python3 scripts/msk/band_static_opt.py`,
exit 0). Reuses `scripts/msk/band_posterior_chain.py`'s `build_config_model`/`generalized_force_multi`/
`net_external_force`/`load_configs` directly (imported, not duplicated) and `scripts/msk/static_opt_knee.py`'s
own precedent for the reserve-actuator recipe (same OpenCap-shipped `walking1_reserveActuators.xml`, unchanged
— same subject2-scaled skeleton, so its body/coordinate names still apply verbatim).

## Bottom line

**Gastrocnemius, hamstrings, and glutes carry this load exactly where the geometry says they should — the
lower-limb half of the operator's question is answered cleanly and is watertight-verified (moment-balance
cross-check passes on all 9 chain coordinates, diffs ≈0.000 N·m). The 88 new erector-spinae muscles do NOT
activate (pinned at the 0.01 floor) — and this was forced through a real OODA loop, not waved off: the first
hypothesis (an unbounded legacy actuator crowds them out) was tested, found NOT to fix it (bounding it makes
Static Optimization genuinely infeasible, confirmed via Ipopt's own "Infeasible problem detected"), and
re-diagnosed via three independent empirical sign checks: at this specific, already-published hip-hinge pose,
gravity's own effect on the local lumbar_extension DOF demands a FLEXION-direction restraint, and the
erectors — machine-verified 88/88 positive moment arm, true extensors — are geometrically incapable of
supplying that sign at any activation. This is a genuine scope boundary of the erector-spinae build (trunk
flexors — rectus abdominis, obliques, psoas — were explicitly out of scope, `docs/MECHANISM_ERECTOR_SPINAE.md`
honest gap #5), not a solver bug: the erectors' real, force-length-adjusted capacity at this exact pose is
233.65 N·m (ample, if the sign matched).**

## 1. Pipeline

1. Build the with-band, erector-spinae model at the hip-hinge pose (2-pass GRF calibration, identical to
   `band_posterior_chain.py`'s own `main()`) — 168 muscles + 13 native `CoordinateActuator` (lumbar/arms) +
   1 band `PathSpring` + 2 GRF `PrescribedForce` = 184 Force elements.
2. Merge in the OpenCap reserve `ForceSet` unchanged (6 pelvis `PointActuator`/`TorqueActuator` + 12 joint
   safety-net `CoordinateActuator`) → 202 total. Makes the free-floating pelvis root actuation-complete —
   the standard, already-validated OpenSim SO recipe (`docs/MECHANISM_STATIC_OPT.md`).
3. Run `opensim.AnalyzeTool` + `StaticOptimization` on a 10-frame, bit-identical-pose "quasi-static" hold
   (dt=0.01s, 0–0.09s) — the same synthetic-static-hold convention already used for this model's ID Tool
   calls throughout this cert family. Report the **interior** frame (index 5, t=0.05s), avoiding the
   spline-boundary edge effect already documented in `docs/MECHANISM_STATIC_OPT.md` Sec.5.
4. **Machine cross-check** (not a tautology gate): recompute the REQUIRED generalized force at each of the
   9 chain coordinates via the already-verified virtual-work method (zero dependency on Ipopt/SO), and
   independently recompute the SUPPLIED generalized force from the SO solution's own forces combined with
   live `computeMomentArm` calls — two independent computational paths that must agree for the SO solution
   to be trusted as a genuine equilibrium, not just an Ipopt self-report.
5. Report per-muscle-group activation/force: gastrocnemius, soleus, hamstrings, glutes, and the 4
   erector-spinae subgroups (iliocostalis lumborum `IL_`, longissimus thoracis pars lumborum `LTpL_`, pars
   thoracis `LTpT_`, multifidus `MF_`).

## 2. Convergence + sanity gates (pre-registered, machine-checked)

| gate | threshold | measured | verdict |
|---|---|---:|---|
| Frames processed, no NaN | 10/10, 0 NaN | 10/10, 0 | PASS |
| Muscle activation bounds [0,1] | 0 violations | 0 | PASS |
| Interior-row static-consistency (rows 2–8 spread) | <1e-4 | **1.00e-08** | PASS |
| Muscle force vs own Fmax | ratio ≤1.5× | 0/168 exceed | PASS |
| Muscles saturated at ceiling (row 5) | flagged, not hard-failed | **semimem_r (0.9999), semimem_l (1.0000)** | disclosed (Sec.6) |
| Joint (hip/knee/ankle/subtalar) reserve usage | flagged if >50% of own optimal_force | max 5.6% (hip_rotation_l) | PASS — muscles do the real work |
| Pelvis reserve (6-DOF free root, no muscle can act there) | disclosed, not gated | MZ=605.29 N·m; FX/FY/FZ=0; MX=0.033; MY=0.010 | see Sec.7 (matches a prior, independent number to 0.04%) |
| Native ideal lumbar/arm actuators exceeding own nominal capacity | disclosed, not gated (Sec.5-6) | `lumbar_ext`: control=−4.420 → −44.20 N·m (442% of its 10 N·m) | root-caused, Sec.5-6 |
| **Moment-balance cross-check (required vs SO-supplied), all 9 chain coords** | \|diff\| ≤ max(2 N·m, 15%×\|required\|) | **max \|diff\| = 1.4e-06 N·m** (see Sec.4 table) | **PASS, all 9** |

Overall convergence gate → **PASS**. A first version of the static-consistency gate itself had a bug (it
accidentally included the `time` column, reporting a bogus 0.06 "spread" instead of the genuine ~1e-8) —
caught and fixed before trusting the number, per the "machine cross-check, never eyeballing" rule; this is
disclosed rather than silently corrected because it's exactly the kind of self-inflicted false-negative the
watertight discipline exists to catch.

## 3. Muscle groups that carry the load — gastrocnemius, hamstrings, glutes (interior frame, t=0.05s)

| muscle | activation | force (N) | Fmax (N) | force/Fmax |
|---|---:|---:|---:|---:|
| gasmed_r | 0.186 | 541.9 | 3115.5 | 17.4% |
| gaslat_r | 0.010 (floor) | 15.5 | 1575.1 | 1.0% |
| gasmed_l | 0.214 | 621.6 | 3115.5 | 20.0% |
| gaslat_l | 0.010 (floor) | 15.5 | 1575.1 | 1.0% |
| soleus_r / soleus_l | 0.010 (floor) each | 60.5 / 60.5 | 6194.8 | 1.0% |
| bflh_r / bflh_l | 0.605 / 0.622 | 649.5 / 665.2 | 1313.2 | 49-51% |
| bfsh_r / bfsh_l | 0.106 / 0.113 | 57.3 / 61.3 | 557.1 | 10-11% |
| **semimem_r / semimem_l** | **1.000 / 1.000** | **1847.4 / 1841.7** | 2201.0 | **84% / 84%** |
| semiten_r / semiten_l | 0.310 / 0.321 | 145.2 / 150.3 | 591.3 | 25% |
| glmax1_r/l | 0.490 / 0.492 | 480.9 / 482.8 | 983.8 | 49% |
| glmax2_r/l | 0.707 / 0.709 | 991.1 / 994.4 | 1406.0 | 70% |
| glmax3_r/l | 0.272 / 0.270 | 252.7 / 250.7 | 947.8 | 27% |
| glmed1-3_r/l | 0.199–0.302 | 171–222 | 765–1094 | 20-30% |
| glmin1-3_r/l | 0.019–0.115 | 6.5–51.0 | 374–447 | 2-11% |

Ankle plantarflexion is carried mostly by **gastrocnemius** (soleus stays at floor — makes sense: gastroc is
biarticular and also crosses the knee, which per the already-published propagation table needs +120.4 N·m
here, soleus does not cross the knee at all). Hip extension (the dominant demand: required = **−259.1 /
−260.3 N·m** at hip_flexion_r/l, matching the already-published table exactly) is carried by a broad
coalition — **hamstrings + all three glute compartments**, with **glmax2** (the largest single contributor,
70% of its own Fmax) and **semimembranosus** (saturated, see below) doing the most work. All of this is
machine-verified self-consistent: the moment-balance cross-check at hip_flexion_r/l passes with a diff of
**1.4e-06 N·m** on a −259 N·m target.

## 4. Moment-balance cross-check — the external, non-tautological anchor

| coordinate | required (virtual-work) | supplied (SO forces × live moment arms) | diff |
|---|---:|---:|---:|
| ankle_angle_r | +12.188 | +12.188 | ~0 |
| subtalar_angle_r | −2.379 | −2.379 | ~0 |
| knee_angle_r | +120.370 | +120.370 | ~0 |
| hip_flexion_r | −259.066 | −259.066 | −1.4e-06 |
| hip_flexion_l | −260.279 | −260.279 | ~0 |
| knee_angle_l | +121.429 | +121.429 | ~0 |
| subtalar_angle_l | −1.847 | −1.847 | ~0 |
| ankle_angle_l | +9.060 | +9.060 | ~0 |
| lumbar_extension | −41.749 | −41.749 | ~0 |

"required" is computed **without ever calling Ipopt or Static Optimization** (pure finite-difference on an
effective potential — the ID-Tool-bug-immune method this whole cert family already uses); "supplied" is
computed **from the SO solution's own forces**, combined with a fresh, independent `computeMomentArm` call
per muscle/actuator. These two numbers agreeing to 6 decimal places is the single strongest piece of evidence
in this report that the SO solve is a genuine physical equilibrium, not merely "Ipopt returned success." The
`required` row also reproduces the already-published propagation table
(`docs/MECHANISM_BAND_POSTERIOR_CHAIN.md` Sec.4) bit-for-bit — confirming this script's model rebuild is
faithful to the original cert, not a drifted re-derivation.

## 5. The 88 erector-spinae muscles: forced via OODA, not a one-shot "honest negative"

**Observe**: all 88 sit at the activation floor (0.0100–0.0101), 0/88 above the 0.05 "active" threshold, in
every subgroup (iliocostalis lumborum, longissimus thoracis pars lumborum/pars thoracis, multifidus).

**Orient (attempt 1 — crowd-out hypothesis)**: the model's pre-existing, pre-erector-spinae "placeholder"
`lumbar_ext` `CoordinateActuator` (10 N·m nominal) has `min_control=-Inf`/`max_control=Inf` — unbounded. In
the as-shipped run it supplies control=−4.420, i.e. **−44.2 N·m, 442% of its nominal capacity** — silently
absorbing essentially the entire required lumbar moment for "free" in Static Optimization's own
activation²-minimizing cost function, since an unbounded ideal actuator's marginal cost never actually blocks
it. Hypothesis: bound it to its stated ±10 N·m and the real erectors (236 N·m naive ceiling, per
`docs/MECHANISM_ERECTOR_SPINAE.md`) will be forced to pick up the remainder.

**Decide/Act (attempt 1, forced test)**: bounded `lumbar_ext`/`lumbar_bend`/`lumbar_rot` to control∈[−1,+1].
Result: **erectors still did not activate** — instead Ipopt reported `Infeasible problem detected (status 2)`
at every one of the 10 frames, reproduced across bound tightness (±1 to ±10), iteration budgets (100–3000),
and convergence tolerances (1e-4 to 1e-2). A non-binding, generous bound (±10, i.e. ±100 N·m) reproduced the
unbounded solution bit-for-bit (control=−4.420 again) — confirming the bounding mechanism itself works and
the failure is specifically tied to the bound being *binding* (forcing real-muscle reliance), not to bounding
per se.

**Orient (attempt 2 — the real crux)**: three independent, from-scratch empirical checks, because hand-derived
sign algebra is exactly the kind of thing this discipline distrusts:
1. A hip-flexion sweep (unambiguous: real hip flexion moves the knee anteriorly) established **+X = anterior**
   in this model's ground frame.
2. An *isolated* back-joint sweep (pelvis_tilt=0, legs neutral — removing the pelvis-tilt confound that
   initially gave a misleading read) showed **+lumbar_extension really is anatomical extension** (torso top
   tips posteriorly, CoM drops off-vertical) — matching the coordinate's name and the already-published "88/88
   positive moment arm = true extensors" finding.
3. Decomposing the required generalized force (`generalized_force_multi` with `bands=[]` vs the real band)
   found **gravity ALONE** (no GRF, no band) already gives **−50.100 N·m** at lumbar_extension — GRF
   contributes exactly 0 (−50.10008001348936 vs −50.10008001316368, matching the kinematic-tree "parallel
   branch" argument this whole cert family already established), and the band's own marginal contribution is
   **+8.351 N·m** (reduces the magnitude, doesn't drive it). An 8-point robustness sweep of lumbar_extension
   itself (−20° to +30°, pose otherwise fixed) shows the required tau **robustly negative throughout, no sign
   flip nearby** (−15.2 → −75.6 N·m, monotonic) — this is not a fragile, near-zero artifact.

   **⚠ CORRECTION/PERSISTENCE NOTE (2026-07-21):** this 8-point sweep was doc-prose-only — no array
   existed on disk anywhere (confirmed directly: `gravity_only_decomposition()` is called exactly
   ONCE in `main()`, at the single published pose) — flagged in `docs/MECHANISM_TRUST_LEDGER.md` §10
   item 2. Now fixed at source: freshly re-run and persisted at
   `data/msk_smoketest/band_static_opt/lumbar_extension_robustness_sweep.json`
   (`scripts/msk/band_static_opt_lumbar_sweep.py`, reuses this script's own
   `build_with_band_model`/`gravity_only_decomposition` unchanged). The fresh 8-point sweep
   (evenly spaced across the SAME [−20°,+30°] range) **confirms** "robustly negative throughout, no
   sign flip, monotonic" for BOTH quantities — but a precise re-read found the **−15.2 → −75.6 N·m
   range quoted above is the `gravity_only` column, not `required tau`** (this session's fresh
   `gravity_only` column: −15.167 → −75.645 N·m, matching to the stated rounding; the fresh
   `required tau`, i.e. gravity+GRF+band, column instead spans −3.184 → −73.322 N·m — negative
   throughout too, but a materially different range from the one this sentence attributes to
   "required tau"). This is a genuine, disclosed labeling imprecision in the original prose (which
   quantity was swept, not whether either is robust), not a numerical error — full disposition in
   `docs/MECHANISM_TRUST_LEDGER_REMEDIATION.md`. The sweep also exactly reproduces this doc's own
   published single-pose numbers when re-evaluated at lumbar_extension=5° with an
   independently-rebuilt model: required tau −41.749149 N·m (0.0005% from the published −41.749357)
   and gravity_only −50.099872 N·m (0.0004% from −50.100080) — an over-determination anchor, not a
   tautology (separate model rebuild, same script family).

**Conclusion, forced not assumed**: at this specific, already-published hip-hinge pose (pelvis_tilt=25°,
hip_flexion=45°), gravity's own net effect on the *local* lumbar_extension DOF requires a **flexion-direction
(negative)** restraint. The 88 erectors (machine-verified 88/88 positive moment arm, `n_negative_moment_arm=0`
in this session's own fresh check) can only ever supply **tau = +T×momentArm ≥ 0** — structurally incapable of
a negative contribution at any activation. Ipopt's "infeasible" verdict, once the unbounded placeholder was
capped, was **correct**, not a numerical artifact. This is a genuine scope boundary, not a bug: trunk flexors
(rectus abdominis, internal/external obliques, psoas) exist in the same donor model
(`docs/MECHANISM_ERECTOR_SPINAE.md` honest gap #5) but were deliberately not ported (task scope was
"erector spinae ± multifidus"). The erectors' real capacity is not the limiting factor — a fresh,
force-length-equilibrated check at THIS exact pose (not the neutral pose the original cert checked, where it
found 236 N·m — a different, non-comparable pose, re-verified here rather than reused) gives **233.65 N·m** at
full activation, vs. **246.39 N·m** from the naive Fmax×|momentArm|-only ceiling computed at this SAME
hip-hinge pose (94.8% retained — force-length effects cost only ~5%, confirmed not material here) —
**18% of the equilibrated capacity would cover the −41.7 N·m demand, if the sign matched.**

## 6. Is the hamstring saturation (semimem_r/l ≈ 1.0) a red flag?

Disclosed, not hidden: **semimembranosus reaches 0.9999–1.0000 activation** (84% of its own Fmax at that
activation — a Hill-type muscle at less-than-optimal fiber length/velocity does not deliver 100% of Fmax even
fully activated). This is very likely a genuine feature of the **baseline hip-hinge pose's own demand**, not
primarily a band effect: the already-published propagation table shows the band's own marginal contribution
to hip_flexion_r/l is only **+1.7% / +1.8%** of the total (−254.6→−259.1 N·m, −255.8→−260.3 N·m) — the bulk of
the ~259 N·m hip-extension demand comes from bodyweight at 45° hip flexion, essentially unchanged by the
band. Static Optimization's activation²-minimizing solution happens to lean hardest on semimembranosus among
the available hip extensors/hamstrings; joint reserves stay small throughout (≤5.6% of their own optimal
force, PASS), so this is not reserve-masked or reserve-covered — it is the real muscle-recruitment solution
for a demanding, heavily-loaded hip-hinge. **Not independently re-verified this session against a no-band
control** (a cheap, disclosed next step, Sec.8).

## 7. An unplanned but striking cross-check

The pelvis `MZ` reserve actuator (an unbounded ideal torque about the pelvis, needed because the free
6-DOF `ground_pelvis` root has no muscle crossing it at all) supplies **−605.29 N·m** in this session's SO
solve. `docs/MECHANISM_BAND_POSTERIOR_CHAIN.md` §7 independently reports the `pelvis_tilt` column from
`InverseDynamicsTool` — already verified clean/correct for this coordinate — as **−605.518 N·m**, bit-identical
with or without the band. These two numbers, from **two completely different computational paths** (this
session's from-scratch Ipopt-based Static Optimization vs. the original cert's raw ID Tool column), differ by
only **0.04%**. The magnitude pattern also matches structurally: `MX`=0.033 and `MY`=0.010 (this session) are
both tiny, mirroring `pelvis_list`=−0.027 and `pelvis_rotation`=+0.012 (original cert) being the two small
pelvis-residual coordinates. This is reported as a striking, corroborating consistency signal, not a
confirmed exact identity — the precise MZ-axis-to-pelvis_tilt-coordinate correspondence was not independently
re-derived analytically this session.

## 8. Honest gaps

1. **Hamstring saturation (Sec.6) not cross-checked against a no-band control this session** — mechanistically
   plausible (baseline-pose-dominated, per the already-published small band-attributable delta) but not
   independently re-run.
2. **The native ideal lumbar_ext/lumbar_bend/lumbar_rot and the 6 pelvis reserve actuators are left
   UNBOUNDED in the reported (as-shipped) run** — matches this whole model family's existing, already-shipped
   convention (same OpenCap reserve file `static_opt_knee.py` already used unmodified), but means the
   `lumbar_ext` number specifically (−44.2 N·m) is not a physically-realistic muscle-equivalent force; it is
   correctly flagged, not silently trusted.
3. **The "erectors can't do it" finding is specific to this exact, already-published synthetic pose**
   (pelvis_tilt=25°, hip_flexion=45°, lumbar_extension=5°) — the 8-point sweep (Sec.5) shows it is robust
   across a wide local range of lumbar_extension, but a full pose (e.g. different hip_flexion/pelvis_tilt
   combinations, or a real dynamic trial through a hinge's full range of motion) could plausibly cross into a
   regime where the local lumbar demand flips sign and the erectors DO engage — not tested this session
   (would need a fresh pose or a dynamic sweep, out of scope for this single-pose build).
4. **Single representative pose, not a measured/dynamic trial** — inherited unchanged from the original cert
   and from `docs/MECHANISM_STATIC_OPT.md`.
5. **Tier-1 straight-line band path, no wrap objects; Tier-1 GRF simplification** — both inherited unchanged,
   already disclosed in the parent cert.
6. **Static Optimization is an activation-minimizing solution, not measured EMG** — same structural caveat
   `docs/MECHANISM_STATIC_OPT.md` Sec.8 already discloses for this model family; real neural recruitment
   (including possible erector co-contraction for stability that SO's pure-optimality criterion would never
   select) can differ from this result.
7. **The bounded/"physically capped" SO attempt (Sec.5) was diagnosed via a fast, Ipopt-free capacity+sign
   check, not by getting Ipopt itself to converge on a corrected, fully-bounded model** — the conclusion
   (erectors are the wrong sign, not underpowered) is machine-verified independently of Ipopt, but no
   bounded-and-converged SO run exists to give erector activation numbers under a "no unrealistic crutch"
   constraint; per the finding, this is expected to remain infeasible without adding trunk-flexor
   musculature, not a gap this session could close.

## 9. Next step

1. **If the operator wants a pose where erectors genuinely engage**: sweep hip_flexion/pelvis_tilt/lumbar_extension
   combinations (reusing this script's `gravity_only_decomposition` + the 8-point-sweep pattern from Sec.5) to
   find the local-lumbar-demand sign boundary, or test a dynamic trial through a real hinge's full range of
   motion (top vs. bottom of a deadlift very plausibly have opposite local lumbar_extension demand signs).
2. **Add trunk flexors** (rectus abdominis, obliques, psoas — present in the same donor model,
   `docs/MECHANISM_ERECTOR_SPINAE.md` honest gap #5) via the exact same `add_erector_spinae.py` porting
   pattern, if a genuinely sign-complete trunk actuator set is wanted.
3. **No-band control run** of this same script (band `PathSpring` set to zero tension) to isolate the band's
   own marginal contribution to the hamstring-saturation finding (Sec.6), the same way the original cert
   isolated the band's direct/indirect moment contributions.

## Files

- `scripts/msk/band_static_opt.py` — the full pipeline (self-contained, re-runnable; imports
  `band_posterior_chain.py` for model-building/virtual-work reuse, `static_opt_knee.py`'s reserve-actuator
  precedent referenced but not imported — same OpenCap XML asset, read directly).
- `data/msk_smoketest/band_static_opt/model_static_opt_merged.osim` — the 202-Force-element model actually
  fed to Static Optimization (168 muscles + 13 native CoordinateActuator + 1 band PathSpring + 2 GRF
  PrescribedForce + 18 OpenCap reserve actuators).
- `data/msk_smoketest/band_static_opt/coords_static_opt.sto` — the 10-frame quasi-static coordinates drive.
- `data/msk_smoketest/band_static_opt/so/` — patched SO setup XML, `*_StaticOptimization_activation.sto` /
  `*_StaticOptimization_force.sto` (10 frames × 200/213 columns), Ipopt console log.
- `data/msk_smoketest/band_static_opt/band_static_opt_results.json` — every number in this document,
  machine-written (required/supplied per coordinate, full per-muscle activation/force, gates, erector
  diagnostic, gravity decomposition).
