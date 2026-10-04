# MECHANISM MSK — hip/ground-anchored band: does the posterior chain load from the ground up? (2026-07-21)

Follow-up to `docs/MECHANISM_MSK_ELASTIC_BAND.md`, which found a CHEST-anchored band loads the
lumbar directly but reaches hip/knee/ankle only INDIRECTLY (via ground-reaction-force whole-body
equilibrium), because legs and torso are PARALLEL branches off the pelvis. The operator asked for
band configurations that load the posterior chain **from the ground up** (ankle→knee→hip→lumbar),
on the erector-spinae-upgraded model. Script: `scripts/msk/band_posterior_chain.py` + hotswappable
dial `scripts/msk/band_posterior_chain_config.json`. Run:
`source_repository/.venv-msk/bin/python3 scripts/msk/band_posterior_chain.py`.
Reuses `scripts/msk/attach_band.py` (imported, not duplicated) + the real Thera-Band table in
`scripts/msk/band_config.json`, unchanged.

## Bottom line

**Yes — but only for 2 of the 3 configs tried, and exactly where the kinematic-tree geometry says
it should work.** The literal reading of "anchors low/front, pulls the pelvis into flexion" **fails
in exactly the same way the chest band failed** (zero direct hip/knee/ankle coupling) — because the
pelvis, like the torso, sits on the *proximal* side of the hip joint, not the distal side. Moving
the SAME anchor from the pelvis to the femur (thigh) gives genuine direct hip coupling. A band
looped under both feet up to both hands gives genuine DIRECT coupling at **every** coordinate in
the chain — ankle, subtalar, knee, hip, AND lumbar — confirmed two independent ways (geometric
moment arm, and a direct/indirect force decomposition), and its net effect is **100% direct** at
every one of the 9 measured coordinates (the band is a fully internal force element; ground-reaction
force is provably, and measurably, unchanged).

## 1. The geometric question, precisely

`docs/MECHANISM_MSK_ELASTIC_BAND.md` established: a `PathSpring`'s moment arm about coordinate q is
nonzero **iff** the path has one point on a body that is a *descendant* of q's joint in the
kinematic tree (serially distal) and one point that is not. Re-verified live on the erector-spinae
model (joint tree is unchanged from the original cert's model — confirmed by direct dump, not
assumed):

```
ground_pelvis: ground -> pelvis            hip_r: pelvis -> femur_r -> [walker_knee_r] -> tibia_r
back:          pelvis -> torso                     -> [ankle_r] -> talus_r -> [subtalar_r] -> calcn_r -> [mtp_r] -> toes_r
                                            hip_l: pelvis -> femur_l -> ... (mirror)
```

`hip_r`, `hip_l`, and `back` are **parallel branches off the pelvis** — pelvis and torso are both on
the *parent* side of every leg joint. This means: **anchoring a band anywhere on {ground, pelvis,
torso} — however "low" or "front" — cannot produce direct hip/knee/ankle coupling**, no matter how
the exercise is colloquially named. "Hip" in "hip-hinge" means the pelvis/waist region anatomically;
the **hip *joint*** requires an attachment on the **femur** side to get direct coupling. Only a path
touching {femur, tibia, talus, calcn, toes} — serially distal to hip/knee/ankle — can couple there
directly. This is the single fact that explains every result below.

## 2. Three configs tested (pre-registered predictions in the config JSON, before measuring)

| tag | anchor | pre-registered prediction |
|---|---|---|
| **pelvis_pullthrough** | ground-fixed low/front (0.45 m anterior of the heel) → **pelvis** (r.ASIS/L.ASIS midpoint) | literal reading of "(a) anchors low/front, pulls the pelvis into flexion." Predict: **zero** direct coupling everywhere in the chain, including lumbar (worse than the chest band). |
| **thigh_pullthrough** | same ground anchor (mirrored per side) → **femur** (proximal thigh markers) | forced-fix / geometric correction of (a): predict **nonzero** hip coupling only, zero knee/ankle/lumbar. |
| **deadlift_feet_to_hands** | **calcn_r** (r_calc marker) → **torso** (sternum midpoint) → **calcn_l** (L_calc marker), one continuous band, one tension | literal reading of "(b) band under both feet to the hands." Predict: **nonzero** coupling at every coordinate ankle→lumbar. |

Pose held **identical** across all 3 (and identical to the original cert's hip-hinge/RDL posture:
pelvis_tilt 25°, hip_flexion 45°, knee 20°, ankle 5°, lumbar_extension 5°) — isolates the anchor
variable, doesn't confound it with a posture change. Same real Thera-Band Blue calibration (secant
through the origin, exact at 100% elongation) as the original cert, applied independently per band
element.

**OODA note on config 2's design:** a first attempt chained the path as a single 3-point
`ground → femur_r → femur_l`. Measured live: hip moment arms of **−0.149 m vs −0.008 m** — grossly
asymmetric for a nominally symmetric pose, because `femur_r` was an unintended "pulley" node in that
routing, not real anatomy. Rejected; rebuilt as **two independent mirrored `PathSprings`** instead.
Re-measured: **−0.161 m vs −0.233 m** (the residual difference is real L/R anthropometric asymmetry
— same phenomenon the original cert noted for the heel markers — not a routing artifact).

## 3. Moment arms — the DIRECT, GRF-independent falsifier (machine-measured via `computeMomentArm`)

| coordinate | pelvis_pullthrough | thigh_pullthrough | deadlift_feet_to_hands |
|---|---:|---:|---:|
| ankle_angle_r | 0 (exact) | 0 (exact) | **−0.0735 m** |
| subtalar_angle_r | 0 (exact) | 0 (exact) | **+0.0148 m** |
| knee_angle_r | 0 (exact) | 0 (exact) | **+0.0323 m** |
| hip_flexion_r | 0 (exact) | **−0.0189 m** | **+0.1800 m** |
| hip_flexion_l | 0 (exact) | **−0.0451 m** | **+0.1823 m** |
| knee_angle_l | 0 (exact) | 0 (exact) | **+0.0284 m** |
| subtalar_angle_l | 0 (exact) | 0 (exact) | **+0.0140 m** |
| ankle_angle_l | 0 (exact) | 0 (exact) | **−0.0668 m** |
| lumbar_extension | 0 (exact) | 0 (exact) | **−0.3393 m** |
| pelvis_tilt (bonus) | −0.0258 m | −0.1191 m (both bands) | 0 (exact) |

"0 (exact)" = below 1e-6 m, the pre-registered noise floor (double-precision arithmetic on this
geometry is exact to ~1e-15; anything above 1e-6 is a real, physically meaningful coupling).
`deadlift_feet_to_hands` is the only config with **zero** empty entries in the target chain.

## 4. Force propagation — the virtual-work method (ID-Tool-bug-immune)

Same workaround as the original cert (`InverseDynamicsTool` corrupts knee/hip; see §5 for an
**expanded** finding this session). All numbers below are `dPE/dq` finite differences on gravity +
band elastic PE + GRF virtual work — never read from ID's own columns.

**Decomposition added this build:** `delta_direct_only` recomputes the with-band generalized force
holding GRF fixed at its **no-band** value — isolating the pure geometric/moment-arm channel from
the GRF-mediated whole-body-equilibrium channel. `delta_indirect = delta_total − delta_direct_only`.

### pelvis_pullthrough (tension 24.615 N)
| coordinate | no_band | with_band | Δtotal | Δdirect | Δindirect | Δtotal as % of baseline |
|---|---:|---:|---:|---:|---:|---:|
| ankle_angle_r | 10.379 | 11.398 | +1.019 | 0.000 | +1.019 | +9.8% |
| subtalar_angle_r | −2.015 | −2.245 | −0.230 | 0.000 | −0.230 | +11.4% |
| knee_angle_r | 121.166 | 118.949 | −2.217 | 0.000 | −2.217 | −1.8% |
| hip_flexion_r | −254.637 | −253.415 | +1.221 | 0.000 | +1.221 | −0.5% |
| hip_flexion_l | −255.792 | −254.660 | +1.132 | 0.000 | +1.132 | −0.4% |
| knee_angle_l | 122.127 | 119.999 | −2.128 | 0.000 | −2.128 | −1.7% |
| subtalar_angle_l | −1.503 | −1.723 | −0.220 | 0.000 | −0.220 | +14.7% |
| ankle_angle_l | 7.416 | 8.368 | +0.952 | 0.000 | +0.952 | +12.8% |
| lumbar_extension | −50.100 | −50.100 | **−0.0000** | 0.000 | −0.0000 | 0.0% |

100% of every nonzero entry is **indirect** (GRF-mediated) — identical failure mode to the original
chest band, and strictly worse (loses the lumbar coupling the chest band had).

### thigh_pullthrough (tension 24.615 N per band)
| coordinate | no_band | with_band | Δtotal | Δdirect | Δindirect | Δtotal as % of baseline |
|---|---:|---:|---:|---:|---:|---:|
| ankle_angle_r | 10.379 | 12.437 | +2.058 | 0.000 | +2.058 | +19.8% |
| subtalar_angle_r | −2.015 | −2.482 | −0.467 | 0.000 | −0.467 | +23.2% |
| knee_angle_r | 121.166 | 115.904 | −5.262 | 0.000 | −5.262 | −4.3% |
| **hip_flexion_r** | −254.637 | −250.075 | **+4.561** | **+0.464** | +4.097 | −1.8% |
| **hip_flexion_l** | −255.792 | −250.763 | **+5.028** | **+1.111** | +3.918 | −2.0% |
| knee_angle_l | 122.127 | 117.046 | −5.082 | 0.000 | −5.082 | −4.2% |
| subtalar_angle_l | −1.503 | −1.951 | −0.449 | 0.000 | −0.449 | +29.9% |
| ankle_angle_l | 7.416 | 9.350 | +1.934 | 0.000 | +1.934 | +26.1% |
| lumbar_extension | −50.100 | −50.100 | **−0.0000** | 0.000 | −0.0000 | 0.0% |

Only hip_flexion_r/l get a genuine direct contribution (10–22% of their own total Δ) — everything
else, including the bulk of hip's own change, is still indirect/GRF-mediated. Matches the moment-arm
prediction exactly.

### deadlift_feet_to_hands (tension 24.615 N, ONE continuous band)
| coordinate | no_band | with_band | Δtotal | Δdirect | Δindirect | Δtotal as % of baseline |
|---|---:|---:|---:|---:|---:|---:|
| **ankle_angle_r** | 10.379 | 12.188 | **+1.808** | **+1.808** | 0.000 | +17.4% |
| **subtalar_angle_r** | −2.015 | −2.378 | **−0.363** | **−0.363** | 0.000 | +18.0% |
| **knee_angle_r** | 121.166 | 120.370 | **−0.796** | **−0.796** | 0.000 | −0.7% |
| **hip_flexion_r** | −254.637 | −259.066 | **−4.430** | **−4.430** | 0.000 | +1.7% |
| **hip_flexion_l** | −255.792 | −260.279 | **−4.487** | **−4.487** | 0.000 | +1.8% |
| **knee_angle_l** | 122.127 | 121.429 | **−0.698** | **−0.698** | 0.000 | −0.6% |
| **subtalar_angle_l** | −1.503 | −1.847 | **−0.344** | **−0.344** | 0.000 | +22.9% |
| **ankle_angle_l** | 7.416 | 9.060 | **+1.644** | **+1.644** | 0.000 | +22.2% |
| **lumbar_extension** | −50.100 | −41.749 | **+8.351** | **+8.351** | 0.000 | −16.7% |

**Every single coordinate: Δindirect = 0.0000 exactly.** 100% of the propagated load at every joint
in the chain is DIRECT — not a GRF equilibrium side-effect. This is the mechanistically cleanest
result of the three.

## 5. Why config 3's GRF doesn't change at all (verified, not assumed)

`net_external_force()` sums the cable-tension force landing on real-body path points (excludes any
ground-fixed point, whose reaction is absorbed by the world). Machine result:

| config | contains a ground-fixed point | band net external force (N) | \|ΔGRF\| from no-band |
|---|---|---|---|
| pelvis_pullthrough | **True** | (21.39, −12.18, 0.08) | 12.31 N |
| thigh_pullthrough | **True** | (44.20, −21.31, 0.14) | 24.53 N |
| deadlift_feet_to_hands | **False** | (3.6e−15, −3.6e−15, 0) ≈ **0** | **0.000000 N** |

`deadlift_feet_to_hands` has both ends on real bodies (both feet + torso) — no fixed external
anchor. A taut multi-point cable under uniform tension exerts **zero net force on the whole body
system** by Newton's third law (a hand derivation, but not trusted blindly — verified numerically
via a general per-point cable-tension formula, confirmed to float-epsilon precision, ~1e-15 against
gravity_force ~767 N). This is *why* every joint's Δ is 100% direct there: there's no GRF channel
left for an "indirect" contribution to travel through. It also shows up independently in the pelvis
residual QC (§7): `pelvis_tilt` residual is **bit-identical** with vs without the band
(−605.518 both) for this config only — configs 1/2 (which DO add external load) shift it by 5.5–9.6 N·m.

## 6. Verdict: which anchoring best loads the full posterior chain

**`deadlift_feet_to_hands`** (band looped under both feet, up to both hands held at the chest).
It is the only config that:
1. Gets a nonzero, machine-verified moment arm at **every** coordinate in the target chain (§3).
2. Shows a nonzero force-propagation Δ at every one of those coordinates, and that Δ is **100%
   direct** (§4) — not contingent on the Tier-1 GRF-split simplification the way configs 1/2 (and
   the original chest band) are.
3. Is the only config whose mechanism generalizes cleanly to a heavier band tier or a dynamic trial
   without re-deriving a GRF assumption (§5) — it's load-path-complete on its own.

`thigh_pullthrough` is a legitimate, partial second place (genuine direct hip loading — real news
relative to the original chest-band cert, which had none) but leaves knee/ankle/lumbar entirely
GRF-mediated. `pelvis_pullthrough` — the most literal reading of "anchors low/front, pulls the
pelvis into flexion" — is a clean, forced, machine-verified **negative**: it fails exactly like the
chest band and is strictly worse (loses lumbar too), because "pulling the pelvis" never crosses the
hip joint itself. This negative was not accepted on a single pass — it was forced through the OODA
loop (design the fix, build it, remeasure) before being reported, per the build discipline; the fix
(`thigh_pullthrough`) is reported alongside it, not just the failure.

## 7. Two bugs found and fixed this session (forced via OODA, not assumed away)

**(a) A NEW numerical bug in the virtual-work method itself.** Repeated calls to
`Coordinate.setValue(state, value, enforceConstraints=True)` on a *reused* State are **not exactly
reversible** in this OpenSim version — a minimal, decisive reproduction: two back-to-back identical
calls (nothing perturbed in between) returned 10.248 then 10.362 N·m for the same quantity (drift
+0.113, no other coordinate touched). Localized via a 3-pass full-loop test: OpenSim's `Assembler`
(invoked internally by `enforceConstraints=True`) has a small, **history-dependent** residual on
every call — present even for coordinates with zero formal constraints. **Fix:** never invoke the
Assembler during finite-difference perturbation. Unconstrained coordinates use
`enforceConstraints=False` (verified bit-for-bit reproducible, spread=0.00e+00, across 3 full passes
through all 9 coordinates). The one real exception — `knee_angle_r/l`, which drive a dependent
`*_beta` coordinate via the patellofemoral `CoordinateCouplerConstraint` — is handled by reading the
constraint's own function directly (verified live: a `LinearFunction`, exactly the identity map,
checked at 0/20/45/90°) and setting the dependent coordinate manually, never via the Assembler.
Result: **byte-for-byte identical output across two independent full script runs** (`diff` clean).
This bug is orthogonal to (b) below and to the original cert's ID Tool bug; `attach_band.py`'s own
`generalized_force()` likely shares this latent issue (same pattern, same model family) — not
patched there (out of scope for this build; its reported numbers were already cross-validated to
<0.5–5% against ID's clean columns, well above this drift's scale) — flagged here for a future editor.

**(b) The original cert's ID Tool bug is BROADER than documented.** Re-verifying it live on this
model (not assumed inherited — a different `.osim` file) found `lumbar_extension_moment` **also**
corrupted at the hip-hinge pose (ID reports −2676.44 N·m vs a verified −50.10 N·m, ~53×) — the
original cert only checked lumbar at a neutral pose, where it happened to be clean (−7.59 vs
−7.61 N·m). Root-caused via a clean A/B, forced rather than assumed: the **original**
(pre-erector-spinae) model gives `lumbar_extension_moment = −7.611692` at the neutral pose (exact
match to the cert); **this** (erector-spinae) model gives `−2906.413542` at the *identical* pose,
with every other coordinate (pelvis_tilt, ankle, hip, knee) bit-identical between the two model
files. Disabling the 88 added erector-spinae muscles (`Force.set_appliesForce(False)`) restores
`lumbar_extension_moment` to exactly `−7.611692` — decisively localizing the cause to the mere
*presence* of the added muscles (even though `InverseDynamicsTool`'s documented algorithm is
supposed to be actuator-set-independent). Not traced to a specific C++ line (same evidence standard
as the original knee/hip diagnosis). **No impact on any number reported above** — lumbar was already
computed via the virtual-work method throughout — but it is essential, disclosed knowledge for any
future script that runs stock ID on this erector-spinae model and expects `lumbar_extension` to be
trustworthy: it is not, regardless of pose.

## 8. Honest gaps

1. **Tier-1 straight-line paths, no wrap objects** — same simplification as the original cert.
   `deadlift_feet_to_hands`'s 3.03 m path (calcn_r→torso→calcn_l) is a straight-line idealization;
   a real band would be redirected by body-surface contact along the way. Hotswappable: add a
   `WrapCylinder` later, same `PathSpring`.
2. **Single representative pose, not a measured/dynamic trial** — inherited from the original cert.
3. **GRF Tier-1 simplification** (symmetric 50/50 split, force-only equilibrium) still applies to
   configs 1/2, with the same large `pelvis_tilt` residual limitation the original cert disclosed
   (5.5–9.6 N·m·N shift here, vs −605 N·m baseline — not independently moment-balanced). Config 3
   sidesteps this specific gap entirely (§5) but inherits the *general* Tier-1 GRF simplification
   for the no-band baseline itself.
4. **Config 2's bilateral asymmetry** (−0.019 m vs −0.045 m hip moment arm) reflects real
   anthropometric L/R asymmetry in this subject's own markers, not a bug — same honest-asymmetry
   convention as the original cert.
5. **The pre-existing ID Tool bug (knee/hip, §4 of the original cert) plus the newly-found lumbar
   corruption (§7b here) are both diagnosed and worked around, not fixed upstream.** Any future
   script reading raw ID columns on this model family must not assume any specific coordinate is
   clean without re-checking live.
6. **Static Optimization / muscle activations not computed** — same descope as the original cert.
   Now more attractive as a next step: the model has 88 NEW erector-spinae muscles whose activation
   response to `deadlift_feet_to_hands`'s verified moments (§4) has never been computed.
7. **Baseline plausibility is a soft, normalized-literature-range check** (hip_flexion baseline
   3.256 N·m/kg bodyweight, within the ~1.5–3.5 N·m/kg ballpark cited in the original cert) — same
   pose, same soft anchor, not re-derived independently for this build.
8. **Anterior offset (0.45 m) for the ground-fixed anchor in configs 1/2 is a disclosed modeling
   choice**, not a measured trial — sized to plausibly sit at/beyond the toes, not derived from a
   specific cable-station geometry.

## 9. Exact next build step

1. **Static Optimization** using `deadlift_feet_to_hands`'s verified moments (§4) as the SO target —
   would show which of the 88 new erector-spinae muscles actually respond to this specific load
   path, the natural next question given the model was upgraded specifically for this.
2. **Add a `WrapCylinder`** around the pelvis/thighs for `deadlift_feet_to_hands`'s path (honest gap
   #1) to test whether the straight-line idealization materially changes the propagation table.
3. **If lumbar's ID corruption (§7b) matters beyond this script**: worth an upstream OpenSim issue
   report — a muscle-presence-triggered corruption of an unrelated coordinate's ID output is a
   general reliability concern for this model family, not specific to band-loading work.

## Files

- `scripts/msk/band_posterior_chain_config.json` — the hotswappable anchor dial (3 configs, each a
  list of body/marker points; change any value to re-parameterize, per the build's own convention).
- `scripts/msk/band_posterior_chain.py` — build + verification script (per-config calibration from
  the real Thera-Band table, moment arms, net-external-force cross-check, virtual-work propagation
  table with direct/indirect decomposition, JSON summary). Imports `attach_band.py` directly
  (`set_pose`/`marker_local`/`midpoint`/`vec3_*`/`read_sto_row`/`find_col`) rather than duplicating.
- `data/msk_smoketest/band_posterior_chain/model_{no_band,pelvis_pullthrough,thigh_pullthrough,deadlift_feet_to_hands}.osim`
  (+ `_probe` variants used for calibration) — the built models.
- `data/msk_smoketest/band_posterior_chain/coords_*.sto`, `id_*/` — ID Tool inputs/outputs (kept for
  the pelvis-residual QC; knee/hip/lumbar columns are the KNOWN-CORRUPTED ones, do not read them
  directly — use `band_posterior_chain_summary.json` instead).
- `data/msk_smoketest/band_posterior_chain/band_posterior_chain_summary.json` — machine-readable
  results (geometry, per-config calibration, moment arms, net-external-force, the full
  direct/indirect propagation table, pelvis residuals, ID-bug reproduction status).
- `data/msk_models/LaiArnoldModified2017_erector_spinae_subject2_scaled.osim` — the model used
  (gitignored; regenerate via `scripts/msk/add_erector_spinae.py` if absent).
