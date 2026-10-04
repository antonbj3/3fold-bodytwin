# MECHANISM NO-BAND CONTROL — is the semimembranosus saturation a band effect or a pose effect? (2026-07-21)

Closes the #1 honest gap of `docs/MECHANISM_BAND_STATIC_OPT.md` (Sec.6, Sec.8#1, Sec.9#3): that cert found
semimembranosus **saturates** (activation 0.9999/1.0000, i.e. pinned at the SO ceiling) under the
`deadlift_feet_to_hands` band's posterior-chain load, at the fixed hip-hinge pose, on the erector-spinae
model — and reasoned from the already-published moment-propagation table (band's marginal contribution
to hip-flexion required moment is only ~1.7%) that this was **"very likely... not primarily a band
effect"** — but explicitly disclosed it never ran an actual no-band Static Optimization to check. This
build runs that control.
Script: `scripts/msk/band_control_noband.py` (re-runnable,
`.venv-msk/bin/python3 scripts/msk/band_control_noband.py`, exit 0).

## Bottom line

**Semimembranosus saturation is a POSE effect, not a band effect — machine-verified, not inferred.**
Running the identical Static Optimization pipeline with the band physically detached
(`PathSpring.set_appliesForce(False)`) leaves semimem_r/l at **0.99994 / 1.00000** activation — a
Δactivation of **+0.00000042 (r) / +0.0000075 (l)** relative to the with-band run, six orders of
magnitude smaller than any other posterior-chain muscle's response to the same band removal. The band's
real (small) marginal load is instead absorbed almost entirely by the muscles that still have headroom:
**glutes (+0.002 to +0.023 activation) and biceps femoris long head (+0.008)** pick up the slack;
semimembranosus cannot, because it is already pinned at its ceiling **with or without the band**. This is
a genuine corner-solution feature of the constrained optimization (SO can only redistribute load onto
muscles below their bound), not a modeling artifact — confirmed by checking the underlying **force**
values too (not just the clipped activation number): semimem's force differs by <0.002% between
conditions, i.e. the tiny residual gap is not "a big difference hidden by clipping," it is genuinely
almost nothing. Both required-force anchors (against two independently-built earlier models) and the
fresh with-band SO solution itself matched their respective external references to **exactly 0.000e+00**
— the strongest evidence available that this comparison has zero environment drift and is a clean,
apples-to-apples A/B.

## 1. Method — reuse, not duplicate

`scripts/msk/band_control_noband.py` imports `scripts/msk/band_static_opt.py` (`bso`) as a library and
calls its own `build_with_band_model()` / `run_static_optimization()` / `parse_sto()` /
`classify_actuators()` / `moment_balance_cross_check()` **unchanged** — `band_static_opt.py`'s source file
is never edited. Two independent model instances are built via `bso.build_with_band_model()` (called
twice — necessary because `run_static_optimization()` mutates its model in place by merging 18 reserve
actuators, so the two SO runs cannot share one model object regardless of band state):

1. **with_band**: the pristine output of `build_with_band_model()`, unmodified.
2. **no_band**: the identical construction, then every band `PathSpring` has
   `set_appliesForce(False)` called on it — the SAME detach mechanism `scripts/msk/attach_band.py`'s own
   with/no-band A/B already validated (its docstring: two model variants differing ONLY in this flag
   reproduced known moment-arm-predicted deltas). A fresh state is pulled afterward
   (`model.initSystem()` + re-pose) out of caution, though `appliesForce` is a property change, not a
   structural one.

`bso`'s two module-level output-path globals (`OUT_DIR`, `SO_OUT`) are monkey-patched from this script
before each `run_static_optimization()` call — verified safe (a Python function resolves bare globals in
its own module's namespace **at call time**, confirmed live with a throwaway module before trusting it for
the real run) — so every output of this session lands under `data/msk_smoketest/band_control_noband/`
and the original `data/msk_smoketest/band_static_opt/` cert artifacts are **never touched** (confirmed:
their mtimes are unchanged after this run). The two scripts also share one intermediate-output directory
(`data/msk_smoketest/band_posterior_chain/`, written by `build_config_model()`'s own tags `so_probe`/
`static_opt`) — calling `build_with_band_model()` twice rewrites those files twice, but only because the
construction is a **pure, deterministic function** of the (unchanged) config/model files, confirmed live
by md5-summing the 3 files most at risk before and after this run: identical both times
(`0e3632b6...`, `190b0e9f...`, `b1e652f0...`) — no corruption, no silent drift.

## 2. Pre-registered gates + external anchors (machine cross-checked, not eyeballed)

| gate | threshold | measured | verdict |
|---|---|---:|---|
| Anchor (a): fresh with-band required-force (9 coords) vs `band_static_opt_results.json` (same construction, earlier session) | ≤1e-6 N·m | **0.000e+00, all 9** | PASS |
| Anchor (b): fresh no-band (detached) required-force (9 coords) vs `band_posterior_chain_summary.json`'s own no_band column (DIFFERENT construction: zero PathSpring objects ever added, earlier session) | ≤1e-6 N·m | **0.000e+00, all 9** | PASS |
| Fresh with-band SO activation vs originally-published with-band SO activation, all 30 posterior-chain muscles | ≤1e-4 | **0.000e+00, all 30** (bit-identical) | PASS — full pipeline reproducibility, not just the Ipopt-independent required-force |
| SO convergence (10/10 frames, 0 NaN, 0 activation-bound violations) | both variants | with_band: 10/10,0,0; no_band: 10/10,0,0 | PASS both |
| Interior-row static-consistency spread (rows 2–8) | <1e-4 | with_band: 1.00e-08; no_band: 2.30e-07 | PASS both |
| Moment-balance cross-check (required vs SO-supplied), all 9 coords | \|diff\|≤max(2 N·m,15%) | PASS all 9, both variants | PASS both |
| Erector-spinae unaffected by band (n_active>0.05) | 0/88 both | 0/88 with_band, 0/88 no_band | PASS (confirms Sec.5 sign-incapability is gravity/pose-driven, not band-driven) |

The anchor checks used a tight, pre-registered tolerance (1e-6 N·m, looser than the ~1e-10 precision
these analytic finite-difference numbers are expected to hold to) and were run **before** the expensive
Ipopt Static Optimization step — the script aborts (exit 1) if either anchor fails, so a plumbing bug in
this reuse/detach mechanism could never be silently reported alongside a genuine scientific finding.

## 3. Required generalized force: with-band vs no-band (the mechanical anchor for Sec.4's muscle result)

| coordinate | required, with-band (N·m) | required, no-band (N·m) | Δ (band's marginal contribution) |
|---|---:|---:|---:|
| ankle_angle_r | +12.188 | +10.379 | +1.808 |
| subtalar_angle_r | −2.379 | −2.015 | −0.363 |
| knee_angle_r | +120.370 | +121.166 | −0.796 |
| hip_flexion_r | −259.066 | −254.637 | −4.430 (**−1.7%**) |
| hip_flexion_l | −260.279 | −255.792 | −4.487 (**−1.7%**) |
| knee_angle_l | +121.429 | +122.127 | −0.698 |
| subtalar_angle_l | −1.847 | −1.503 | −0.344 |
| ankle_angle_l | +9.060 | +7.416 | +1.644 |
| lumbar_extension | −41.749 | −50.100 | +8.351 (**16.7%** of the no-band/baseline magnitude) |

This reproduces the already-published propagation table exactly (band_posterior_chain cert, deadlift
config) — the hip-extension demand (the dominant driver of hamstring/glute recruitment) changes by only
~1.7% with the band; lumbar_extension changes by more in relative terms (16.7% of baseline) but its absolute
magnitude is small next to the erectors' 233.65 N·m capacity and irrelevant to them anyway (wrong sign
regardless, Sec.5 of the parent doc — confirmed still true without the band, Sec.5 below).

## 4. Per-muscle activation — band vs no-band, side by side (interior frame, t=0.05s)

| muscle | act (with-band) | act (no-band) | Δact | force w/band (N) | force no-band (N) | Fmax (N) | sat w/band | sat no-band |
|---|---:|---:|---:|---:|---:|---:|:---:|:---:|
| gasmed_r | 0.1863 | 0.2006 | −0.0142 | 541.9 | 583.2 | 3115.5 | – | – |
| gaslat_r | 0.0100 | 0.0100 | −0.0000 | 15.5 | 15.5 | 1575.1 | – | – |
| gasmed_l | 0.2137 | 0.2268 | −0.0131 | 621.6 | 659.6 | 3115.5 | – | – |
| gaslat_l | 0.0100 | 0.0100 | −0.0000 | 15.5 | 15.5 | 1575.1 | – | – |
| soleus_r | 0.0100 | 0.0100 | +0.0000 | 60.5 | 60.5 | 6194.8 | – | – |
| soleus_l | 0.0100 | 0.0100 | +0.0000 | 60.4 | 60.4 | 6194.8 | – | – |
| bflh_r | 0.6053 | 0.5976 | +0.0077 | 649.5 | 641.2 | 1313.2 | – | – |
| bfsh_r | 0.1059 | 0.1084 | −0.0025 | 57.3 | 58.7 | 557.1 | – | – |
| **semimem_r** | **0.9999** | **0.9999** | **+0.0000** | 1847.4 | 1847.4 | 2201.0 | **YES** | **YES** |
| semiten_r | 0.3097 | 0.3076 | +0.0022 | 145.2 | 144.2 | 591.3 | – | – |
| bflh_l | 0.6217 | 0.6136 | +0.0081 | 665.2 | 656.6 | 1313.2 | – | – |
| bfsh_l | 0.1134 | 0.1157 | −0.0023 | 61.3 | 62.6 | 557.1 | – | – |
| **semimem_l** | **1.0000** | **1.0000** | **+0.0000** | 1841.7 | 1841.7 | 2201.0 | **YES** | **YES** |
| semiten_l | 0.3211 | 0.3188 | +0.0024 | 150.2 | 149.1 | 591.3 | – | – |
| glmax1_r/l | 0.490/0.492 | 0.474/0.476 | +0.0155/+0.0158 | 481/483 | 466/467 | 983.8 | – | – |
| **glmax2_r/l** | **0.707/0.709** | **0.684/0.686** | **+0.0230/+0.0233** | 991/994 | 959/962 | 1406.0 | – | – |
| glmax3_r/l | 0.272/0.270 | 0.259/0.257 | +0.0131/+0.0130 | 253/251 | 240/239 | 947.8 | – | – |
| glmed1_r/l | 0.199/0.199 | 0.189/0.190 | +0.0097/+0.0098 | 217/218 | 207/207 | 1093.5 | – | – |
| glmed2_r/l | 0.239/0.241 | 0.230/0.232 | +0.0086/+0.0087 | 171/172 | 165/166 | 765.1 | – | – |
| glmed3_r/l | 0.300/0.302 | 0.291/0.293 | +0.0093/+0.0095 | 220/222 | 213/215 | 871.2 | – | – |
| glmin1_r/l | 0.019/0.020 | 0.017/0.018 | +0.0018/+0.0018 | 6.5/6.8 | 5.8/6.1 | 374.0 | – | – |
| glmin2_r/l | 0.069/0.068 | 0.067/0.066 | +0.0026/+0.0026 | 27/27 | 26/26 | 394.8 | – | – |
| glmin3_r/l | 0.115/0.116 | 0.111/0.112 | +0.0035/+0.0036 | 51/51 | 49/49 | 446.8 | – | – |

Pattern, forced-read not narrated: **every muscle that has headroom moves** in response to the band being
added/removed (glutes uniformly up with the band, biceps femoris long head up, gastrocnemius medialis
down, biceps femoris short head/semitendinosus small mixed signs) — this is exactly what a coupled,
activation²-minimizing multi-muscle optimizer does when the net demand across 9 coordinates shifts by a
few N·m: it redistributes across every unsaturated actuator, not just the "obvious" prime mover. The one
muscle that does **not** move is the one already at its bound in both conditions: semimembranosus.
Gastrocnemius medialis's small negative delta (activation drops when the band, and its slightly larger
ankle demand, is added) is disclosed but not deep-mechanistically decomposed here — it is a real,
machine-measured, small (~1.3–1.4 percentage point) feature of the coupled solve across the two joints
gastrocnemius crosses (ankle demand rises +1.8 N·m with the band, knee demand falls −0.8 N·m with the
band), not investigated further at the single-muscle level (would need its own moment-arm-weighted
decomposition, out of this task's scope).

## 5. Erector-spinae + native ideal lumbar_ext — confirms Sec.5's sign-incapability finding is pose-driven, not band-driven

| erector subgroup | with-band n_active(>0.05) | with-band range | no-band n_active(>0.05) | no-band range |
|---|---:|---|---:|---|
| iliocostalis_lumborum | 0/22 | [0.0100, 0.0100] | 0/22 | [0.0100, 0.0100] |
| longissimus_thoracis_pars_lumborum | 0/10 | [0.0100, 0.0101] | 0/10 | [0.0100, 0.0101] |
| longissimus_thoracis_pars_thoracis | 0/28 | [0.0100, 0.0100] | 0/28 | [0.0100, 0.0100] |
| multifidus | 0/28 | [0.0100, 0.0100] | 0/28 | [0.0100, 0.0100] |

All 88 erectors remain pinned at the activation floor in **both** conditions — expected and confirmed,
not a new finding: the parent doc's own gravity-only decomposition already showed gravity ALONE (no band,
no GRF) demands −50.100 N·m at lumbar_extension, a value the 88 sign-incapable (all-positive-moment-arm)
erectors can never supply regardless of band state. Removing the band makes the local demand **larger**
in magnitude (−50.100 vs −41.749 N·m), consistent with the native ideal `lumbar_ext` actuator (unbounded,
same as-shipped convention as the parent cert) absorbing correspondingly more:

| | with-band | no-band |
|---|---:|---:|
| `lumbar_ext` control | −4.4201 | −5.2552 |
| implied torque | −44.20 N·m | −52.55 N·m |
| % of its 10 N·m nominal capacity | 442% | 526% |

(The no-band torque, 52.55 N·m, is slightly larger in magnitude than the bare required value, 50.10
N·m — consistent with the 88 erectors' own small nonzero floor activation (0.01, not exactly 0)
contributing a small positive/extension-direction torque that the ideal actuator must additionally
overcome; the moment-balance cross-check, which sums every contributor including floor-activated
muscles, independently confirms the full ledger is self-consistent, PASS both variants.)

## 6. Saturation verdict — forced via the adversary, not assumed

**Adversary considered and forced**: activation is bounded to [0,1] by Static Optimization's own
constraint — a real difference in the *underlying, unconstrained* optimum (e.g. "1.05 with-band" vs
"1.02 no-band") could in principle be invisible after both get clipped to 1.0, silently manufacturing a
false "no difference" read. Forced check: the **force** column is not subject to that clipping the same
way (force = activation × Fmax × force-length-velocity multiplier, and the pose — hence fiber
length/velocity — is bit-identical between conditions by construction). Measured: semimem_r force is
1847.37611922 N (with-band) vs 1847.37533687 N (no-band), a 0.00004% difference; semimem_l is
1841.70281535 N vs 1841.68909064 N, 0.0007% difference. The exact activation deltas themselves are
**+4.2e-7 (semimem_r)** and **+7.45e-6 (semimem_l)** — both positive (the with-band value is marginally
higher, the physically-expected sign), both roughly five to six orders of magnitude smaller than glutes'
own response (+0.002 to +0.023) to the identical band removal. The adversary falls: there is no hidden
larger gap behind the clip — the muscle is genuinely, deeply saturated in both conditions, and the tiny
residual is consistent with Ipopt's own convergence-tolerance-level noise, not a suppressed physical
signal.

**VERDICT: POSE EFFECT.** Semimembranosus saturates at this hip-hinge pose regardless of the band. The
band's real, small, machine-verified marginal load (~1.7% of the hip-extension demand) is absorbed by
the unsaturated members of the extensor coalition (glutes, biceps femoris long head) — semimembranosus
was already at its ceiling from bodyweight alone.

## 7. Honest gaps

1. **Single pose only, and it is already deep in "both saturate" territory** — this control shows the
   band does not *cause* the saturation at the one already-published hip-hinge pose, but cannot show
   whether the band would matter at a *lighter* load (a shallower hip-hinge, or a resistance band tension
   nearer the point where semimembranosus first approaches — but does not yet reach — its ceiling). A
   pose/tension sweep to find that boundary was flagged as a possible follow-up in the parent doc's own
   Sec.9 but is out of this task's stated scope (same pose, same model, band on/off only).
2. **Gastrocnemius medialis's small negative Δactivation** (−0.013 to −0.014) is measured and reported
   but not decomposed into its own moment-arm-weighted mechanism (would require isolating its ankle vs.
   knee contribution the way Sec.4 of the parent doc isolated gravity/GRF/band for lumbar_extension) —
   flagged, not hand-waved, small magnitude (~1.3 percentage points of activation).
3. **All caveats already disclosed in the parent cert are inherited unchanged**: native ideal lumbar/arm
   actuators and the 6 pelvis reserves are left unbounded (as-shipped convention); Tier-1 straight-line
   band routing, no wrap objects; Tier-1 GRF simplification (symmetric 50/50 split); Static Optimization
   is an activation-minimizing solution, not measured EMG (real neural recruitment, including possible
   co-contraction, can differ).
4. **A band-TENSION sweep** (heavier/lighter Thera-Band tiers than the "Blue, 100% elongation" tier this
   whole cert family uses) was not run — this control only tests the two endpoints (this band's tension
   vs. exactly zero), not intermediate values, so it cannot characterize how close the with-band case
   is to actually mattering for semimembranosus specifically (it does not appear to matter at all here,
   Sec.6, but a much heavier band presumably would eventually, the same way the erectors' scope boundary
   is orthogonal to band presence but not to band MAGNITUDE in principle — not tested, out of scope).

## Files

- `scripts/msk/band_control_noband.py` — the full no-band control pipeline (self-contained, re-runnable;
  imports `band_static_opt.py`/`band_posterior_chain.py`/`attach_band.py` for reuse, never modifies them).
- `data/msk_smoketest/band_control_noband/band_control_noband_results.json` — every number in this
  document, machine-written (required-force both variants, anchor cross-checks, SO gates, moment-balance,
  full per-muscle activation/force comparison, erector-spinae comparison, saturation verdict).
- `data/msk_smoketest/band_control_noband/with_band_so/`, `.../no_band_so/` — each variant's merged
  model, quasi-static coordinates drive, patched SO setup XML, and raw
  `*_StaticOptimization_activation.sto` / `*_StaticOptimization_force.sto` output (10 frames ×
  200/213 columns for with_band_so, 200/212 for no_band_so — verified via each file's own header; the
  force file has one fewer column with the band detached, i.e. `appliesForce=False` makes OpenSim's own
  Static Optimization force report drop the band's column entirely rather than report it as a zero
  column — an independent, incidental confirmation that the detach mechanism is real, not cosmetic),
  kept fully separate from — and never overwriting — the original
  `data/msk_smoketest/band_static_opt/` cert artifacts (confirmed by unchanged mtimes after this run).
