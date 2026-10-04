# MECHANISM WOBBLING-MASS — soft-tissue layer for shank/thigh, first outer-layer step (2026-07-21)

Executes the operator's ask: a first-step "wobbling mass" (soft-tissue) layer, toward the outer
layers ("allt under huden" + skin) — split a segment's existing rigid mass into a bone core +
spring-damper-coupled soft-tissue mass, for shank/thigh, driven by REAL subject2 mocap, quantify
the effect on peak/impact segment force vs the rigid assumption, and connect it explicitly to the
skin-marker soft-tissue-artifact (STA) this repo's own audit already flagged as a live,
silently-uncorrected error source in every IK number it produces
(`docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md` row "Dermis/epidermis"). Script:
`scripts/msk/wobbling_mass.py`. Run: `.venv-msk/bin/python3 scripts/msk/wobbling_mass.py`
(no args; writes `data/msk_smoketest/wobbling_mass/wobbling_mass_results.json`).

## Status: MECHANISM VERIFIED (every internal/analytic check passes exactly). REAL-DATA
## APPLICATION IS PARTIAL, NOT A CLEAN WIN — disclosed candidly, not oversold. The 100Hz-IK
## real-data path is resolution-limited (own validity gate FAILS for the drop-landing trial); the
## strongest real evidence comes from a supplementary, real-measured, higher-rate driver, and it
## is parameter-sensitive in a way this session DERIVED and CONFIRMED the mechanism for (the
## `crossover_r(mu)` geometric result below), not merely observed as noise.

---

## 1. Model & geometry (derived, not assumed)

Per segment, the EXISTING rigid-body mass `M` (unchanged, model-measured — no new mass invented)
splits into a bone core `m_r` and a soft-tissue "wobbling" mass `m_w` (`M = m_r + m_w`, conserved,
checked in code), coupled by a spring `k` + damper `c`, one relative translational DOF `z(t)` along
the segment-vertical (dominant impact) axis. This is Gruber (1998) / Pain & Challis's own
"wobbling mass model" class, at lumped 1-segment/1-relative-DOF granularity.

Base-excitation EOM: `m_w·z'' + c·z' + k·z = -m_w·a(t)`, where `a(t)` is the bone/skeleton
acceleration (what "rigid" already assumes is the whole segment's acceleration).

```
F_rigid(t)  = (m_r+m_w)·a(t)                      -- the existing rigid-body assumption
F_wobble(t) = m_r·a(t) - (k·z(t) + c·z'(t))        -- the true bone-core reaction force
EXACT IDENTITY:  F_wobble(t) = F_rigid(t) + m_w·z''(t)     (proved once, asserted every run)
```

Three analytic, non-tautological limits (geometric consequences of the same EOM, not hand-tuned):

| limit | prediction | measured (synthetic battery) | tol | verdict |
|---|---|---|---|---|
| `k→∞` (weld) | ratio→1 | 1.0006 | 1% | PASS |
| `m_w→0` (no wobbling mass) | ratio→1 | 1.0000 | 1% | PASS |
| sinusoidal `ω→∞` | ratio→`m_r/M = 1-μ` | 0.3012 (want 0.3000, μ=0.7) | 1% | PASS |

A fourth, new geometric result this session **derived and then verified** (not assumed): solving
`|F_wobble/F_rigid| = 1` for the frequency ratio `r = ω/ω_n` gives

```
r*(μ) = sqrt(2 / (2-μ))          -- INDEPENDENT of ζ (verified: 9/9 (μ,ζ) combos, ratio(r*)=1.000000
                                      to 1e-6, for μ∈{0.5,0.7,0.9}×ζ∈{0.3,0.5,1.0})
```

Below `r*(μ)`: the wobbling mass **amplifies** peak segment force (ratio>1). Above `r*(μ)`: it
**attenuates** (ratio<1). This is the μ-dependent analogue of the classical vibration-isolation
"√2 crossover" (a different transfer function than the textbook one, since `F_wobble/F_rigid` is
not the classical displacement/force transmissibility — re-derived here for this specific
quantity, then checked). This single formula is what later **explains** (not just discloses) why
the real-data parameter sweep below shows a mix of positive and negative "attenuation."

**Adversary forced (OODA) before trusting any real number**: "this is a numerical-integration
artifact, not real physics." Forced via a synthetic battery run before any real data is touched:
(1) the exact identity above, machine-asserted to `1.27e-16` relative error; (2) the time-domain
RK4 integrator cross-checked against the independent closed-form frequency-response `H(ω)`
(over-determination — two different derivations of the same number) across 6 frequencies
(3–45 Hz), all agreeing to ≤0.02%; (3) the three void-floor limits above; (4) the `crossover_r(μ)`
geometry. **All pass — the battery gates the script itself: if any of these fail, `main()` aborts
before touching real data (measured: it did abort on a first draft, see §4).**

## 2. Literature anchor (live-verified this session — PMID/DOI resolved via EuropePMC/PubMed, not recalled)

| ref | PMID | DOI | finding used |
|---|---|---|---|
| Gruber K, Ruder H, Denoth J, Schneider K. *A comparative study of impact dynamics: wobbling mass model versus rigid body models.* J Biomech. 1998;31(5):439-44. | 9727341 | 10.1016/S0021-9290(98)00033-5 | Qualitative, founding result: rigid-body inverse dynamics gives "completely incorrect internal torques and forces" during impact vs a wobbling-mass ground truth. No numeric table in the abstract. |
| Pain MT, Challis JH. *The role of the heel pad and shank soft tissue during impacts: a further resolution of a paradox.* J Biomech. 2001;34(3):327-33. | 11182123 | 10.1016/S0021-9290(00)00199-8 | Qualitative: the shank's wobbling mass makes a "significant" contribution to impact energy dissipation. |
| Pain MT, Challis JH. *The influence of soft tissue movement on ground reaction forces, joint torques and joint reaction forces in drop landings.* J Biomech. 2006;39(1):119-24 (epub 2004-12-30). | 16271595 | 10.1016/j.jbiomech.2004.10.036 | **THE quantitative anchor.** 4-segment wobbling-mass model, 43cm drop landing: peak-vGRF model error 1.2%; joint forces/torques "up to nearly 50% lower" for wobbling vs rigid. Used as an order-of-magnitude plausibility anchor, never a tautology gate. |
| Peters A, Galna B, Sangeux M, Morris M, Baker R. *Quantification of soft tissue artifact in lower limb human motion analysis: a systematic review.* Gait Posture. 2010;31(1):1-8. | 19853455 | 10.1016/j.gaitpost.2009.09.004 | STA-magnitude review, cited as the connective pointer (§6). Abstract-only access did not yield an extractable mm figure this session (paywalled body text) — an honest gap, not papered over. |

A live search for the other commonly-cited numeric wobbling-mass parameter source (Liu & Nigg,
~2000) did not resolve a verifiable record this session — not used, not guessed from recall.

## 3. Data (real subject2 lab mocap, not synthetic — reused, not re-derived)

Same subject, same model family `validate_joint_force.py`'s own PASS-gated joint-force cert
already uses (`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`), same
`parse_mot`/`savgol_smooth_and_derivs` pipeline (reused directly, not reimplemented).

Real, live-queried model masses (not a generic anthropometry table):

| body | mass (kg) |
|---|---|
| femur_r | 9.65488 |
| tibia_r | 3.84839 |
| patella_r | 0.08948 |
| **total body** | **78.2000** (matches `MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md`'s independently-stated 78.2 kg) |

Segment convention (exact match to `force_transmission_scene.py:547-548`): **thigh_r** =
femur_r+patella_r = 9.744 kg; **shank_r** = tibia_r = 3.848 kg.

Two real trials, both with real 2000Hz force-plate GRF:

- **walking1** (gait): 100Hz IK, 158 samples, 1.57s.
- **DJ1** (real drop-landing): 100Hz IK, 100 samples spanning [1.25,2.24]s; real measured peak
  vertical GRF 465%BW at contact — the SAME activity class (drop landing) as the Pain & Challis
  (2006) anchor above.

## 4. Two bugs found and fixed via forced OODA (not hidden — this is the actual work)

**Bug #1 — activity-mismatched onset detector.** First version detected "contact onset" as a
rising edge crossing 5%BW from a near-zero baseline: correct for DJ1 (real flight phase, GRF
genuinely ≈0 pre-contact) but **measured to fail-open on walking1** — its own real vy_total never
drops below 83.6%BW anywhere in the 1.57s trial (both legs continuously share body-weight support;
no flight phase in gait). Fix: detect the instant of steepest loading rate (max `d(vy_total)/dt`,
biomechanics' own "rate of force development" peak) — one detector for both activity classes.

**Bug #2 — convolution edge artifact, caught by re-running the fix on real data before trusting
it.** The smoothed-then-differentiated onset detector's `mode="same"` convolution zero-pads at the
array boundary; on walking1 this produced a spurious apparent spike at t=0.0015s (RFD=143,720 N/s)
that the naive argmax picked up instead of the real event. Verified directly: the RAW (unsmoothed)
signal's own steepest rise is at **t=0.579s** (RFD=20,287 N/s, ~7× smaller — the smoothed edge
spike was numerical, not physical), and inspecting the raw GRF there confirms a genuine loading-
response rise (107%BW trough at t=0.55s → 124%BW peak by t=0.69s, the classic double-hump walking
vGRF). Fix: exclude the smoothing kernel's own edge margin from the argmax search.

**A third, non-bug finding, diagnosed not just observed**: the OTHER validated cert in this repo's
own SavGol window (11 samples = 110ms, tuned for low-frequency muscle-driven joint kinetics) is
comparable to or LONGER than the period of every frequency this wobbling-mass question sweeps
(8-35Hz → 125-29ms) — i.e. it structurally smooths away exactly the content this model needs. Fix:
a dedicated window-sensitivity sweep (5,7,9,11 samples, uniform poly-order=3 for fair comparison),
smallest-that-is-computed taken as primary (same "smallest window that's usable" convention
`force_transmission_scene.py`/`validate_joint_force.py` already established for a different
question) — reported as a full sensitivity table, not silently picked.

## 5. Results (representative point: μ=0.7, f_n=15Hz, ζ=0.5 — one generic midpoint, same for both
## segments; full sweep μ∈[0.5,0.9]×f_n∈[8,35]Hz×ζ∈[0.3,1.0], 100 combos, reported alongside)

| trial | segment | validity gate (whole-body residual, ≤8%BW) | attenuation % @ representative | sweep min/median/max % | timing shift | crossover sign-match |
|---|---|---|---|---|---|---|
| walking1 | thigh_r | **PASS** (1.4–4.7%BW) | −3.58% | −18.0/−1.1/+3.3 | 0.0ms | 62% |
| walking1 | shank_r | **PASS** | +3.55% | −2.7/+2.5/+6.8 | 0.0ms | 9% |
| DJ1 | thigh_r | **FAIL** (14.6–15.6%BW x, 18.7–28.2%BW y, all windows) | +10.04% | −20.9/−2.6/+14.5 | +10.5ms | 84% |
| DJ1 | shank_r | **FAIL** | +9.37% | −17.0/−2.1/+14.3 | +10.0ms | 82% |

**Supplementary** (real 2000Hz measured-GRF whole-body proxy, first ~40ms of impact — sidesteps the
100Hz validity concern entirely by never using IK-derived kinematics, at the cost of a different,
disclosed approximation: whole-body acceleration standing in for each segment's own):

| trial | attenuation % @ representative | sweep median |
|---|---|---|
| walking1 | −4.50% | −0.10% |
| **DJ1** | **+19.79%** | −4.40% |

**Reading this honestly**: DJ1's 100Hz-IK PRIMARY numbers (+10.0%/+9.4%, both "look" literature-
consistent) are NOT counted as a confirmed claim — their own whole-body Newton-residual validity
gate fails at every window tried, 2-3.5× over the pre-registered 8%BW ceiling (reused from
`validate_joint_force.py`), elevated in **all three axes** (not just vertical — consistent with a
broad differentiation-resolution ceiling for this faster/larger motion, not an isolated bug).
Trusting a number from a signal that fails its own validity check is exactly the "unforced
adversary" this repo's method exists to catch. walking1's PRIMARY passes validity cleanly but its
own detected event (a mid-stance-trough-to-loading-response transition, peak accel only
1.6-2.4 m/s² — genuinely gentle, not a sharp heel-strike-from-swing) gives a weak, sub-floor,
sign-unstable-across-windows signal — an honest negative for THIS gentle event.

The most trustworthy real-data number in the whole analysis is the **DJ1 supplementary +19.8%**:
real, 2000Hz, directly measured (no IK-derived differentiation at all), same direction and same
order of magnitude as Pain & Challis (2006)'s independently-published "up to nearly 50% lower"
(smaller — expected: different subject, different drop height, generic swept parameters never fit
to reproduce 50%). Its full-sweep median is negative (−4.4%) — **mechanistically explained, not
hand-waved**: this window's own dominant frequency (~10Hz, machine-measured via FFT) sits close to
the `crossover_r(μ)` boundary for several swept `f_n` values, so roughly half the generic grid
predicts amplification and half attenuation — confirmed by the crossover-sign-prediction check
matching the ACTUAL sign of 82-84% of DJ1's 100 swept combinations (vs the much weaker 9-62% match
on walking1's shorter, gentler window, where a single-dominant-frequency proxy is itself a poor
description of a 15-sample, non-impulsive signal — disclosed, not hidden).

## 6. Claim / adversary / anchor (watertight form)

**Claim A — mechanism**: a spring-damper-coupled wobbling mass measurably attenuates-or-amplifies
peak segment force relative to the rigid assumption, with sign/magnitude governed by the geometric
crossover `r` vs `r*(μ)`, not by chance. **Adversary**: numerical-integration artifact. **Forced**:
exact identity (1.27e-16), closed-form cross-check (≤0.02%, 6 frequencies), 3 void-floor limits
(≤1%), `crossover_r(μ)` geometry (9/9 exact). **Anchor**: Gruber 1998 (qualitative, same paradigm),
Pain & Challis 2006 (quantitative, ~50%-lower direction/magnitude). **Verdict: C — CONFIRMED.**

**Claim B — present in THIS repo's real subject2 data for shank/thigh**: **Adversary**: the
measured "attenuation" is a 100Hz-IK differentiation/smoothing artifact, not real physics.
**Forced**: whole-body Newton-residual validity gate (reused, not invented, ≤8%BW) swept across 4
SavGol windows, for both trials. **Result**: the adversary does NOT fall for DJ1's primary path
(gate fails at every window) — that evidence is honestly excluded, not counted. It DOES fall for
walking1 (gate passes) but the resulting signal is sub-floor and sign-unstable — inconclusive, not
a confirm. A structurally different, real, higher-rate path (measured-GRF-driven, immune to this
specific adversary) gives a real, sizable, literature-direction-consistent number for DJ1.
**Verdict: C — PARTIALLY CONFIRMED, first-step, resolution-limited** (not NOT-CONFIRMED: the
mechanism is real and the best-available real evidence supports it; also not a clean CONFIRM: the
segment-specific 100Hz-native path itself does not yet clear its own validity bar).

## 7. STA connection (the mocap-artifact side of the same physics)

`z(t)` — the SAME relative-displacement DOF that filters/attenuates the bone-core force — is, by
construction, this model's own first-step PREDICTION of how far a skin marker on that segment
moves relative to the true bone during the same impact: exactly the soft-tissue-artifact this
repo's own audit already flagged as unmodeled and "silently-uncorrected" in every IK number it
produces (`docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md`, "Dermis/epidermis" row). At DJ1's
representative point: thigh_r `z` peak=12.7mm/rms=6.0mm at ~10.5Hz; shank_r peak=25.2mm/rms=12.3mm
at ~5.2Hz — magnitudes and frequencies squarely in the range that would corrupt marker-derived
joint angles for these segments during landing. This is a MECHANISTIC connection (the same
mathematical object drives both effects), not a numerically cross-anchored one: Peters et al.
(2010, PMID 19853455) is the standard STA-magnitude review, but abstract-only access this session
did not yield an extractable independent mm figure to check these numbers against — an honest gap
for a follow-up session, not a claimed validation.

## 8. Honest gaps (first step — see the script's own module docstring for the full list)

1. Lumped 2-mass-per-segment, 1 relative DOF — not FEM soft tissue; muscle/fat/skin are one
   combined "wobbling mass," not separated.
2. μ, f_n, ζ are generic and swept (not a live-verified precise per-segment table — paywalled body
   text of the cited papers), turned into a robustness envelope rather than an invented point value.
3. 100Hz IK materially limits this analysis — measured, not just anticipated (§4-5).
4. The supplementary measured-driver applies whole-body acceleration as a per-segment proxy, valid
   only for the initial ~40ms.
5. STA connection is mechanistic, not numerically cross-anchored to an independent measurement.
6. Only the right-side thigh/shank of two real trials (one gait, one drop-landing) — no
   left/right, no additional subjects, no jump take-off (only landing).

## 9. Files

- `scripts/msk/wobbling_mass.py` — the full model, synthetic adversary battery, and real-data
  pipeline (new file, ~900 lines).
- `data/msk_smoketest/wobbling_mass/wobbling_mass_results.json` — full machine-readable output
  (pre-registered thresholds, literature, battery, both trials' segment + supplementary results,
  headline verdict table).
- Reused, not re-derived: `scripts/msk/validate_joint_force.py` (`parse_mot`,
  `savgol_smooth_and_derivs`, `G`, `MODEL_FILE`, `SUBJECT_DIR`) and the same real subject2 lab-mocap
  data those PASS-gated certs already use.

ISOLATION (COORDINATOR.md §1): bodytwin only; `.venv-msk`; external lab-mocap drive read in place, never
written; no git commit/push; both files listed above are new.
