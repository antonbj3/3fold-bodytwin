# MECHANISM MUSCLE DAMAGE — eccentric-injury accumulation, a first-step extension of fatigue (2026-07-21)

**Scope, stated up front:** a **threshold/accumulation PROXY** for eccentric (lengthening-contraction)
muscle damage — explicitly **NOT** a sarcomere-popping mechanistic model (no titin, no half-sarcomere
population, no simulated myofibril failure; contrast `docs/MECHANISM_CROSSBRIDGE_MODEL.md`, which DOES
model sub-cellular cross-bridge kinetics, but for one muscle, decoupled from the twin). **Single bout**:
no repeated-bout protective adaptation is modeled. **Generic**: reuses this twin's already-published,
not subject-fitted, model/SO parameters. **Shares the SO common-mode**: walking1's force comes from the
same Static-Optimization solve `docs/MECHANISM_MUSCLE_FATIGUE.md`/`MECHANISM_METABOLIC_COST.md` already
flagged as an effort-minimizing under-estimate of true co-contraction. **First step.** Code:
`scripts/msk/muscle_damage.py`. Evidence: `data/msk_smoketest/muscle_damage/muscle_damage_results.json`
(regenerated every run — a HYPOTHESIS artifact, not a certified cell). ISOLATION (bodytwin `COORDINATOR.md`
§1): bodytwin only; all real mocap/SO files read in place, never mutated; no git commit/push; new files.

---

## 1. Why this, why now

`docs/MECHANISM_MUSCLE_FATIGUE.md` added the motor-unit fatigue history-state Hill/Millard structurally
cannot have. Eccentric **damage** is a *different* history state — non-recoverable within one bout,
triggered not by cumulative effort but specifically by **active force production while the fiber
LENGTHENS** (Proske & Morgan 2001, live-verified below — sarcomeres on the descending limb of the
length-tension curve stretch past filament overlap under load). This script also **directly connects
to `docs/MECHANISM_CROSSBRIDGE_MODEL.md`'s own finding**: that sub-cellular cross-bridge model
over-predicts eccentric force **2.57x** relative to the twin's real `Millard2012EquilibriumMuscle`
`ForceVelocityCurve` (its concentric-branch agreement is good, RMS diff 0.05; eccentric RMS diff 1.88).
Section 7 below reuses that model's own saved evidence to ask: if the cross-bridge model's higher
eccentric tension is closer to sub-cellular reality, how much would this twin's own damage estimate
change?

## 2. Sign convention — MACHINE-VERIFIED, and a disclosed divergence from the task's own prose

`Muscle.getFiberVelocity()` is **positive when the fiber LENGTHENS** (eccentric) and **negative when it
SHORTENS** (concentric) — re-derived here as this script's own G0 check (Pearson correlation against an
*independent* finite-difference of `fiber_length(t)`, **all 80 muscles**, not just the 2 this repo's
`muscle_spindle.py` originally checked): **median r = 1.0000, min r = 0.9897** (gaslat_r, still far
above the 0.9 gate) — **PASS**.

**Honest, disclosed divergence**: the task brief's own prose describes eccentric lengthening as
"negative fiber velocity" — the **opposite** of this machine-verified OpenSim convention (also
independently consistent with `crossbridge_contraction_model.py`'s own documented convention: "norm_v...
+1 = max lengthening"). This script uses the **verified** convention throughout (positive = lengthening
= eccentric), flagged here explicitly rather than silently switched without disclosure — exactly the
class of bug this project's own memory flags ("a sweep-all-fails may be a shared convention bug — test
with a synthetic/independent-method control," which is precisely what this G0 check is).

## 3. Citations — every PMID/DOI live-verified this session (NCBI eutils; WebSearch was session-quota-exhausted, confirmed live, the same pre-existing repo fallback)

| # | Reference | PMID / DOI | Role |
|---|---|---|---|
| 1 | Proske U, Morgan DL (2001). *Muscle damage from eccentric exercise: mechanism, mechanical signs, adaptation and clinical applications.* J Physiol 537(Pt 2):333-45. | 11731568, `10.1111/j.1469-7793.2001.00333.x`, PMC2278966 (OA, fetched) | **Primary mechanism review** — sarcomere-popping on the descending limb; direct quotes fetched: *"only one form of exercise, eccentric exercise... leaves us stiff and sore the next day"*; *"one commonly encountered example of eccentric exercise is downhill walking"*; soreness *"sets in several hours later and peaks at about 48h."* Gives **no quantitative** eccentric/concentric ratio (checked, not assumed absent). |
| 2 | Brooks SV, Zerba E, Faulkner JA (1995). *Injury to muscle fibres after single stretches of passive and maximally stimulated muscles in mice.* J Physiol 488(Pt 2):459-69. | 8568684, `10.1113/jphysiol.1995.sp020980`, PMC1156684 (OA, fetched) | **Quantitative magnitude anchor**: *"the work done to stretch the muscle was the best predictor of the magnitude of injury"* (R²=0.76 active / 0.85 passive); **30% strain = significant-force-deficit threshold** in maximally activated muscle (used below). Passive muscle needs *more* strain at moderate levels but shows *greater* deficits above ~50% strain than active — a genuine, disclosed non-monotonic reversal. |
| 3 | Lieber RL, Fridén J (1993). *Muscle damage is not a function of muscle force but active muscle strain.* J Appl Physiol 74(2):520-6. | 8458765, `10.1152/jappl.1993.74.2.520` (no PMC, paywalled; abstract fetched live) | **Forced adversary** for the force-weighted proxy: identical 25% (then 12.5%) strain at very different force (via activation timing) → equivalent damage; strain-magnitude effect P<0.001, force/timing effect P>0.7. Motivates the force-*blind* "active strain rate" variant computed alongside, not swept aside. |
| 4 | McHugh MP (2003). *Recent advances in the understanding of the repeated bout effect.* Scand J Med Sci Sports 13(2):88-97. | 12641640, `10.1034/j.1600-0838.2003.02477.x` | Justifies the **single-bout-only** scope limit — a second exposure would genuinely damage less; not modeled here. |
| 5 | Rajagopal A et al. (2016). *Full-Body Musculoskeletal Model for Muscle-Driven Simulation of Human Gait.* IEEE Trans Biomed Eng 63(10):2068-79. | 27392337, PMC5507211 | **Phase-timing anchor, REUSED** (already live-verified in this repo's `validate_emg_timing.py`, not re-searched): `vet.EMG_LITERATURE`'s per-muscle expected %GC active windows, direct quotes fetched there previously. |

Two of three newly-pulled PMIDs this session were **wrong before verification** (Lieber/Fridén's real
PMID 8458765 vs. a recalled 8458772; McHugh 2003's real PMID 12641640 required disambiguating 3
esearch candidates by title, not the top hit) — consistent with this project's own previously-measured
~62% recalled-citation drift rate.

## 4. Model — two decorrelated proxy variants, forced against each other

- **PRIMARY ("eccentric work")**: `damage_rate(t) = F_norm(t) · max(0, v_norm(t))`, gated to zero
  whenever activation ≤ 0.05 (see §6.1's forced fix) — normalized active force × normalized lengthening
  velocity, the mechanical-power/negative-work quantity ref [2] found was its **own** best empirical
  predictor.
- **ADVERSARY ("active strain rate"), forced by ref [3]**: force-*blind* — `activation(t)>0.05 ?
  max(0, fiber_velocity_ms(t)) : 0` — directly motivated by Lieber & Fridén's finding that force
  magnitude does not predict damage, active strain does.
- Both integrate to a **monotone, non-recovering** accumulator (no repair term — single bout, §8).
- **Robustness check**: Spearman rank-correlation between the two variants' per-muscle totals, n=80
  muscles: **r = 0.976** (gate >0.5: **PASS** — the two literature positions agree on which muscles rank
  highest here, though they need not have).

## 5. Walking1 results — per muscle/phase (primary, quantitative; all 80 muscles reconstructed)

Reuses the already-validated SO activation+force (`static_optimization/so/`, zero re-solve) plus a
rigid-tendon fiber-velocity reconstruction (the same forced-fix this repo's `muscle_spindle.py`/
`metabolic_cost.py` already established: compliant tendon makes `getFiberVelocity()` silently return
exactly 0 at every independently-solved frame).

**%GC phase-binning vs. the REUSED Rajagopal-2016 activation-timing anchor** (`vet.EMG_LITERATURE`,
`vet.MUSCLE_GROUPS`, zero re-derivation):

| group | peak ecc.-work %GC | expected window(s) | verdict | D_work (1 cycle) |
|---|---:|---|---|---:|
| Vasti (quad) | 7 | (0,30) early stance | **PASS** | 0.000921 |
| Rectus femoris | 58 | (0,30), (50,73) | **PASS** | 0.010520 |
| Hamstrings | 91 | (0,30), (85,100) terminal swing | **PASS** | 0.002545 |
| Gastrocnemius | 39 | (30,60) late stance | **PASS** | 0.004755 |
| Soleus | 22 | (30,60) late stance | **MISS** (see §6.2) | 0.002093 |
| Tibialis anterior | 4 | (0,30), (60,100) | **PASS** | 0.001750 |
| Gluteus maximus | 10 | (85,100), (0,50) | **PASS** | 0.000019 |

**6/7 PASS** (gate ≥5/7). Top-5 muscles by D_work: recfem_r, recfem_l, tfl_r, tfl_l, tibant_l — **the
quadriceps' biarticular member (rectus femoris) and TFL dominate**, not the monoarticular vasti, exactly
the biarticular-compounding pattern §6.1 diagnoses in detail.

**Reading the confirmed pattern**: quads (vasti) peak eccentric loading in **early stance** (loading
response — the classic "controlled knee flexion" weight-acceptance phase); hamstrings peak in
**terminal swing** (decelerating the shank before heel-strike — the textbook hamstring-strain-injury
phase in sports medicine); rectus femoris peaks in **pre-swing/initial-swing** (decelerating rapid knee
flexion at toe-off); gastrocnemius peaks in **late stance**. All four match this repo's own
already-live-verified EMG-timing literature, independently, with zero new phase-window citation needed.

## 6. Two forced, symmetric surprises — both diagnosed via OODA, neither hidden nor hand-waved

### 6.1 Brooks-threshold gate FAILS — semimembranosus/rectus femoris exceed 30% active strain (real, geometric, not a bug)

Pre-registered gate: peak single-episode ACTIVE (mean activation >5%) strain, all 80 muscles, should
stay under ref [2]'s 30% mouse-EDL single-stretch injury threshold. **Measured: 55.33% (semimem_l),
51.86% (semimem_r), 50.47% (recfem_r), 46.88% (recfem_l) — gate FAILS.**

**A one-shot "surprise" was not accepted at face value (per the mandate: honest-negative is not a free
pass) — forced through OODA before being reported:**

1. **Observe**: is this a sign/numerical bug? The same muscle's own G0 sign-QC is **r=0.9999995**
   (near-perfect) — the finite-difference and the analytic `getFiberVelocity()` agree almost exactly, so
   this is not a divergent/glitching reconstruction.
2. **Orient (hypothesis 1 — rigid-tendon inflation)**: forcing rigid tendon (necessary for non-zero
   velocity) dumps the whole MTU excursion onto the fiber for long-tendon muscles — already a disclosed
   bias elsewhere in this repo. **Tested directly**: a *separate compliant-tendon* (model-default) pass,
   length-only, for the same episode gives semimem_l ptp = **50.2%** of Lopt — barely smaller than the
   rigid-tendon 55.3%. **Hypothesis 1 is insufficient on its own** — refuted as the *primary* cause, not
   silently kept.
3. **Orient (hypothesis 2 — genuine biarticular compounding, geometrically derived)**: instrumenting the
   raw IK trace over the episode (t=0.13–0.53s) shows `knee_angle_l` falling **63°→4°** (rapid extension)
   while `hip_flexion_l` rises **~9°→26°** (hip flexing) *simultaneously*. Semimembranosus is hip
   extensor + knee flexor: hip flexion lengthens it at the hip, knee extension lengthens it at the knee
   — **both joint actions compound rather than cancel**. **Machine-checked geometric cross-check** (live
   `computeMomentArm` at the model's own configuration, not a heuristic): knee moment arm ≈40–47mm ×
   59° ROM ≈ **44mm** + hip moment arm ≈50–56mm × 17° ROM ≈ **16mm** = **≈60mm predicted**, vs. **53–58mm
   measured** (compliant/rigid) — independent geometric confirmation, not asserted.
4. **Bilateral reproducibility**: semimem_r shows the *same* ~52% peak, at its own later terminal-swing
   instant (t≈0.85s, this leg's own phase within the clip) — the finding is not a one-leg sampling fluke.
5. **Verdict**: the 30–55% fiber-strain FAIL is a **real, geometrically-explained** finding — the
   classic biarticular hip-flexion + knee-extension compounding stretch of the hamstring during terminal
   swing (and analogously, rectus femoris' hip-extension + knee-flexion compounding during early
   stance/pre-swing) — **not** a code bug, but also **not** proof that walking is injurious in reality.
   Two genuine, disclosed uncertainties remain open (§8): (a) ref [2]'s 30% threshold is a mouse
   EDL/soleus *ex-vivo* single-stretch number — whether it transfers quantitatively to human
   hamstring/rectus-femoris *in-vivo* gait strain is untested here; (b) this model's own compliant-tendon
   parameters for these two muscles specifically provide little buffering (52.6%→50.2%, barely
   different from rigid), which may or may not match the true tendon/aponeurosis compliance. Read
   positively: this is precisely the mechanism sports medicine already identifies as the **most common
   site of real hamstring/rectus-femoris strain injury** (at running/sprinting demand, not gentle
   walking) — the twin's own measurement recovers the right mechanism and the right two muscles.

### 6.2 Phase-pattern gate — a genuine void-floor bug, found, diagnosed, fixed at the source (not tuned to pass)

First run: **4/7** groups matched (vasti and gluteus maximus MISSED). Forced diagnosis before accepting
the FAIL: instrumenting the exact per-frame values at vasti's "peak" %GC bin (63) showed activation
pinned at **SO's numerical floor (~0.010)** — not real engagement — while fiber velocity was large (fast
passive swing-phase motion). The primary metric's `F_norm × max(0,v_norm)` product is technically
non-zero there but is a **void-floor artifact** (floor-activation × high-velocity), not genuine
eccentric *loading* — exactly the failure mode this project's own methodology names explicitly
("run a void-floor sweep on big margins"). **Fix, forced not tuned**: gate the primary metric by the
*same* activation floor (>0.05) already used in the adversary variant — a symmetric correction applied
identically to both proxies, verified to actually resolve the diagnosed cases before being kept: vasti
now peaks at %GC=7 (PASS, early stance) and gluteus maximus at %GC=10 (PASS) — **6/7 PASS**.

**Soleus remains a MISS after the fix** (peak at %GC=22 vs. the expected (30,60) late-stance window) —
checked, and this is **not** the same floor artifact (soleus activation at %GC=22 is well above the
floor). This instead reflects a genuine, well-known **sub-phase** distinction the coarse activation-only
literature window cannot resolve: soleus/gastrocnemius are active through most of stance, but do
**eccentric** work (controlled ankle dorsiflexion as the tibia advances over the planted foot) in
**early-to-mid** stance, then switch to **concentric** push-off in **late** stance — ref [5]'s own window
is an *activation*-timing anchor, not an eccentric-specific one (disclosed in §8), so this MISS is a
genuine limitation of the anchor's resolution, not of the measurement.

## 7. Cross-bridge connection (Section required by task) — does the 2.57x eccentric over-prediction change the estimate?

Reuses `crossbridge_contraction_model_evidence.json`'s own saved `norm_v`/`crossbridge_fv`/`millard_fv`
arrays (eccentric branch, norm_v∈[0.05,0.95]) — **not recomputed**. Building the ratio
`crossbridge_fv/millard_fv` and applying it (interpolated at each frame's own measured v_norm) as a
sensitivity-scaled force multiplier on top of walking1's SO force, for the landing-relevant muscles:

| muscle | D_work (primary) | D_work (crossbridge-scaled) | scale |
|---|---:|---:|---:|
| vaslat_r | 0.000689 | 0.001055 | ×1.53 |
| vasmed_r | 0.000231 | 0.000341 | ×1.47 |
| recfem_r | 0.010520 | 0.021808 | ×2.07 |
| gasmed_r | 0.002147 | 0.003198 | ×1.49 |
| gaslat_r | 0.002608 | 0.003902 | ×1.50 |
| soleus_r | 0.002093 | 0.003157 | ×1.51 |

Because real walking's peak `v_norm` is small (0.17–0.24, i.e. 17–24% of Vmax) relative to the
cross-bridge model's tested eccentric-branch domain, the applicable ratio sits at the **low** end of its
own curve (1.45–2.1×, not the full 2.57× peak, which occurs near v_norm→0.95) — meaning the cross-bridge
finding matters **more** at higher velocities than walking ever produces. **Honest gap**: this sensitivity
could **not** be applied to DJ1 (no per-muscle force there, §8) — precisely where DJ1's much higher
measured v_norm (0.40–1.43, §9) would put it closer to the cross-bridge model's own strongly-diverging
high-velocity range, and where the correction would matter most.

## 8. DJ1 (real drop-landing) — kinematic + real-joint-power only, no Static Optimization solved

No SO exists for DJ1 in this repo; none was solved here (out of scope for a first step, against this
repo's "reuse, never re-solve" discipline). Instead: **rigid-tendon fiber velocity is a pure function of
joint geometry, independent of activation** — verified directly (not assumed): reconstructing DJ1's
fiber kinematics at two wildly different constant activations (0.02 vs 0.90) gives IDENTICAL fiber
length/velocity to machine precision (max diff 2.4e-11 m, 1.9e-11 m/s) — **PASS**. DJ1's own real,
precomputed inverse-dynamics moments (`OpenSimData/Mocap/ID/DJ1.sto`, produced by the original 2021
pipeline bundled with the dataset, reused read-only) give a **convention-free** net joint power
(moment × angular velocity — invariant to the coordinate's own flexion/extension sign convention).

Landing instant detected via the real force-plate GRF (reusing `wobbling_mass.py`'s own forced-fixed
`detect_contact_onset`, peak dGRF/dt): **t=1.5855s, RFD=141,975 N/s.**

**Peak lengthening v_norm, DJ1 landing window vs. walking1's own whole-trial peak** (pre-registered:
DJ1 > walking1 for every landing-relevant muscle):

| muscle | walking1 peak v_norm | DJ1 landing peak v_norm | ratio |
|---|---:|---:|---:|
| vasint_r | 0.206 | 0.520 | **2.53×** |
| vaslat_r | 0.199 | 0.504 | **2.53×** |
| vasmed_r | 0.196 | 0.497 | **2.54×** |
| recfem_r | 0.217 | 0.405 | **1.87×** |
| gasmed_r | 0.237 | 0.761 | **3.22×** |
| gaslat_r | 0.207 | 0.661 | **3.20×** |
| soleus_r | 0.167 | 1.431 | **8.59×** |

**All 7 ratios >1 — gate PASS.** Soleus shows the largest jump (the ankle is the first shock absorber at
impact). **Joint-power coincidence**: most-negative (energy-absorbing) net knee power in the landing
window is **−1407 W at t=1.620s**, 30ms from vaslat_r's own peak lengthening (t=1.59s); most-negative
ankle power is **−836 W at t=1.590s**, 10ms from soleus_r's peak lengthening (t=1.58s) — **PASS**
(high joint-level loading and peak fiber lengthening coincide in time, a convention-free physical check,
not assumed).

This matches Proske & Morgan (2001)'s own qualitative anchor — "downhill walking" as the classic everyday
eccentric-damage exercise, i.e. an *exaggerated* eccentric-lengthening task, not ordinary level gait — and
here a real drop-landing shows 1.9–8.6× the peak eccentric fiber velocity that a full stride of level
walking ever produces, for exactly the muscles (quads, plantarflexors) classically implicated in
landing/plyometric-training-induced DOMS.

## 9. Pre-registered gates — 7/8 PASS

```
G0_sign_convention_median_r_gt_0.9:            PASS  (median r=1.0000, min r=0.9897)
walking_active_strain_under_brooks_threshold:  FAIL  (55.3% vs 30% -- diagnosed real, §6.1, not a bug)
force_vs_strain_rank_robustness:               PASS  (Spearman r=0.976)
linear_accumulation_selfcheck:                 PASS  (ratio=1.0008, N=20 tiles)
phase_pattern_majority_match:                  PASS  (6/7, after the forced void-floor fix, §6.2)
dj1_activation_invariance_selfcheck:           PASS  (max diff ~2e-11)
dj1_eccentric_velocity_exceeds_walking:        PASS  (1.87-8.59x, all 7 muscles)
dj1_moment_velocity_coincidence:               PASS  (10-30ms)
```

**Overall: FAIL (1/8 gates), reported honestly, not forced to a hidden pass.** The one FAIL is a
genuine, geometrically-diagnosed, symmetric finding (§6.1), not a defect — kept as a pre-registered FAIL
per the mandate that a rejection/negative carries the same forced-steelman burden as a confirmation, and
an unjustified after-the-fact loosening of the threshold would be exactly the dishonest move the
methodology forbids.

## 10. Honest gaps (disclosed, not hidden)

- **Proxy, not mechanistic.** `D_work`/`D_strain` are relative loading indices with no physical units —
  useful for ranking muscles/phases, not an absolute clinical damage measure (contrast
  `crossbridge_contraction_model.py`'s own sub-cellular, one-muscle model).
- **Single bout only** — no repeated-bout protective adaptation (ref [4]) is modeled; a second exposure
  would genuinely damage less; this accumulator is monotone/non-recovering by construction.
- **Generic parameters** — this repo's already-published, not subject-fitted, SO/model parameters.
- **Shares the SO common-mode**: walking1's F_norm inherits SO's known effort-minimizing
  under-estimate of true co-contraction/force (same root cause `metabolic_cost.py` already flagged,
  externally verified against Koelewijn et al. 2019) — true eccentric loading is likely
  under-estimated here too, in the same direction, for the same reason.
- **DJ1 is kinematic + joint-power only, not per-muscle force** — no Static Optimization exists for this
  trial and none was solved here; DJ1's own D_work/D_strain are not computed, only fiber lengthening
  velocity (activation-invariance verified) and real net joint power.
- **Cross-species/cross-preparation threshold anchor**: ref [2]'s 30% strain threshold is a mouse
  EDL/soleus *ex-vivo* single-stretch number applied to human hamstring/rectus-femoris *in-vivo* gait
  strain — a genuine extrapolation, disclosed, same class this repo's `muscle_spindle.py` already made
  with cat afferent data.
- **The %GC phase anchor (ref [5]) is activation-timing, not eccentric-specific** — this script adds the
  eccentric/concentric split on top of that anchor via its own fiber-velocity-sign measurement; a
  window match is therefore necessary-not-sufficient (confirms the muscle is active in the right window;
  soleus's §6.2 MISS shows this distinction is real, not just a technicality).
- **DJ1's 100Hz IK** differentiates a very fast impact event (RFD≈142,000 N/s over a 10ms sample
  interval) — this repo's own `MECHANISM_WOBBLING_MASS.md` already found DJ1's 100Hz-IK second derivative
  (acceleration) has known validity problems; this script's velocity is only a first derivative
  (comparatively more robust) but the same underlying sampling-rate limitation applies, disclosed not
  hidden — DJ1's large v_norm values (up to 1.43) should be read as directionally strong, not
  millisecond-precise.
- **Cross-bridge sensitivity (§7) could not be applied to DJ1** — exactly where the 2.57x over-prediction
  would matter most (DJ1's velocities are far higher, closer to the cross-bridge model's own strongly
  diverging high-velocity range) is the one place this correction could not be computed.

## 11. Files

- `scripts/msk/muscle_damage.py` — the model + all 9 steps + self-checks + gates. Run with
  `.venv-msk/bin/python3 scripts/msk/muscle_damage.py` (~3s wall time).
- `data/msk_smoketest/muscle_damage/muscle_damage_results.json` — full evidence (citations, per-muscle
  episodes, damage proxies, phase report, cross-bridge sensitivity, DJ1 comparison, gates).
- This doc.
