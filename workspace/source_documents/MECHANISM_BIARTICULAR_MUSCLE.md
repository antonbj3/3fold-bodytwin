# MECHANISM BIARTICULAR MUSCLE — Lombard's paradox from measured moment-arm ratios, on this twin's own geometry (2026-07-22)

**Status: HYPOTHESIS awaiting independent QC.** Confidence tier: **mixed** — (a) this twin's own live
OpenSim moment-arm geometry [high confidence, machine-computed, 4/4 regression-checked against
already-validated numbers to <0.11%]; (b) the equal-moment-arm counter-model falsification
[proof-grade — an algebraic identity, numerically confirmed to machine epsilon, not a fitted result];
(c) live-verified primary jumping/cycling literature with real quantitative content [high confidence,
verbatim NCBI abstracts, PMID+DOI cross-checked]; (d) the gastrocnemius energy-transfer *percentage*
itself [literature-sourced only this session — no jump trial exists in this repo, disclosed gap].

## 0. What this is / is not

Biarticular (two-joint) muscles span two joints, so their moment at each joint depends on **both**
joint angles. Two classic puzzles: **Lombard's paradox** (hamstrings + rectus femoris co-contract in
rising from sitting / cycling — apparent antagonists — yet both simultaneously assist hip **and** knee
extension) and **gastrocnemius's proximal-to-distal energy transfer** (ankle push-off power partly
sourced from knee-extension power via the biarticular gastrocnemius, not the monoarticular soleus
alone). This document builds a small, geometry-first, machine-checked model of both, **on this twin's
own already-validated OpenSim leg** (`docs/MECHANISM_MOMENT_ARM_VALIDATION.md`'s exact model), not a
toy example — and forces the naive "equal moment arm" counter-model as the pre-registered falsifier.

This is **not** a re-run of `docs/MECHANISM_MOMENT_ARM_VALIDATION.md` — that document swept **one**
coordinate at a time (hip **or** knee, with the other frozen at 0°) to validate peak magnitudes against
cadaveric literature. This document queries **both coordinates together** at real combined sit-to-stand
postures (the actual quantity Lombard's paradox is about — a muscle's moment arm at joint A while joint
B is simultaneously flexed), adds a moment arm that prior document never measured (rectus femoris at
the hip), and asks a genuinely different question: not "is the peak magnitude right" but "does the
real 2×2 moment-arm matrix admit a co-contraction solution that extends both joints, and does the naive
equal-arm matrix not."

## 1. Pre-registration

- **C** (task's own framing): with the model's *real* moment arms, there exists a co-contraction force
  pair `F_hamstrings>0, F_rectus_femoris>0` producing **net extension at both the hip and the knee**
  simultaneously. **¬C**: no such pair exists (the system is "locked" — any co-contraction either does
  nothing net or flexes at least one joint).
- **Falsifier / forced adversary**: a counter-model where **each muscle is given one equal-magnitude
  moment arm applied at both joints** (the naive "single lever arm" simplification) must **fail** to
  produce a dual-extension solution — geometrically, this construction makes the muscle-to-joint-torque
  matrix **rank-1** (singular, σ_min=0), so C must hold for the real matrix (σ_min bounded away from 0)
  and ¬C must hold with **certainty, not just on average**, for the naive matrix, across the same pose
  set — a false positive here would be accepting C from a construction that also "works" for the
  falsifier; a premature negative would be declaring ¬C from one unlucky pose without orienting on why.
- **External anchor (decorrelated, never a tautology gate)**: Visser et al. 1990 (PMID 2079066,
  cadaveric, n=6 limbs — an anchor this twin's model lineage was **not** tuned to, see §6) for the
  moment-arm-ratio asymmetry itself; Andrews 1987 (PMID 3611133) and Gregor et al. 1985 (PMID 4008501)
  for EMG/kinetics-based confirmation (and a genuine nuance) of the paradox in real cycling; Prilutsky &
  Zatsiorsky 1994 (PMID 8106533), Bobbert & van Ingen Schenau 1988 (PMID 3379084), and Gregoire et al.
  1984 (PMID 6511147) for quantitative gastrocnemius/rectus-femoris energy-transfer numbers in jumping.
- **Threshold for gastrocnemius**: literature reports a specific proximal→distal transfer magnitude
  (Prilutsky & Zatsiorsky: 178.6±45.7 J, 97.1±27.2% of hip-joint work, squat jump). This session does
  **not** re-derive that percentage (no jump trial exists in this repo) — pre-registered as an honest,
  disclosed gap, not silently patched.

## 2. Method

Script: `scripts/msk/biarticular_lombard_paradox.py` (new, self-contained, re-runnable). Model: the
**same** `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` (subject2) already validated in
`docs/MECHANISM_MOMENT_ARM_VALIDATION.md`, read-only. Sign convention **reused, not re-derived**: that
document already established, live, on this exact model (via a finite-difference check on `vasint_r`),
that `computeMomentArm(state, coord)` is positive iff the muscle's tension would increase the
coordinate, and that both `hip_flexion_r` and `knee_angle_r` are flexion-positive.

1. **Regression/spot-check gate (run first, before trusting anything new)**: reproduce 4 numbers
   already in the validated CSV/JSON using this new script's own code path.
2. **A genuinely new measurement**: `recfem_r`'s moment arm about `hip_flexion_r` — absent from the
   prior sweep (its `HIP_FLEXION_MUSCLES` dict never included rectus femoris).
3. **The 2×2 model**: `R = [[ma_H_hip, ma_R_hip], [ma_H_knee, ma_R_knee]]` — hamstrings (`H`) lumped via
   an Fmax-weighted resultant across `semimem_r/semiten_r/bflh_r` (same construction the prior document
   used for adductor magnus's 4 bundles, §8 there); rectus femoris (`R`) is a single muscle. Extension
   generalized-force vector `Ext(F_H,F_R) = -R @ [F_H,F_R]` (negated because both coordinates are
   flexion-positive). C holds at a pose iff an open region of `F_H,F_R>0` makes both components of `Ext`
   positive — solved **two independent ways** (closed-form ratio-interval algebra; a 20,001-point
   angular numeric sweep, sign-assumption-free) and required to agree.
4. **The geometric governor**: `σ_min(R)`, the smallest singular value — R must be full rank for the two
   muscles to reach two independent torque directions at all.
5. **The falsifier**: an equal-moment-arm counter-model, `R_equal`, built by giving each muscle **one**
   |moment-arm| magnitude (mean of its own hip/knee magnitudes) at both joints, same sign/direction.
6. Repeated across **13 poses**: 7 "diagonal" sit-to-stand-like poses (hip and knee flexing together,
   0° to 100°) + 6 "off-diagonal" posture-adversary poses (hip and knee flexed independently) — the
   task's own named perturbation ("posture changes could flip the paradox") tested directly, not assumed
   away.
7. **Gastrocnemius**: live-measured `gasmed_r`/`gaslat_r` moment arms at **both** the knee and the ankle
   simultaneously near a push-off-representative pose (knee 20° flexed, ankle 20° plantarflexed) — a
   geometric necessary-condition check, explicitly **not** a re-derivation of the literature's dynamic
   power-transfer percentage (§7).

## 3. Regression/spot-check gate — PASS, 4/4, before any new number is trusted

| quantity | this script | already-validated value | rel. diff |
|---|---:|---:|---:|
| `recfem_r` knee moment arm @ 24.1° | −5.6368 cm | −5.64 cm (`MOMENT_ARM_VALIDATION.md` §5) | 0.056% |
| `semimem_r` hip moment arm @ 0° | −4.4756 cm | −4.48 cm (`MOMENT_ARM_VALIDATION.md` §4) | 0.099% |
| `gasmed_r` knee moment arm @ 86.1° | 3.3565 cm | 3.36 cm (`MOMENT_ARM_VALIDATION.md` §5) | 0.103% |
| `gaslat_r` knee moment arm @ 88.1° | 3.5334 cm | 3.53 cm (`MOMENT_ARM_VALIDATION.md` §5) | 0.096% |

All 4 within 0.11% — this script's state handling, units, and sign convention reproduce the prior,
already-validated computation before any new claim (recfem's hip moment arm, the combined-pose 2×2
matrix, the gastrocnemius ankle sweep) is built on it.

## 4. New measurement — rectus femoris hip-flexion moment arm (never measured in this repo before)

| hip flexion (deg) | 0 | 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `recfem_r` moment arm (cm) | 3.62 | 4.21 | 4.68 | 5.04 | 5.26 | **5.33** | 5.24 | 4.97 | 4.54 | 3.93 | 3.18 |

Positive throughout (confirming, by direct measurement rather than assumption, that `recfem_r` is a hip
**flexor** on this flexion-positive coordinate) and peaked at ~50° — a non-monotonic, physiologically
plausible shape (mirroring the non-monotonic hip-flexor moment-arm curves already established for
iliopsoas-type muscles in this model family).

## 5. Does the real 2×2 matrix produce Lombard's paradox? — pose-by-pose, machine-checked

| hip° | knee° | ρ_H=\|ma_H,hip/ma_H,knee\| | ρ_R=\|ma_R,hip/ma_R,knee\| | det(R) | σ_min (cm) | feasible? | ratio band F_H/F_R |
|--:|--:|--:|--:|--:|--:|:--:|---|
| 100 | 100 | 0.849 | 1.518 | −0.000598 | 0.916 | **NO** | — |
| 90 | 90 | 0.972 | 1.413 | −0.000547 | 0.705 | **NO** | — |
| 70 | 70 | 1.162 | 1.158 | +0.000009 | 0.009 | yes (knife-edge) | [0.935, 0.938] |
| 50 | 50 | 1.306 | 1.027 | +0.000674 | 0.633 | **yes** | [0.878, 1.117] |
| 30 | 30 | 1.402 | 0.894 | +0.001239 | 1.175 | **yes** | [0.827, 1.297] |
| 10 | 10 | 1.321 | 0.752 | +0.001282 | 1.340 | **yes** | [0.780, 1.370] |
| 0 | 0 | 1.180 | 0.691 | +0.001028 | 1.167 | **yes** | [0.767, 1.310] |
| 100 | 40 | 0.754 | 0.605 | +0.000358 | 0.426 | **yes** | [0.970, 1.209] |
| 40 | 100 | 1.659 | 2.257 | −0.000511 | 0.559 | **NO** | — |
| 90 | 20 | 0.996 | 0.706 | +0.000668 | 0.745 | **yes** | [0.976, 1.375] |
| 20 | 90 | 1.395 | 1.532 | −0.000169 | 0.188 | **NO** | — |
| 60 | 100 | 1.543 | 2.300 | −0.000651 | 0.728 | **NO** | — |
| 100 | 60 | 0.723 | 0.683 | +0.000089 | 0.108 | yes (knife-edge) | [0.978, 1.034] |

**8/13 poses feasible** (both solution methods — closed-form and the 20,001-point numeric sweep — agree
at **13/13 poses**). Feasibility tracks the geometric condition ρ_H > ρ_R exactly (det(R)>0 in every row
where and only where feasibility holds) — the Lombard's-paradox condition is confirmed to be precisely
"hamstrings' own hip:knee moment-arm ratio must exceed rectus femoris' own hip:knee ratio," not a vaguer
"moment arms differ" statement.

### 5.1 The unexpected result was ORIENTED, not shrugged off (OODA, not a one-shot negative)

The initial grid's 5 infeasible poses **(100,100), (90,90), (40,100), (20,90), (60,100)** all share one
thing: **knee ≥ 90°**. This was not accepted as a vague "sometimes it fails" — it was traced to a
mechanism via two decisive follow-up sweeps:

- **Fixed hip=90°, knee swept 10°→130° every 5°**: feasible for knee∈[10°,65°], **flips to infeasible**
  for knee∈[70°,115°] (det crosses zero between 65° [+0.000016] and 70° [−0.000107]), **flips back** to
  feasible for knee≥120° (det crosses zero between 115° [−0.000175] and 120° [+0.000081]).
- **Fixed knee=90°, hip swept 0°→120° every 10°**: infeasible at **every single hip angle tested**
  (det from −0.000098 to −0.000547, never crosses zero) — confirming the gating variable is **knee**
  angle, not hip angle, at least across this range.

**Mechanism, measured not asserted**: `recfem_r`'s own knee (extensor) moment-arm magnitude is
non-monotonic across this range — 5.63 cm at knee=20° → 2.30 cm at knee=100° → 3.09 cm at knee=120°
— a shrinking-then-recovering quadriceps mechanical advantage at deep flexion, the same qualitative
shape as the independently well-documented quadriceps/patellofemoral moment-arm-vs-flexion literature
(plausibly consistent with this repo's own `docs/MECHANISM_PATELLOFEMORAL_FORCE.md`, **not
cross-checked against that doc's own numbers this session** — disclosed gap). This shrinkage is what
flips ρ_R above ρ_H in the 70-115° band.

**Correct scope of this finding (symmetric QC, stated plainly)**: this is a fact about whether
**hamstrings+rectus-femoris ALONE** can produce dual extension — it is **not** a claim that real
sit-to-stand or squatting fails to extend both joints in that range (it empirically does — people rise
from deep chairs routinely). Real dual extension in that range is supplied by vasti (knee) and glutes
(hip), the actual prime movers; the biarticular pair's job there is not necessarily to supply net
extension unassisted. Read §5's table as: *the geometric condition for the "textbook" two-muscle-alone
resolution of Lombard's paradox holds across the majority of tested postures (extension to moderate
flexion, and very deep flexion) but has a real, mechanistically-explained, measured gap at
intermediate-to-deep knee flexion* — a sharper and more informative result than a flat yes/no.

**A suggestive, disclosed-as-unconfirmed convergence**: Gregor, Cavanagh & LaFortune 1985 (PMID
4008501) — using real cycling EMG + pedal-reaction kinetics, a wholly different method/era/data type —
independently found a genuine **knee flexor** net moment appearing "approximately half way through the
propulsive phase" of cycling, i.e., dual-extension-via-co-contraction does not hold at every crank angle
in real human cycling either. Cycling's knee-flexion ROM plausibly overlaps this session's measured
70-115° gap zone, but exact cycling joint-angle time series were **not** pulled this session to confirm
the overlap numerically — reported as a suggestive qualitative convergence, not a proven match.

## 6. The falsifier: equal-moment-arm counter-model — FAILS at every pose, exactly as predicted

Construction: each muscle gets **one** |moment-arm| magnitude (mean of its own hip/knee magnitudes),
applied at both joints with its own real sign preserved.

**Algebraic identity (proven before running anything, then confirmed numerically)**: forcing
`|ma_hip|=|ma_knee|` for both muscles makes `det(R_equal) = 0` **exactly**, whenever each muscle's own
action is the textbook-opposite pattern at its two joints (hamstring: hip−/knee+; RF: hip+/knee−) — this
is not a fitted or approximate result, it is a rank-1 (singular) matrix by construction.

| | measured across all 13 poses |
|---|---|
| `det(R_equal)` | 2.9×10⁻¹⁹ to ~1×10⁻¹⁷ — machine-epsilon zero, every pose |
| `σ_min(R_equal)` | 5.3×10⁻¹⁸ to ~1×10⁻¹⁷ — machine-epsilon zero, every pose |
| feasible dual-extension poses | **0 / 13** |

**Physical reading**: with equal moment arms, `M_knee ≡ −M_hip` for **any** `(F_H,F_R)` — the naive model
can only *transfer* moment between joints (extend one exactly as much as it flexes the other), never
*add* net extension at both simultaneously, for any co-contraction ratio, at any tested posture. This is
the precise, geometry-derived version of the task's own falsifier statement ("a model with equal moment
arms... predicts a locked/indeterminate system") — confirmed with certainty, not just on average, and
contrasted directly against the real model's pose-dependent-but-frequently-feasible pattern (§5).

## 7. Sharper test: does "equal activation" actually land inside the feasible band?

A feasible ratio *interval* existing is necessary but not sufficient — does the simplest physiologically
naive hypothesis (both muscles at the same fraction of their own Fmax) actually fall inside it?

`Fmax_hamstrings_total / Fmax_recfem` = **1.873** (native model Fmax) or **1.635** (Handsfield-PCSA-
corrected Fmax, `data/msk_models/subject2_scaled_handsfield_fmax_corrected.osim` — a sensitivity check,
since `docs/MECHANISM_SPECIFIC_TENSION.md`/`docs/MECHANISM_FMAX_PCSA_VALIDATION.md` already found the
native Fmax somewhat over-strong).

**Finding: 0/8 feasible poses have either ratio inside that pose's own feasible band** (e.g., pose
(0°,0°) band is [0.767, 1.310]; pose (70°,70°) band is [0.935, 0.938] — both ratios sit above every
tested band's upper bound). Equal relative activation is **not** the geometric solution — real
co-contraction achieving dual extension needs relatively **more** rectus-femoris force than an
equal-%Fmax assumption gives. This is consistent with (not proof of) real EMG evidence: `docs/
MECHANISM_EMG_CONSTRAINED_FORCE.md`/`docs/MECHANISM_EMG_TIMING.md` already measured real, asymmetric,
phase-specific hamstring-vs-quadriceps EMG timing for this same subject's gait trial (reused as a
pointer, not re-derived here); Andrews 1987 and Gregor 1985 (both above) likewise describe genuinely
asymmetric, phase-varying activation, never a matched-%MVC co-contraction.

## 8. Gastrocnemius — knee-to-ankle energy transfer: geometric prerequisite + literature anchor

**This twin's own geometry** (push-off-representative pose: knee 20° flexed, ankle 20° plantarflexed):

| | knee moment arm (cm) | ankle moment arm (cm) | ankle/knee gear ratio |
|---|--:|--:|--:|
| `gasmed_r` | 2.65 | −5.00 | 1.88 |
| `gaslat_r` | 2.83 | −5.15 | 1.82 |

Both gastrocnemius heads carry a **substantial, same-order-of-magnitude** moment arm at both joints
simultaneously (neither near zero) — the **necessary** (not sufficient) structural/geometric
prerequisite for proximal-to-distal power transfer to operate at all. This is a **new**, this-twin-
specific geometric measurement; it does **not** reproduce the literature's own dynamic power-transfer
*percentage* (no jump kinematics/kinetics/GRF trial exists in this repo's corpus — an honest, disclosed
gap, not silently patched).

**The literature anchor, live-verified this session, real quantitative content**:

- **Prilutsky & Zatsiorsky 1994** (PMID 8106533, n=5 subjects, inverse-dynamics 4-link 8-muscle leg
  model): squat vertical jump, maximum proximal→distal energy transfer via rectus femoris +
  gastrocnemius = **178.6 ± 45.7 J (97.1 ± 27.2% of the work done by the hip joint moment)**; landing,
  maximum distal→proximal transfer = 18.6 ± 4.2 J (38.5 ± 36.4% of ankle-joint-moment work). **This is
  the single quantitative anchor for the "monoarticular-only model under-predicts ankle
  power/jump-height" claim** — it is not independently re-derived here.
- **Gregoire et al. 1984** (PMID 6511147): "a high power output of **3000-4000 W** was delivered in the
  ankle joints during plantar flexion... attributed to a sequential energy flow from hip to knee and
  ankle joints" via glutes+RF+vasti+gastrocnemius; biarticular muscles reach a high-force/
  low-contraction-velocity regime during transfer "without losses due to eccentric contractions."
- **Bobbert & van Ingen Schenau 1988** (PMID 3379084): EMG plateau-activation sequence *"m.
  semitendinosus, long head of m. biceps femoris, m. gluteus maximus, m. vastus medialis, m. rectus
  femoris, m. soleus, m. gastrocnemius"* (proximal-to-distal); mechanical-energy optimality "can be
  satisfied by transportation of energy via the biarticular m. rectus femoris and m. gastrocnemius."

## 9. Decorrelated external anchors — verified live this session (NCBI eutils + Crossref)

| # | citation | PMID / DOI | role |
|---|---|---|---|
| 1 | Lombard WP (1903). *The Action of Two-Joint Muscles.* Am Phys Educ Rev 8(3):141-145 | DOI [10.1080/23267224.1903.10649915](https://doi.org/10.1080/23267224.1903.10649915) (Crossref — not in PubMed, pre-MEDLINE era, confirmed via a 0-result esearch first) | the original historical description |
| 2 | Andrews JG (1987). *Lombard's Paradox revisited.* J Biomech 20(6):565-75 | PMID [3611133](https://pubmed.ncbi.nlm.nih.gov/3611133/) | EMG/kinetics confirmation, cycling |
| 3 | Gregor RJ, Cavanagh PR, LaFortune M (1985). *Knee flexor moments... a creative solution to Lombard's Paradox.* J Biomech 18(5):307-16 | PMID [4008501](https://pubmed.ncbi.nlm.nih.gov/4008501/) | real net-knee-flexor-moment nuance, cycling |
| 4 | Visser JJ, Hoogkamer JE, Bobbert MF, Huijing PA (1990). *Length and moment arm of human leg muscles...* Eur J Appl Physiol 61(5-6):453-60 | PMID [2079066](https://pubmed.ncbi.nlm.nih.gov/2079066/) | cadaveric (n=6 limbs), independent moment-arm-ratio confirmation |
| 5 | Bobbert MF, van Ingen Schenau GJ (1988). *Coordination in vertical jumping.* J Biomech 21(3):249-62 | PMID [3379084](https://pubmed.ncbi.nlm.nih.gov/3379084/) | EMG sequence + energy-transportation argument, jumping |
| 6 | Gregoire L, Veeger HE, Huijing PA, van Ingen Schenau GJ (1984). *Role of mono- and biarticular muscles in explosive movements.* Int J Sports Med 5(6):301-5 | PMID [6511147](https://pubmed.ncbi.nlm.nih.gov/6511147/) | quantitative ankle power (3000-4000 W), jumping |
| 7 | Prilutsky BI, Zatsiorsky VM (1994). *Tendon action of two-joint muscles...* J Biomech 27(1):25-34 | PMID [8106533](https://pubmed.ncbi.nlm.nih.gov/8106533/) | quantitative energy-transfer J and %, jumping/landing |

All 7 confirmed via NCBI eutils `esearch`→`esummary` (title/journal/author match)→`efetch` (verbatim
abstract text), or (#1) Crossref, live this session — not recalled from memory. Reused without
re-verification (already machine-verified in this repo, see `docs/MECHANISM_MOMENT_ARM_VALIDATION.md`
§3): Nemeth & Ohlsen 1985 (PMID 3988782), Spoor & van Leeuwen 1992 (PMID 1733995), Buford et al. 1997
(PMID 9422462).

**Why Visser 1990 is a genuinely decorrelated anchor, not a tuned-to one**: `MOMENT_ARM_VALIDATION.md`
§3/§9 already established, live, which papers Rajagopal et al. 2016 (this twin's model lineage) cited
*by name* when validating its own geometry — Buford 1997, Spoor & van Leeuwen 1992, Nemeth & Ohlsen
1985, Dostal 1986, Arnold 2000. Visser 1990 is **not** on that list — an independent cadaveric dataset
this twin's geometry was never fit to. Its qualitative finding — *"the effect of changing the hip angle
on the biceps femoris muscle length is much larger than that of changing the knee angle. For the rectus
femoris muscle the reverse was found"* (length-angle sensitivity ∝ moment arm) — is **exactly** the
ρ_hamstring>1, ρ_rectus-femoris<1 pattern this session's live OpenSim measurement also finds at most
poses (§5), from a wholly independent data source and era.

## 10. Symmetric QC — what would make this wrong, checked

- **Could the "feasible" result be a methodology artifact?** No — 4/4 regression spot-checks vs
  already-validated numbers pass to <0.11% (§3), and the feasibility test itself is solved two
  independent ways (closed-form + numeric sweep) that agree at 13/13 poses.
- **Could the equal-arm falsifier's "always fails" be a coincidence of the specific construction
  chosen?** No — `det(R_equal)=0` is an algebraic identity for *any* equal-magnitude construction given
  the muscles' own (measured, not assumed) sign pattern, proven before running code, then confirmed to
  machine epsilon at all 13 poses.
- **Could the "8/13 feasible" result be cherry-picked toward a positive story?** No — the 5 infeasible
  poses (including a "knife-edge" pair at det≈0.00001-0.00009, i.e. two of the 8 "feasible" poses are
  themselves barely inside the boundary) are reported in full (§5 table), the gating mechanism was
  actively hunted down (§5.1) rather than left as an unexplained mixed result, and the honest scope
  limit (the 2-muscle-alone system, not real in-vivo dual extension) is stated plainly rather than
  implied away.
- **Could the "equal activation ratio never lands in-band" finding reverse under the corrected Fmax?**
  Checked directly (§7): native-model ratio 1.873, Handsfield-PCSA-corrected ratio 1.635 — both give the
  same qualitative answer (always above the feasible bands found), the conclusion is robust to this
  specific sensitivity check.

## 11. Honest gaps (full list)

1. The 2-muscle (hamstrings+RF) system is the *strongest/purest* form of the paradox test — real
   sit-to-stand/cycling also recruits vasti and glutes simultaneously, the actual prime movers; the
   §5.1 "gap zone" is a fact about the biarticular pair alone, not a claim that real dual-joint
   extension fails there (it does not).
2. Static/quasi-static geometry only — no forward dynamics, no force-length/velocity effects, no real
   EMG-driven force magnitudes for hamstrings/`recfem_r` during an actual sit-to-stand trial in this
   repo (`recfem_r` has **no** right-side EMG channel in this subject's dataset at all — an
   already-disclosed gap in `docs/MECHANISM_EMG_CONSTRAINED_FORCE.md` §2, reused not rediscovered). A
   real STS EMG trial (`STS1_EMG.sto`) exists in this repo and could exploit the hamstring channels in
   a follow-up; not attempted this session.
3. The gastrocnemius energy-transfer *percentage* (Prilutsky's 97.1±27.2%, Gregoire's 3000-4000 W) is
   literature-sourced, not independently re-derived — no jump trial exists in this repo's corpus.
4. Single subject (subject2), right leg only, one model lineage (LaiArnold/Rajagopal) — the scope
   caveat this entire document family discloses every time.
5. The equal-arm construction (mean of \|hip\|,\|knee\| magnitudes) is one reasonable operationalization
   of "naive equal moment arm," not the only possible one — though the qualitative result (det=0
   exactly) holds for any equal-magnitude construction by the same algebraic identity.
6. The §5.1 mechanism (non-monotonic `recfem_r` knee moment arm) was linked qualitatively to the
   well-known quadriceps/patellar moment-arm literature but **not** cross-checked against this repo's
   own `docs/MECHANISM_PATELLOFEMORAL_FORCE.md` numbers this session.
7. Gregor 1985's cycling knee-flexor-moment finding is reported as a suggestive convergence with this
   session's gap zone, **not** a numerically confirmed match (cycling joint-angle time series not
   pulled this session).
8. Only 2 axis-aligned slices (hip=90° fixed, knee=90° fixed) were finely swept to locate the crossover
   angles; the crossover's hip-angle-invariance was checked only across hip∈[0°,120°] at knee=90°, not
   exhaustively across the full pose space.

## 12. Verdict

- **C (real moment arms can produce Lombard's paradox) is NOT falsified, and is confirmed across the
  majority of tested postures (8/13, including the entire extension-to-moderate-flexion range that
  covers most of gait, cycling, and the terminal phase of sit-to-stand)** — with a real, honestly
  disclosed, mechanistically-explained gap at intermediate-to-deep knee flexion (~70-115°) for the
  two-muscle-alone system specifically.
- **¬C for the naive equal-moment-arm counter-model is confirmed with certainty — 0/13 poses feasible,
  σ_min ≈ machine-zero at every single pose** — a clean, algebraically-forced, universally-observed
  falsification exactly matching the task's own predicted "locked/indeterminate" failure mode.
- **A sharper, more honest secondary finding**: even where the paradox resolves geometrically, the
  naive "equal relative activation" hypothesis never lands inside the feasible band (0/8) — real
  co-contraction must be asymmetric, consistent with (not proven by) this repo's own already-measured
  asymmetric EMG timing and with the cycling EMG literature.
- **Gastrocnemius**: the geometric necessary-condition (substantial, simultaneous knee+ankle moment
  arms) is confirmed on this twin's own geometry; the literature's quantitative transfer numbers
  (Prilutsky & Zatsiorsky 1994, Gregoire et al. 1984, Bobbert & van Ingen Schenau 1988) are live-verified
  but not independently re-derived — honestly scoped as a literature anchor, not a twin-native result.

## Files

- `scripts/msk/biarticular_lombard_paradox.py` — the full, self-contained, re-runnable pipeline (loads
  the same read-only external model `docs/MECHANISM_MOMENT_ARM_VALIDATION.md` validated; regression
  spot-checks; new `recfem_r` hip-moment-arm sweep; 13-pose 2×2 matrix analysis; equal-arm falsifier;
  gastrocnemius ankle-transfer geometry).
- `data/msk_smoketest/biarticular_lombard/biarticular_lombard_results.json` — full raw output.
- `docs/MECHANISM_BIARTICULAR_MUSCLE_evidence.json` — curated evidence summary (citations, all key
  numbers, honest gaps, couples_to).
- Read-only inputs: `/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/
  OpenSimData/Mocap/Model/LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`;
  `data/msk_models/subject2_scaled_handsfield_fmax_corrected.osim` (Fmax sensitivity check only);
  `data/msk_smoketest/moment_arm_validation/moment_arm_raw_curves.csv` (regression-check reference).

Isolation respected throughout: bodytwin only; no git commit/add/push; no existing file modified; all
new files created this session only.
