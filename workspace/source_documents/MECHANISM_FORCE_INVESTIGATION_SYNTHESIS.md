# MECHANISM — Joint-Force Over-Prediction: Investigation Synthesis (the "connection")

**Status: the knee/hip contact-force over-prediction is understood, localized, decomposed, and partially
corrected — bounded honestly.** This doc connects ~20 decorrelated threads into one picture. Every claim is
pinned to an external anchor the model did not get to choose; confidence tiers are explicit. It supersedes the
scattered per-thread headlines and records where I corrected myself.

---

## 1. The headline, honestly bounded
The twin over-predicts in-vivo joint contact force during walking. Cohort range (subjects 2–4, walking1):
- **Knee 1.51× (subj2) → 1.99× (subj3)**, **hip 1.41× → 1.82×**, vs OrthoLoad in-vivo medians (knee 258 %BW,
  hip 274 %BW). Every ratio > 1.4×. `MECHANISM_CROSS_SUBJECT.md`.
- **It is real, not a solve artifact**: a decorrelated forward-dynamic solver (CMC) reproduces it within ~10%,
  slightly higher — the divergence guard passed (crossing-muscle sets byte-identical). `MECHANISM_CMC_SECOND_SOLVE.md`.
- **subject2 is the MILDEST case** — my earlier "~1.5×" headline was optimistic; the honest number is a range.
- ⚠️ **SAME-SUBJECT confirmation RETRACTED (units bug, 2026-07-22)** — the DM flagship compared the twin's
  muscle-free Tier-1 (104 %BW) against a **units-broken** anchor: the raw eTibia CSV is in **lbf**, read as N,
  so DM's "58–68 %BW" is actually **≈296 %BW** (437.9 lbf → 1948 N; a 0.67×BW peak knee force is physically
  impossible, being *below* the ground reaction force). Corrected, the muscle-free Tier-1 floor *under*-predicts
  DM — as a muscle-free bound should — and whether the *full* twin over-predicts DM is **unknown** (Tier-2
  blocked on DM's TKA knee). So the cross-cohort escape hatch is **NOT** closed at the knee by same-subject
  in-vivo force. The cohort measured peaks (JW 271, DM ≈294, SC 262 %BW; `MECHANISM_GRAND_CHALLENGE_MULTISUBJECT.md`)
  are consistent with the instrumented-knee literature. My "independent re-parse" reproduced the bug because it
  shared the same hidden units assumption — a caution that re-verification sharing an assumption is not
  decorrelated. **The cross-cohort over-prediction below (twin Tier-2 vs OrthoLoad, consistent units) is
  unaffected and still stands on its own** — but it is no longer independently confirmed by same-subject force.
  `MECHANISM_GRAND_CHALLENGE_DM.md` (corrected).
- ⚠️ **RECONCILIATION that reframes this whole synthesis (2026-07-22):** the twin's *full deformable-contact*
  knee model (OpenSim-JAM **COMAK**, the `KNEE-CELL` flagship in `MECHANISM_STATE.md`) predicted **2.503 BW vs DM's
  measured 2.584 BW = ~3% peak (r=0.802)** — i.e. on same-subject in-vivo ground truth the twin's best knee model
  **MATCHES**, it does not over-predict.
  ⚠️ **2026-07-24 — that number needs its caveats here too** (see `MECHANISM_STATE.md` §1 for the full three): it is
  **best-of-3** (honest range 3–10% peak / r 0.66–0.82, single subject/limb), it rests on **n_eff=2**, and a void-floor
  test shows **r and peak-error do not discriminate** this pipeline from a generic cross-cohort waveform given only real
  stance timing. The "MATCHES rather than over-predicts" reframing above still stands — it is about *quantity identity*
  (deformable contact vs a static-optimization free-body cut), and that conclusion is unaffected — but the 3% must not be
  read as evidence of subject-specific mechanistic accuracy. The 1.5× over-prediction documented below is the **LaiArnold
  static-optimization free-body-cut** reaction (a *different quantity* than deformable contact) on *subject2* vs a
  *cross-cohort* OrthoLoad median. So the over-prediction is most likely a **property of the SO free-body-cut
  method**, not a twin-physics error at the knee. A clean COMAK-vs-SO comparison on the *same* subject is the open
  work that would settle it; until then, read every "twin over-predicts" claim below as **SO-method-scoped**.

## 2. What generalizes vs what is subject/joint-specific (a correction I made)
- **Generalizes 8/8 (subjects × joints): the *directional* asymmetry** — the over-prediction concentrates at
  **push-off**, always more than at weight-acceptance. `MECHANISM_WAVEFORM_CROSS_SUBJECT.md`.
- **Does NOT generalize: the absolute "matches in-vivo at weight-acceptance."** That is specifically **knee + SO**
  and holds for 2 of 3 subjects (subj2 Z=0.97, subj3 Z=1.2), **not** subj4 (Z=2.77). The **hip is elevated across
  the whole cycle** — even subject2's hip weight-acceptance was Z=1.72. I overstated "the twin is correct at
  early stance"; corrected here.

## 3. WHY — the decomposition, joint by joint (this is the "connection")
Both major joints' over-predictions are **localized to a single over-loaded muscle, majority-correctable by
rerouting to synergists — NOT a diffuse everywhere-error.** The corrective *lever* differs:

| joint | localized cause | correction (independent anchor) | closes | reroutes to | residual |
|---|---|---|---:|---|---:|
| **Knee (push-off)** | over-strong **biarticular gastrocnemius** (Fmax too high vs MRI-PCSA) | MRI-PCSA Fmax correction — a *prediction*, grep-verified never fit to OrthoLoad | ~17% (push-off Z 2.67→1.23); **selective** — weight-acceptance match preserved | **soleus** (mono-articular, doesn't load the knee); total ankle force ~unchanged (+1.4%) | ~21% |
| **Hip** | **gluteus-medius over-recruitment** by SO vs real EMG (585N vs 100N) | EMG-informed redistribution (moment-budget-preserving LP/QP) | 55–82% (1.41→1.18→1.08) | other hip **abductors** (glmin, tfl) | 9–19% |

Threads: `MECHANISM_CONTACT_WAVEFORM.md`, `MECHANISM_CONTACT_MUSCLE_DECOMP.md`, `MECHANISM_FMAX_PCSA_VALIDATION.md`,
`MECHANISM_FMAX_CORRECTION_WAVEFORM.md`, `MECHANISM_PUSHOFF_PLANTARFLEXOR.md`, `MECHANISM_HIP_STRUCTURAL_MECHANISM.md`.

**The knee push-off chain is the tightest connection** — five independent threads all pointed at the biarticular
gastroc: (1) over-prediction is push-off-only → (2) push-off is gastroc-dominated (48.7% of knee contact) →
(3) gastroc is over-strong vs MRI-PCSA → (4) correcting its strength selectively fixes push-off → (5) the fix
works by rerouting to soleus, and the total Achilles force was only modestly high (Z=1.06 vs a direct in-vivo
transducer). Cause, mechanism, and validated fix, cross-checked five ways.

**Spine** corroborates the theme: the pure-kinematic (no-muscle) lumbar load ≈ VBR (56 vs 53 %BW); the muscle
contribution drives the raw 3.94× — but the VBR device's load-share fraction is an unmeasured confound.
`MECHANISM_SPINE_GAIT_VBR.md`.

## 4. Ruled OUT as the source (each an independent thread)
- Wobbling/soft-tissue mass (walking +0.02%) — `MECHANISM_WOBBLING_JOINT_FORCE.md`
- *Missing* co-contraction — imposing measured EMG makes the knee **worse** — `MECHANISM_EMG_CONSTRAINED_FORCE.md`
- Prime-mover moment-arm geometry — within cadaveric spread (but semi-circular: the model was fit to it) — `MECHANISM_MOMENT_ARM_VALIDATION.md`
- *Uniform* strength scaling and the SO objective exponent (p=2 already optimal) — `MECHANISM_OVERPREDICTION_DECOMP.md`
- The contact model **alone**: a heat-rate metabolic model that never touches contact geometry over-predicts
  calorimetry by the same ~1.5–1.6× → the excess is **upstream in the muscle forces**. `MECHANISM_METABOLIC_CALORIMETRY.md`.
  The deformable-contact (JAM) comparison is now a **DIAGNOSED-GAP**: the JAM invocation was proven to freeze
  every muscle at activation=0.05 regardless of input (two different inputs → identical frozen output), so the
  earlier "relieves some" / 0.7–0.9× number is **retracted** until the invocation is fixed. `MECHANISM_JAM_CONTACT_DECORR.md` + `MECHANISM_RESIDUAL_DECOMPOSITION.md`.

**Two corroborations that sharpen the localization (added):**
- **Ankle = control.** The ankle (talocrural) contact force *matches* its cadaveric anchor at push-off (1.02–1.09×,
  0/101 %GC points over the bar) → the push-off excess is not the ankle-plantarflexion demand but specifically the
  biarticular gastroc's **knee-crossing** component. `MECHANISM_ANKLE_WAVEFORM.md`.
- **Patellofemoral does NOT over-predict** (77 %BW, inside the 50–150 %BW band) → the excess is specific to the
  tibiofemoral + hip, not "everything." `MECHANISM_PATELLOFEMORAL_FORCE.md`.
- **Gastrocnemius is the recurring culprit on 3 axes**: over-strong (Fmax), fiber runs long vs ultrasound
  (85–95 vs ~57 mm), tendon-decoupling too weak — the same muscle, three linked errors. `MECHANISM_GASTROC_FASCICLE.md`.

## 5. The irreducible uncertainty (do not "fix" by tuning to the anchor)
Tendon-slack-length is poorly measured; within its ±10% band the knee contact ratio spans **0.78×–1.59×**
(`MECHANISM_TENDON_SLACK_SENSITIVITY.md`). The single-point headline hides this. Re-tuning l_TS to hit OrthoLoad
would be fitting-to-the-anchor, not validation.

**The residual decomposition makes this sharper (`MECHANISM_RESIDUAL_DECOMPOSITION.md`):** after the Fmax
correction, the knee push-off residual (1.214×) **dissolves into l_TS uncertainty** — a joint re-sweep on the
corrected model swings it 0.62×–1.31× (a band that *includes* in-vivo). So there is **no clean "irreducible
structural floor" at the knee** — the residual's very size is not well-determined. The **hip's 9–19% residual is
the one firm number** (l_TS cannot explain it; survives the Fmax + glute-med corrections) — but it is still
anchored to a cross-cohort OrthoLoad median (no same-subject in-vivo hip force exists) and most likely reflects
the SO-then-subtract method's known structural bias (Moissenet), cleanly isolated rather than a new phenomenon.
**And the Fmax fix itself is subject2-strong** (`MECHANISM_FMAX_FIX_CROSS_SUBJECT.md`): it closes 54% of
subject2's gap but only ~22% for subjects 3/4 — subject2 is the best case on every axis (mildest, most
correctable, cleanest solve), so the flagship numbers are its most-favorable. Bone stress does **not** coherently inherit the contact
over-prediction (de-inflating contact force *raises* it, via the Moissenet muscle–reaction coupling) —
`MECHANISM_BONE_STRESS_CONSEQUENCE.md`.

## 6. Why this is trustworthy without domain expertise (the structural handles)
Decorrelated external anchors, none of which the model chose the answer to: **OrthoLoad** in-vivo implants
(knee/hip/spine-VBR), indirect-**calorimetry** COT, **MRI-PCSA** muscle volumes (Handsfield), direct **Achilles
transducer** (Finni), the subject's own measured **EMG**, cadaveric **moment-arm** & **bone-strain** literature.
The evidence the structure works is that it **caught my errors** (below), not that I was right.

## 7. Corrections the structure forced (self-catching)
- "~1.5×" was the mildest subject; honest range 1.4–2.0×.
- "Twin matches at weight-acceptance" → knee+SO+2-of-3-subjects only; hip elevated throughout.
- MCL "58% strain" → **ε_r-fragile, 50–74%** — not robust as a number (`MECHANISM_LIGAMENT_EPSR_SENSITIVITY.md`).
- "0/1397 tension-only violations" → **tautological** (a force-law identity), not a validation.
- JAM "doesn't relieve" verdict → rested on a **donor-BW vs subject-BW normalization** mismatch; consistent
  units flip it toward borderline.
- Four wrong PMIDs/anchor-bands in my own subagent briefs, each caught by the agent (Miller2014, Handsfield PMID,
  the FE-vs-transducer Achilles band, a tier label).

## 8. Open / data-gated (what would raise resolution next)
- **Subject-specific in-vivo knee** (Grand Challenge) — fetch underway; would resolve the cross-cohort + JAM
  donor-geometry confounds at once.
- Clean **OrthoLoad hip waveform** (current copy corrupt) — re-fetch underway.
- The **residuals** (knee ~21% push-off, hip 9–19%) — attributable to l_TS + rigid-cut + structural, not yet
  cleanly split.
- Subject-specific muscle volumes (currently a 3rd-party MRI-PCSA transcription), n=3 subjects, SO-only for
  most (CMC only for subj2), no subject-specific cartilage/foot models.

---
*Confidence tiers per claim live in the individual docs and commit messages. This synthesis adds no new
measurement — it connects the measured ones. `MECHANISM_TRUST_LEDGER.md` remains the row-level ledger; this is
the force-investigation narrative over it.*
