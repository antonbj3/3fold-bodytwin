# MECHANISM BONE STRESS CONSEQUENCE — does the 1.4–2.0x joint-contact-force over-prediction propagate into inflated bone stress? (2026-07-21)

**Question.** `docs/MECHANISM_CROSS_SUBJECT.md` measured the twin's hip/knee joint **contact** force
over-predicting the OrthoLoad in-vivo anchor by **1.412–1.990×** (subject2's own JR-convention ratios:
hip 1.4129×, knee 1.5147×; cross-subject range across subjects 2–4: 1.412–1.990×). `scripts/msk/bone_stress.py`
consumes exactly those two joint-contact-force vectors (labels `hip_joint_contact`/`knee_joint_contact` in its
own free body) plus 26 femur-attaching muscle forces to predict peak femur cortical stress. **If the input is
over-predicted, is the output stress correspondingly inflated — possibly past physiological bounds?**

**Headline answer: NO, not coherently — and the reason is mechanistic, not a modeling artifact.** Across
**46 distinct configurations** (8 force-scale variants × primary geometry; the AS-IS/matched/1.4×/2.0× force
variants swept across all 5 literature offsets; AS-IS/matched swept across all 9 literature cross-sections),
**every single one** predicts peak femur strain **at or above 1573 µε** — comfortably above the verified
400–800 µε typical-walking band (Lanyon 1975, Burr 1996 — Sec. ANCHORS) — and the AS-IS (inflated-input)
configuration already sits **inside** the task's stated 2500–4000 µε damage-threshold band (3233 µε). Critically,
**de-inflating the joint contact force by subject2's own measured ratio does not bring stress down toward
physiological levels — it makes it WORSE** (3233→3493 µε, **+8.0%**) at the primary geometry, and this direction
holds at **all 9 swept cross-sections** and at **3 of 5 swept neck-shaft offsets** (0/30/40mm; the direction
flips modestly at 50/60mm — Sec. 3.3). Full ablation of both joint contact forces (holding muscles fixed)
roughly **doubles** stress (3233→7559–7631 µε). An exact, closed-form decomposition (Sec. 4) shows why: at
the peak section, the knee's own large individual bending contribution (189.2 N·m) is **nearly cancelled** by
the muscle+patellofemoral contribution (192.5 N·m combined, opposite sign) — de-inflating knee's force
**shrinks the cancelling term, not the thing being cancelled**, so the net bending moment at that location
*grows*. **The elevated bone-stress finding is real but is decoupled from the specific 1.4–2.0× joint-contact
over-prediction ratio** — it traces to the (already-flagged, literature-swept, not subject-measured) 40mm
neck-shaft offset assumption and to the muscle-force magnitudes, which this counterfactual does not (and, per
Sec. 6, structurally cannot cleanly) correct.

| | AS-IS (inflated input) | de-inflated (subject2's own ratio) | falsifier check |
|---|---:|---:|---|
| Peak stress, primary config | −54.97 MPa | −59.37 MPa | de-inflation **increased** magnitude |
| Peak microstrain | 3233 µε | 3493 µε | **+8.0%**, not a drop |
| Classification (task bands) | AT/ABOVE damage-threshold band | AT/ABOVE damage-threshold band | unchanged classification |
| Min µε across all 46 configs | 1573 µε (2.0× the 800 µε typical-walking edge) | — | **falsifies** "coherent propagation" |

---

## 1. Method — reuse discipline, and what is genuinely new

`scripts/msk/bone_stress_consequence.py` imports `scripts/msk/bone_stress.py`'s own geometry/free-body/
section functions **unedited** (`femur_axis_geometry`, `frame_loads_at`, `internal_force_moment`,
`combined_stress`, `region_label`, `pure_sided_bend_magnitudes`, `run_synthetic_selftests` — all module-level,
imported directly, never copy-pasted) and `validate_joint_force.py`'s `parse_mot`/Savitzky-Golay machinery
(same functions `bone_stress.py` itself uses). md5 of all three reused scripts is **identical before and
after** this session's runs (Sec. 5) — zero bytes edited.

**The counterfactual.** `bone_stress.py`'s own free body labels its two OrthoLoad-validated loads
`"hip_joint_contact"` and `"knee_joint_contact"` — these ARE, by direct code inspection, the exact
JointReaction-derived vectors whose peak magnitudes (387.04 %BW hip, 391.11 %BW knee) are the numbers
`docs/MECHANISM_CROSS_SUBJECT.md` compared against OrthoLoad. This script scales **only** those two labels'
force vectors by `1/ratio` at every frame, holding all 26 femur-attaching muscle forces, the patellofemoral
joint force, segment weight, and translational d'Alembert inertia **unchanged** — the literal, minimal
operationalization of "what if the joint contact force were not over-predicted."

**What is new** (cannot pre-exist, since `bone_stress_results.json` stores only swept peaks, not raw per-frame
loads): (1) re-running the identical per-frame free-body loop to recover raw per-frame loads; (2) the
scale/ablation counterfactual itself; (3) a vectorized, numpy-batch reimplementation of the method-of-sections
formula (`internal_force_moment_batch`) — derived directly from bilinearity of the cross product
(`cross(p−x_s, f) = cross(p,f) − cross(x_s,f)`, splitting each load's moment into an s-independent
precomputed term and a rank-1 s-dependent correction) so that 46 full offset/geometry sweeps run in ~4 seconds
instead of ~80 seconds each; verified point-by-point against the reused, unedited `bs.internal_force_moment`
to float precision (max diff 1.1e-12, Sec. 5) before being trusted for anything; (4) an **exact closed-form
decomposition** of the internal (N, M_bend) at a fixed section into additive per-load-category contributions
(hip / knee / patellofemoral / muscles / weight / inertia) — exact because the distal/proximal-averaged
method-of-sections sum is a **linear functional** of the load set (each load's own distal-vs-proximal side
assignment doesn't depend on what else exists, so per-category sums add to the total by linearity — proved
algebraically, then verified numerically as a hard gate, diff 1.6e-12, Sec. 5).

**Pre-registered thresholds** (stated before any new number was computed):
- *"coherent downstream propagation"*: de-inflated microstrain drops by a factor comparable to the applied
  force ratio **and** moves from the damage-threshold band toward the typical-walking band.
- *"physiological, decoupled from the inflation"*: de-inflated microstrain lands in 400–800 µε while AS-IS
  does not.
- **Falsifier for "coherent propagation"**: AS-IS already inside the typical-walking band (it is not — 3233
  µε), OR ablating the joint contact force entirely leaves stress roughly intact or lower (it does not — it
  roughly **doubles**, Sec. 3.2). **This falsifier fired.**

Neither pre-registered bucket cleanly fits what was actually measured — the finding is a third, better-evidenced
outcome (Sec. 3–4), reported as such rather than forced into either bucket.

---

## 2. Verification (machine-checked before any downstream number was trusted)

| gate | result |
|---|---|
| Ratio literals (1.4129×, 1.5147×, subject3's 1.745×/1.990×) match the LIVE `cross_subject_summary.json` | PASS |
| Reused `bs.run_synthetic_selftests()` (bone_stress.py's own beam-sectioning closed-form tests) | PASS |
| This script's own hip/knee %BW re-derivation vs the pre-existing validated certs | PASS (387.03722/391.11478 %BW, diff <1e-5) |
| `internal_force_moment_batch` vs reused `bs.internal_force_moment`, 130 sampled (frame,s) points | PASS, max diff 1.1e-12 |
| AS-IS primary peak reproduces `bone_stress_results.json` bit-for-bit | PASS, diff 4.26e-14 MPa |
| AS-IS full 5-point offset sweep reproduces upstream bit-for-bit | PASS, all 5 |
| Decomposition sums to total (N, M_bend) | PASS, diff 1.6e-12 / 5.7e-14 |
| Linearity gate: hip/knee category contribution scales exactly by hip_scale/knee_scale | PASS, diff <2.3e-13 |
| Geometry-scaling orthogonal to force-scaling (symmetric QC, Sec. 4.2) | PASS, 0.33% relative (see honest note below) |

**All 11 gates PASS** (`data/.../bone_stress_consequence_results.json`, key `gates_all_pass: true`).

One gate needed a mid-session correction, disclosed rather than hidden: a first pass gated geometry/force
orthogonality at an exact (<1e-6 absolute) tolerance and it FAILED (diff 6.84e-3, i.e. AS-IS's geometry-sweep
max/min ratio 2.0889 vs matched's 2.0958). Forced Orient before accepting either a false pass or a false
alarm: geometry (A, I, c) enters **only** the final stress conversion, never the (N, M_bend) computation
itself — exact orthogonality holds *at a fixed (frame, s)*, but the reported peak is a `max()` over a
discrete (frame, s) search whose ranking metric depends on geometry, and `max()` is a nonlinear operator.
Confirmed directly: AS-IS's peak lands at the identical (s, t) for all 9 swept geometries, while matched's
peak lands at a *different* (s, t) for exactly one (the smallest, od24×th8mm) of the 9 — a genuine, tiny,
expected discrete argmax shift, not a confound. Re-registered the gate at a realistic 5% relative tolerance
(actual: 0.33%) — PASS. This is exactly the kind of self-audit the method requires: the first tolerance was
wrong, not the physics.

---

## 3. Results

### 3.1 Primary configuration (offset=40mm, OD=27mm/th=6.5mm — `bone_stress.py`'s own primary point)

| variant | hip scale | knee scale | peak stress (MPa) | peak µε | classification | location |
|---|---:|---:|---:|---:|---|---|
| AS-IS (no correction) | 1.000 | 1.000 | −54.97 | 3233 | AT/ABOVE damage-threshold band | s=0.156, t=0.56s |
| matched (subject2's own ratio) | 0.708 | 0.660 | **−59.37** | **3493** | AT/ABOVE damage-threshold band | s=0.150, t=0.52s |
| uniform 1.4× de-inflation | 0.714 | 0.714 | −55.73 | 3278 | AT/ABOVE damage-threshold band | s=0.156, t=0.53s |
| uniform 2.0× de-inflation | 0.500 | 0.500 | −71.12 | 4184 | ABOVE damage-threshold upper bound | s=0.150, t=0.52s |
| worst measured cross-subject ratio (subject3's, applied to subject2's own geometry) | 0.573 | 0.503 | −71.76 | 4221 | ABOVE damage-threshold upper bound | s=0.150, t=0.52s |
| ablation: zero BOTH contact forces | 0 | 0 | **+128.50** | 7559 | ABOVE damage-threshold upper bound | s=0.150, t=0.51s |
| ablation: zero hip only | 0 | 1.0 | −103.13 | 6067 | ABOVE damage-threshold upper bound | s=0.850 ("distal shaft"), t=0.08s |
| ablation: zero knee only | 1.0 | 0 | −129.73 | 7631 | ABOVE damage-threshold upper bound | s=0.150, t=0.51s |

Every de-inflation/ablation variant moved stress **away from**, not toward, the 400–800 µε physiological band.
The monotonic trend (more aggressive "correction" → higher predicted stress: 3233→3278→3493→4184→4221→
6067–7631 µε) is the opposite of the naive "fix the inflated input, the output improves" expectation.

### 3.2 Ablation isolates the driver: muscles + geometry, not joint contact force

Zeroing **both** contact-force terms (hip AND knee) while holding the 26 muscle forces, patellofemoral force,
weight, and inertia fixed gives **+128.50 MPa (7559 µε)** — roughly **2.3×** the AS-IS magnitude, not a
reduction. This directly falsifies "the joint contact force is what's driving the elevated stress"; it
supports "muscle forces + the geometric offset assumption dominate, and the joint contact force partially
*counteracts* rather than adds to the bending" (mechanism in Sec. 4).

### 3.3 Offset sweep — the direction is not universal, reported honestly

| offset | AS-IS (MPa / µε) | matched (MPa / µε) | matched > AS-IS (magnitude)? |
|---:|---:|---:|:---:|
| 0mm (null eccentricity) | −26.75 / 1573 | −52.53 / 3090 | **True** (+96.4%, the largest swing found anywhere) |
| 30mm | −47.19 / 2776 | −56.86 / 3345 | **True** (+20.5%) |
| 40mm (primary) | −54.97 / 3233 | −59.37 / 3493 | **True** (+8.0%) |
| 50mm | −62.82 / 3695 | −62.56 / 3680 | False (−0.4%, essentially flat) |
| 60mm | −70.72 / 4160 | −66.19 / 3893 | False (−6.4%) |

The sign flips between 40mm and 50mm: at large assumed offsets, hip's own eccentricity-scaled moment
contribution grows large enough that de-inflating hip's force finally wins out over the knee-muscle
cancellation-breaking effect (Sec. 4). **Even in the direction that helps, the effect is small** (−6.4% at
best, 60mm) and never approaches the 4–5× reduction that would be needed to reach the 400–800 µε band from
any of these starting points (2460–4814 µε at minimum across the full geometry×offset grid, Sec. 3.4).

### 3.4 Geometry sweep — direction confirmed across all 9 literature cross-sections (symmetric QC)

At every one of the 9 swept cross-sections (24/27/30mm OD × 5/6.5/8mm thickness), matched > AS-IS in
magnitude (checked directly, not assumed): od24×th5 (4814→5208 µε) through od30×th8 (2304→2485 µε). **The
counter-intuitive direction is a property of the force decomposition, not an artifact of which literature
cross-section happens to be selected** — this is the direct answer to the task's own flagged symmetric-QC
concern ("a wrong bone geometry could mask or manufacture the effect independent of the force"): the
geometry is the SAME literature-swept range `bone_stress.py` already discloses as unmeasured (Sec. 6), it is
not being cherry-picked here, and the qualitative finding is invariant to it.

---

## 4. The mechanism — an exact decomposition, not a guess

At the AS-IS primary peak (s=0.1556, t=0.56s, subtrochanteric shaft), the internal bending-moment vector
splits **exactly** (verified, Sec. 2) into:

| load category | axial N (N) | \|M_bend\| (N·m) | M_bend, X-component | M_bend, Z-component |
|---|---:|---:|---:|---:|
| hip_joint_contact | −1443.78 | 35.74 | −35.61 | +3.02 |
| knee_joint_contact | −1237.16 | **189.21** | −26.15 | **−187.39** |
| patellofemoral_joint_contact | −8.31 | 41.84 | +0.61 | +41.84 |
| muscles (26, summed) | +1329.27 | 150.70 | −32.19 | **+147.23** |
| femur segment weight | +46.28 | 1.16 | +0.44 | −1.07 |
| d'Alembert inertia (translational) | +3.41 | 1.78 | +0.13 | +1.78 |
| **TOTAL** | **−1310.29** | **92.93** | **−92.77** | **+5.40** |

Two different things are happening in the two bending components, and averaging them into one scalar
`|M_bend|` per category obscures it:
- **X-component** (which dominates the 92.93 N·m total): hip (−35.6), knee (−26.2), and muscles (−32.2) are
  all the **same sign** and add up coherently. This is where the joint contact force genuinely contributes to
  bending, as expected.
- **Z-component**: knee's own contribution is **large** (−187.4 N·m) but is **almost exactly cancelled** by
  muscles + patellofemoral (+147.2 + 41.8 = +189.1 N·m, opposite sign) — net Z is only +5.4 N·m. De-inflating
  knee's force by its measured ratio (×0.660) shrinks knee's Z-contribution to −123.7 N·m, but muscles +
  patellofemoral's Z-contribution is **untouched** (they are not scaled in this counterfactual) — so the
  cancellation breaks, and net Z jumps to +68.2 N·m, more than 12× larger. This is what drives the
  matched-config's total `|M_bend|` up despite hip+knee individually shrinking (verified exactly: the
  linearity gate confirms each category's own contribution scales exactly by its own ratio — the surprise is
  in how the vector **sum** behaves, not in any individual term).

This is a direct, quantified illustration of the Moissenet mechanism already named in
`docs/MECHANISM_JOINT_FORCE_SCORECARD.md` (re-verified live this session, PMID 24210475: *"joint reaction
forces are typically estimated using a two-step method, computing first the musculo-tendon forces by [SO]
and then deducing the joint reaction forces from the force equilibrium... the joint reaction forces are
usually overestimated"*). The joint contact force is **not an independent input** to this beam model in the
way a clean "fix the input, the output improves" intuition assumes — it is entangled, via the free body's own
force/moment balance, with the (here, unrescaled) muscle forces. Scaling one term of an equilibrium down while
holding the others fixed does not, in general, shrink the net result; whether it does depends on which
direction that term's own contribution was pointing relative to everything else, which is exactly what the
decomposition (not a guess) reveals.

A single scalar "contact-force fraction of total bending" was also computed (0.548, i.e. hip+knee project to
about 55% of the total M_bend's own direction) but is flagged as **not** a reliable predictor of scaling
direction — it is a projection-based summary of a single vector sum, and (as shown above) it is the *vector*
cancellation structure, not this scalar fraction, that determines whether de-inflating the contact force
raises or lowers the net result.

---

## 5. Cross-check with a concurrent sibling finding (informational, not independently audited here)

A concurrent sibling investigation this session (`data/msk_smoketest/subject2_walking1/
metabolic_calorimetry_anchor/metabolic_calorimetry_anchor_results.json`, spot-checked read-only, not
re-derived or audited by this analysis) reports a metabolic heat-rate model over-predicting a net-cost-of-
transport calorimetry anchor (Koelewijn et al. 2019, PMID 31532796) by a mean ratio of **1.79×** (range
1.55–2.11× across 4 model combinations) — comparable in order of magnitude to this repo's own force
over-prediction ratios (mean 1.54×, the same `knee_so=1.5146`/`hip_so=1.4119` etc. figures this document
also uses), and that sibling analysis's own gate reads `h_broad_systemic_over_prediction_supported: True`.

This is **consistent with, not proof of**, a shared upstream cause: if Static Optimization's muscle-force
solution is itself systemically elevated (which would inflate both metabolic heat rate, which scales with
muscle activation, AND joint contact force, via the Moissenet mechanism above), then this document's own
finding — that muscle-force contributions (never scaled down in this counterfactual) dominate and sustain
the elevated bone-stress prediction regardless of how the *derived* joint contact force is corrected — is
exactly what that shared-cause picture would predict. **This is not tested here** (would require re-scaling
the muscle-force term itself, a different and larger counterfactual than this document's scope) and is
flagged as the concrete next step (Sec. 7), not claimed as resolved.

---

## 6. Symmetric QC on geometry (the task's own flagged concern, answered directly)

The task asks: *"a wrong bone geometry could mask or manufacture the effect independent of the force — verify
the geometry is anchored (cite) before attributing stress level to the force."* Answered in two parts:

1. **Is the geometry anchored?** No, and this is not new to this document — `bone_stress.py`'s own honest
   gaps (inherited unchanged, `docs/MECHANISM_BONE_STRESS.md` Sec. 2/6) already disclose the 27mm OD / 6.5mm
   thickness circular annulus and the 40mm neck-shaft offset as a **literature-swept range** (boneandspine.com
   clinical reference; THA-literature "native offset" ~40–44mm target), not this subject's own CT/DXA-measured
   geometry. This document does not change or improve that — it inherits the gap.
2. **Does that unanchored geometry choice manufacture the de-inflation-increases-stress finding?** No —
   verified directly (Sec. 3.4): the finding holds at **all 9** swept cross-sections, and the geometry-vs-
   force-scale orthogonality gate (Sec. 2) confirms geometry and force-scaling act independently to within
   0.33% (the residual explained exactly by a discrete argmax shift, not a systematic interaction). So while
   the **absolute** stress numbers in this document inherit the upstream geometry uncertainty (as they must),
   the **qualitative, decisive finding** — de-inflation does not bring stress toward the physiological band —
   is not an artifact of which point in that literature range was chosen.

---

## 7. Anchors — live-verified this session (NCBI eutils + Europe PMC REST), not recalled

| citation | PMID / DOI | verified content | used as |
|---|---|---|---|
| Lanyon LE, Hampson WG, Goodship AE, Shah JS. "Bone deformation recorded in vivo from strain gauges attached to the human tibial shaft." *Acta Orthop Scand.* 1975;46(2):256-68. | 1146518 / 10.3109/17453677508989216 | Abstract (fetched verbatim): WALKING peak strain **~400 µε** (compression, before toe-off); RUNNING peak **~850 µε** (tension). First human in-vivo tibial strain-gauge study. | Primary anchor for the 400–800 µε "typical walking" band (its own low/high ends). |
| Burr DB, Milgrom C, Fyhrie D, et al. "In vivo measurement of human tibial strains during vigorous activity." *Bone.* 1996;18(5):405-10. | 8739897 / 10.1016/8756-3282(96)00028-2 | Abstract (verbatim): even **vigorous** activity stayed <2000 µε, "about three times higher than... walking" (implies walking ~650–700 µε); *"strains in the human tibia during walking and running are well below the fracture threshold."* | Corroborates the walking band; corroborates the damage/fracture ordering. |
| Milgrom C et al. "Do high impact exercises produce higher tibial strains than running?" *Br J Sports Med.* 2000;34(3):195-9. | 10854019 / 10.1136/bjsm.34.3.195 | Verified abstract; running/impact > walking strain ordering. | Directional corroboration. |
| Milgrom C et al. "In-vivo strain measurements to evaluate the strengthening potential of exercises on the tibial bone." *J Bone Joint Surg Br.* 2000;82(4):591-4. | 10855890 / 10.1302/0301-620x.82b4.9677 | Verified abstract: *"walking, an activity which has been found to have only a minimal effect on bone mass."* | Confirms walking as the low-strain reference activity in this literature. |
| Cezayirlioglu H, Bahniuk E, Davy DT, Heiple KG. "Anisotropic yield behavior of bone under combined axial force and torque." *J Biomech.* 1985;18(1):61-9. | 3980489 / 10.1016/0021-9290(85)90045-4 | Verified (NCBI efetch + Europe PMC core, both queried) to exist and to study **human AND bovine tibial AND femoral** cortical bone yield loci under combined tension/compression/torsion — the topically right reference class for a "cortical yield stress" figure — but neither source's abstract restates an exact MPa number. | Confirms the RIGHT reference class exists for the task's stated 150–180 MPa cortical yield band; does **not** independently re-verify the exact figure (same honest-gap tier as `bone_stress.py`'s own unresolved Reilly & Burstein 1975 E=17 GPa citation, PMID 1206042 — re-checked again this session via both NCBI efetch and Europe PMC core, still no abstract on file anywhere). |
| Frost HM 1987 (PMID 3688455) | reused, not re-fetched | 1500–3000 µε formation-onset threshold (already live-verified in `docs/MECHANISM_BONE_REMODELING.md`) sits *below* the task's stated 2500–4000 µε "damage threshold" band. | Confirms a formation < damage < fracture ordering; the task's specific 2500–4000 figure itself was not independently pinned to an exact primary source this session (same tier as that document's own already-flagged "4th pathologic-overload bin" hedge). |
| Moissenet F, Chèze L, Dumas R. "A 3D lower limb musculoskeletal model..." *J Biomech.* 2014;47(1):50-8. | 24210475 / 10.1016/j.jbiomech.2013.10.015 | Re-verified live this session (already cited in `docs/MECHANISM_JOINT_FORCE_SCORECARD.md`): *"the joint reaction forces are usually overestimated"* by the two-step (SO-then-subtract) method architecture. | The named mechanism Sec. 4's decomposition illustrates directly with real numbers. |

**Cross-bone honest gap, stated plainly**: every walking-strain anchor above is a **tibial** in-vivo
strain-gauge study (the standard human in-vivo site, since the tibia is subcutaneously accessible; genuine
human in-vivo **femoral** cortical strain-gauge data during unassisted gait essentially does not exist in the
literature, since it would require surgical, non-prosthesis-linked gauge placement). This document's
prediction is **femoral** (subtrochanteric/diaphyseal). Frost's mechanostat thresholds (already used upstream,
reused here) are bone-tissue-generic, not tibia-specific, which partially — not fully — justifies the
cross-bone comparison. Flagged, not hidden.

**Consistency cross-check (internal, both unit systems agree):** converting the task's 150–180 MPa yield band
through the SAME E=17 GPa linear-elastic conversion this document uses gives ≈8824–10588 µε — well above
every microstrain value this document computed (max 7631 µε, the knee-ablation case) — consistent with this
document's own direct MPa comparison (max computed magnitude 129.73 MPa, still below 150 MPa). **No tested
configuration reaches yield in either unit system**, even though several exceed the (µε-domain) damage/
microdamage-risk threshold. These are different physiological failure modes (cumulative fatigue-microdamage
risk under repeated cyclic loading vs. acute overt yield from one loading event) and are not interchangeable —
reported as two distinct comparisons, not conflated into one.

---

## 8. Confidence tiers and the falsifier, stated explicitly

- **Claim**: this beam model's predicted femoral bone stress during gait is robustly above the physiological
  typical-walking strain range, and this elevation is decoupled from (not a coherent downstream propagation
  of) the specific 1.4–2.0× joint-contact-force over-prediction.
- **Confidence tier**: **in-vivo-anchored** for the external physiological bands themselves (Lanyon 1975,
  Burr 1996 — live PMID-verified, primary source, this session) crossed with **method-only** for the beam
  model's own absolute stress level (inherits every honest gap `bone_stress.py` already disclosed: generic,
  not CT/DXA-measured, cross-section and neck-shaft offset; circular-annulus idealization; uncorrected
  rotational inertia; bending-dominated result in tension with the verified Taylor 1996 anchor, PMID 8673318).
  The **decoupling** claim itself (Sec. 3–4) is machine-verified to float precision (an internal, closed-form
  property of this specific model), independent of whether the model's absolute MPa level is itself accurate.
- **Falsifier, pre-registered and fired**: "coherent propagation" required de-inflated stress to drop
  proportionally to the applied ratio and move toward the physiological band. Measured: it moved further away
  (+8.0% at the primary config; up to +96.4% at zero offset; full ablation +134–136%), and the global minimum
  across all 46 tested configurations (1573 µε) never came within 2× of the 800 µε physiological-band edge.
  **This falsifies coherent propagation as the explanation.**
- **What would falsify the decomposition's own causal story** (Sec. 4): if the same Z-component
  near-cancellation did NOT reproduce with an independently-recomputed free body, or if the linearity gates
  failed. Neither happened (Sec. 2) — checked, not assumed.

---

## 9. Honest gaps

1. **FIXED-MUSCLE-FORCE counterfactual, not a re-solved simulation.** If Static Optimization's own muscle
   recruitment is itself a root cause of the contact-force over-prediction (the Moissenet mechanism, Sec. 4),
   the 26 muscle forces used unchanged in every de-inflated variant here still carry that same unremoved
   bias. This counterfactual isolates the joint-CONTACT-force-specific contribution only; it does not
   simulate "what SO would have produced with a corrected contact force." Given the concurrent sibling
   metabolic finding (Sec. 5, informational), the true physiologically-corrected bone stress plausibly sits
   lower than anything computed here — untested, flagged as the concrete next step.
2. **Uniform-ratio-across-the-whole-trial assumption.** The measured 1.4–2.0× ratio is a peak-vs-peak
   comparison; applying it as a constant scalar to the entire per-frame force trajectory is a simplifying
   assumption, not independently verified to hold at every gait-cycle phase.
3. **Cross-bone comparison** (Sec. 7): tibial in-vivo strain anchors vs. a femoral model prediction.
4. **Inherits every upstream `bone_stress.py`/`bone_remodeling.py` honest gap unchanged**: generic (not
   subject-measured) cross-section and offset; circular-annulus idealization; only translational (not
   rotational) inertia corrected; the bending-dominated-vs-Taylor-1996-compression-dominated tension.
5. **The 2500–4000 µε "damage threshold" and 150–180 MPa yield-stress bands are the task's own stated
   figures.** Live verification this session confirmed the right reference class for each (Frost 1987's own
   1500–3000 µε formation threshold sits below the damage band, a consistent ordering; Cezayirlioglu 1985
   studies exactly human tibial+femoral cortical yield behavior) but did not independently re-extract the
   exact numeric figures from primary text — reported as citation-exists/topically-correct, not as an
   independently re-verified precise number (same tier as the upstream E=17 GPa gap).
6. **The concurrent sibling metabolic-calorimetry finding (Sec. 5) was spot-checked, not audited.** Its own
   internal numbers were read directly from its JSON and are internally consistent with the hip/knee ratios
   this document also uses, but its own literature-verification rigor was not independently re-derived here.
7. subject2/walking1, right leg, single trial only — same scope boundary as every MSK cert in this repo family.

---

## 10. Files

- `scripts/msk/bone_stress_consequence.py` (new) — the full counterfactual pipeline: re-derives the per-frame
  free body (verified bit-for-bit against `bone_stress.py`'s own output), the vectorized batch method-of-
  sections solver (verified against the reused, unedited scalar function), the 8-variant scale/ablation sweep,
  the offset and geometry sweeps, the exact closed-form decomposition, and the literature-band classification.
- `data/msk_smoketest/subject2_walking1/bone_stress_consequence/bone_stress_consequence_results.json` — every
  number in this document, machine-written, `gates_all_pass: true`.
- **Read-only, not modified** (md5-verified identical before/after, Sec. 2): `scripts/msk/bone_stress.py`
  (90cc3ffa6c601d1dffcef776779c5727), `scripts/msk/bone_remodeling.py` (3481c9b029f40d7a4cfa54b743a5f2dc),
  `scripts/msk/validate_joint_force.py` (e2bb15958cfe8d72f44a3e0ac3d95950).
- **Read-only inputs**: `data/msk_smoketest/subject2_walking1/bone_stress/bone_stress_results.json`,
  `data/msk_smoketest/cross_subject_validation/cross_subject_summary.json`,
  `data/msk_smoketest/subject2_walking1/static_optimization/{jr,so}/*.sto`.
- **Spot-checked, not depended on**: `data/msk_smoketest/subject2_walking1/metabolic_calorimetry_anchor/
  metabolic_calorimetry_anchor_results.json` (concurrent sibling instance's own work, Sec. 5).
- No git commit, no git push performed (isolation respected). No file outside this task's own new outputs was
  touched.
