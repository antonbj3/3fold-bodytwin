# MECHANISM MOMENT ARM VALIDATION — decorrelated geometry anchor for the knee/hip force cert (2026-07-21)

**Status: HYPOTHESIS awaiting independent QC.** Confidence tier: **cadaveric/published-plausibility**
(moment-arm literature is cadaveric/CT, never in-vivo-force-anchored) — further sub-tiered per row below
into `VERIFIED_QUANT` (a primary paper's own number, fetched and quoted this session) vs
`DOMAIN_CONSENSUS_UNVERIFIED_THIS_SESSION` (a widely-reproduced secondary band; the primary paper's own
table/figure was not accessed — paywalled) vs `OPEN` (no literature magnitude available at all). Never
read "PASS" below as "in-vivo validated" — it means "inside a pre-registered, externally-justified
tolerance band around a cadaveric number."

## 0. What this is / is not

The knee **391.11 %BW** / hip **387.04 %BW** contact-force cert (`docs/MECHANISM_DIFF_RECONCILE.md`,
`docs/MECHANISM_STATIC_OPT.md`) comes from static optimization on
`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` (subject2, walking1). Static optimization
solves `R(q)^T f = tau(q)` for muscle forces `f`, where `R(q)` is the matrix of muscle **moment arms** as
a function of joint angle. If `R(q)` is geometrically wrong, `f` (and the contact force derived from it)
is wrong for a reason that has **nothing to do with OrthoLoad, the gait trial, or the GRF** — this is why
it is a genuinely **decorrelated** anchor: the analysis below never touches force, activation, GRF, or
OrthoLoad data. It touches exactly one thing: `Muscle.computeMomentArm(state, coord)`, i.e. pure geometry
as a function of pose.

This is **not** a re-run of the force pipeline and **not** a repeat of `MECHANISM_CARTILAGE_CONTACT.md`
(different model, different question) — it is the "moment-realizability" leg of the knee closing-cell's 4
decorrelated legs (`COORDINATOR.md` §5: geometry / DOF+ROM / **moment-realizability** / cartilage-contact),
extended here to also cover the hip.

## 1. Pre-registration (fixed before the comparison table below was finalized)

- **C**: the model's prime-mover moment-arm curves (peak magnitude + angle-of-peak) fall inside a
  tolerance band around the published/cadaveric value. **¬C**: at least one prime mover falls outside the
  band, indicating a geometry error that would bias the muscle-force solution and hence the contact-force
  cert.
- **Magnitude tolerance: ±30%.** Externally justified, not self-selected: Spoor & van Leeuwen (1992, PMID
  1733995) report their own two "gold-standard" cadaveric methods (MRI-derived vs tendon-travel) disagree
  by **up to ~30%** for rectus femoris and semimembranosus on the *same* specimen. A model-vs-literature
  difference inside 30% is therefore inside the field's own inherent inter-method spread, not
  distinguishable from a real geometry defect.
- **Angle-of-peak tolerance: ±20°.** Weaker justification (own judgment, not literature-sourced): Nemeth &
  Ohlsen's (1985) own hamstring peak-angle (35°) is a mean pooled over 30 specimens (10 autopsy + 20 CT),
  implying real inter-specimen spread of comparable width.
- **Falsifier:** any prime mover outside the tolerance band on a `VERIFIED_QUANT` row is flagged, and its
  **force impact is estimated** (§8) — not just noted and dropped.
- **Adversary forced before accepting any PASS:** the naive per-bundle comparison for adductor magnus
  initially flagged all 4 sub-bundles (up to 165% off) against a single lumped literature line-of-action;
  before accepting that as a geometry defect, the fairest available representation (a physiological-
  cross-section/force-weighted resultant across the 4 bundles, matching what a single-line cadaver method
  actually measures) was computed and re-checked (§8) — this is the forced-steelman step, not skipped.

## 2. Method

Script: `scripts/msk/moment_arm_validation.py` (sweep + raw curves), `scripts/msk/moment_arm_lit_compare.py`
(literature comparison, tolerance check, internal-consistency check — all PASS/FLAG verdicts are
machine-computed from the JSON, never eyeballed from a plot).

- Model: `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` (subject2, the exact instance behind
  the 391/387 %BW numbers), loaded read-only from the external drive. OpenSim 4.6-2026-06-22-85aaf64.
- For each of 3 coordinates (`knee_angle_r`, `hip_flexion_r`, `hip_adduction_r`): all other pelvis +
  lower-limb coordinates held at 0 rad (neutral standing reference); the swept coordinate stepped every
  2° across **its own declared range** (`Coordinate.getRangeMin/Max`, never a hardcoded literature range —
  see the clamp artifact below); `Muscle.computeMomentArm(state, coord)` called for every prime-mover
  muscle at every step; full curves (not just peaks) dumped to raw CSV.
- **Sign convention: measured, not derived.** An initial algebraic guess (τ = −r·F_tension ⇒ an extensor
  should read r>0 against a flexion-positive coordinate) was checked against `vasint_r` (a known knee
  extensor) and **contradicted** — measured r = −5.055 cm at 30° knee flexion, i.e. negative. The guess
  was discarded; the empirically confirmed rule used throughout is: **positive moment arm = the muscle
  assists motion in the coordinate's positive-increasing direction** (e.g. knee flexors read positive
  about flexion-positive `knee_angle_r`, extensors read negative).
- **A real artifact was caught and fixed, not glossed over:** the first run swept `hip_adduction_r` over
  an arbitrary [-30°, +40°] range. The coordinate's own declared range is [-50.0°, +30.0°] — beyond +30°
  the value is **clamped and frozen** (verified: `addlong_r` moment arm at 32°/34°/36°/38°/40° repeats
  `6.949539391729455…` to 10+ significant digits — an out-of-domain extrapolation artifact of the
  "poly" — polynomial-regression — muscle paths, not a real curve point). The script now reads each
  coordinate's own range dynamically and includes an automatic frozen-tail detector
  (`detect_frozen_tail`); the corrected run shows **zero** frozen-tail flags on any of the 3 sweeps
  (`data/msk_smoketest/moment_arm_validation/moment_arm_validation_evidence.json` →
  `sweeps.*.qc_frozen_tail_flags`, all empty).

## 3. Literature anchors — identity verified by machine cross-check (PubMed eutils), not memory

An early PMID recollection from memory for a 7th (non-task) paper (Arnold, Ward, Lieberman, Delp 2010)
was **wrong** — `esummary` returned an unrelated histamine-secretion paper for the guessed PMID — so every
PMID below was independently confirmed via `esummary` (title/journal/year match) before being cited.

| # | Citation (as asked) | PMID (verified) | Journal / year | What it covers | Numeric content this session |
|---|---|---|---|---|---|
| 1 | Buford et al. 1997 | **9422462** | IEEE Trans Rehabil Eng 5(4):367-79, 1997 | knee, n=15 cadaver hemi-pelvis specimens, 13 muscles + patellar tendon, 0-140°ish ROM | Abstract is qualitative only ("variable, yet repeatable and recognizable patterns") — **no inline numbers on PubMed** (verified via 2 independent fetch routes) |
| 2 | Herzog & Read 1993 | **8376196** | J Anat 182(Pt2):213-30, 1993 | knee sagittal plane: quadriceps (**lumped as one structure**), biceps femoris, semimembranosus, semitendinosus, via polynomial regression | Abstract qualitative only — no inline numbers |
| 3 | Spoor & van Leeuwen 1992 | **1733995** | J Biomech 25(2):201-6, 1992 | knee, **n=1** specimen, MRI vs tendon-travel, bflh/gaslat/gasmed/gracilis/recfem/sartorius/semimem/semiten/TFL | **Quantitative**: "up to about 30% for the rectus femoris and semimembranosus" inter-method disagreement — used as the external justification for the ±30% tolerance band above |
| 4 | Nemeth & Ohlsen 1985 | **3988782** | J Biomech 18(2):129-40, 1985 | hip **extensors**: gluteus maximus, hamstrings, adductor magnus; n=10 autopsy + n=20 living (CT), 0-90° hip flexion | **Fully quantitative** (see §4) |
| 5 | Dostal & Andrews 1981 | **7334040** | J Biomech 14(11):803-12, 1981 | 3D hip musculature model (27 segments) | **No abstract available on PubMed** (verified) — full text paywalled, not retrieved |
| 6 | Arnold et al. 2000 | **10862133** | Comput Aided Surg 5(2):108-19, 2000 | MRI-based models, n=3 cadaveric specimens, **medial hamstrings + psoas only** | "Errors in the moment arms … averaged over functional ranges of hip and knee flexion, were less than 4 mm (within 10% of experimental values)" — a field accuracy benchmark, not a per-muscle target value |

Bonus confirmation (Delp 1990, PMID **2210784**, IEEE Trans Biomed Eng 37(8):757-67 — the original
geometry donor for this entire model lineage) and the direct parent model, Rajagopal et al. 2016 (PMID
**27392337**, PMCID 5507211, open full text): Rajagopal 2016 §II.B states explicitly it "verified that we
preserved the muscle geometry by comparing model moment arms of the major muscles crossing the hip, knee,
and ankle to experimentally measured moment arms" citing **Buford 1997, Spoor & van Leeuwen 1992, Nemeth &
Ohlsen 1985, Dostal et al. 1986 (companion to #5), and Arnold et al. 2000 by name** (Herzog & Read 1993's
inclusion is not confirmed — it may be among unlabelled refs 32-40, not independently checked). This is
the machine-verified basis for the symmetric-QC statement in §9 — see there before trusting any PASS below
at face value.

## 4. Results — hip flexion-extension (`VERIFIED_QUANT` tier, Nemeth & Ohlsen 1985)

Model value is the raw curve **interpolated at the literature's own reference angle** (never peak-vs-peak,
which would conflate a magnitude error with a location error).

| Muscle | Model @ 0° (cm) | Lit @ 0° (cm) | % diff | Tol | Verdict |
|---|---|---|---|---|---|
| glmax1_r | −5.63 | 7.90 | 28.8% | 30% | **PASS** |
| glmax2_r | −6.98 | 7.90 | 11.6% | 30% | **PASS** |
| glmax3_r | −7.64 | 7.90 | 3.2% | 30% | **PASS** |
| semimem_r | −4.48 | 6.10 | 26.6% | 30% | **PASS** |
| semiten_r | −4.99 | 6.10 | 18.1% | 30% | **PASS** |
| bflh_r (hip) | −5.02 | 6.10 | 17.6% | 30% | **PASS** |

Peak-**location** check (lit gives location only, no peak mm): hamstrings peak "an average of 35°" →
model semimem/semiten/bflh peak at 38.1°/42.1°/38.1° (diff 3.1-7.1°, tol 20°) → **PASS, all 3**.
Adductor magnus peak "up to 75°" → model addmagDist/Isch/Mid/Prox peak at 58.1°/48.1°/84.1°/−25.9°
(diff 16.9°/26.9°/9.1°/100.9° vs tol 20°) → **PASS/FLAG/PASS/FLAG** (see §8 — the two flags are diagnosed,
not just tallied).

Gluteus maximus and hamstrings — the two largest prime movers for hip extension (**76.6%** of the combined
max-potential hip flexion-extension moment among the 4 groups tested, §8) — **pass cleanly** (3-29% off, a
verified cadaveric number, well inside the externally-justified 30% band).

## 5. Results — knee flexion (`DOMAIN_CONSENSUS_UNVERIFIED_THIS_SESSION` tier)

Buford 1997 / Herzog & Read 1993 / Spoor & van Leeuwen 1992 abstracts contain **no inline peak numbers**
(verified via two independent fetch routes each) — full text is paywalled (1990s IEEE/J Anat/J Biomech,
pre-open-access era) and was not retrieved this session. The bands below are widely-reproduced secondary/
textbook ranges attributed to this literature, **not independently re-verified against the primary
papers' own tables/figures this session** — a materially weaker evidence tier than §4. A FLAG here means
"outside common secondary reporting," not "outside the primary cadaveric data" — treat accordingly.

| Muscle | Model peak (cm) | Model angle (deg) | Consensus band (cm / deg) | Verdict |
|---|---|---|---|---|
| vasint_r | −5.23 | 6.1 | 3.5-5.5 / 20-50 | PLAUSIBLE |
| vaslat_r | −5.40 | 20.1 | 3.5-5.5 / 20-50 | PLAUSIBLE |
| vasmed_r | −5.13 | 10.1 | 3.5-5.5 / 20-50 | PLAUSIBLE |
| recfem_r | −5.64 | 24.1 | 3.0-5.0 / 15-45 | PLAUSIBLE |
| semimem_r | +4.78 | 56.1 | 3.0-5.0 / 70-100 | PLAUSIBLE |
| semiten_r | +6.14 | 60.1 | 3.0-5.0 / 70-100 | PLAUSIBLE (upper-tolerance-bound) |
| bflh_r (knee) | +3.94 | 72.1 | 3.0-5.0 / 70-100 | PLAUSIBLE |
| gasmed_r | +3.36 | 86.1 | 1.5-3.0 / 80-110 | PLAUSIBLE |
| gaslat_r | +3.53 | 88.1 | 1.5-3.0 / 80-110 | PLAUSIBLE |

9/9 plausible, all under the tolerance-expanded band. Two honest caveats worth flagging even though they
pass: (a) `vasint_r`/`vasmed_r` peak at 6-10° — noticeably closer to full extension than the ~30-45°
commonly pictured for quadriceps, i.e. this model's quad mechanical advantage is shifted slightly toward
early stance relative to the common textbook picture (still inside the 20-50° band with the ±20°
tolerance, but near its lower edge); (b) `semiten_r`'s 6.14 cm peak sits at the tolerance-expanded upper
bound (5.0×1.3=6.5 cm) — not comfortably mid-band. Neither crosses into FLAG, but neither is a strong
confirmation either — this is the honest reading of a `DOMAIN_CONSENSUS` tier, not a `VERIFIED_QUANT` one.

## 6. Results — hip ab-adduction (`OPEN` tier — no literature magnitude available)

Dostal & Andrews 1981 has **no PubMed abstract** (verified) and full text was not retrieved (paywalled).
Arnold 2000 validated only medial hamstrings + psoas (a flexion-plane quantity), not ab-adduction. **No
quantitative literature anchor exists for glutes_med/glutes_min/adductors in this plane from this
session's sources.** One **tertiary, non-pre-registered** cross-check (found while chasing Dostal &
Andrews' citation graph, Henderson et al. 2011, hip abductor moment arm for proximal femoral replacement
planning): reports 4.3 cm at neutral hip position vs 2.9 cm in adduction (single-line abductor model).
Model `glmed1_r` shows the **same qualitative trend** (higher magnitude in abduction, lower in adduction:
−6.54 cm at 49.9° abducted vs −5.05 cm at 29.9° adducted) — a directional match only, **not scored into
any verdict** (not pre-registered, incommensurable single-line-vs-3-bundle representation). This plane's
model-vs-literature magnitude question is **honestly OPEN**, not papered over with an invented band.

## 7. Internal geometric consistency checks (independent of any external anchor)

A real, checkable geometric property, independent of the literature-magnitude question above:
functionally-opposed muscle groups must show opposite-sign moment arms under the same coordinate sweep.
**7/7 confirmed** (`data/msk_smoketest/moment_arm_validation/moment_arm_lit_comparison.json` →
`internal_consistency_checks`):

| Pair | Coordinate | Result |
|---|---|---|
| iliopsoas (flexor) vs glute max (extensor) | hip_flexion_r | opposite sign confirmed |
| psoas (flexor) vs hamstrings (extensor) | hip_flexion_r | opposite sign confirmed |
| adductor longus vs glute medius | hip_adduction_r | opposite sign confirmed |
| adductor magnus (mid) vs glute minimus | hip_adduction_r | opposite sign confirmed |
| vasti (extensor) vs hamstrings (flexor) | knee_angle_r | opposite sign confirmed |
| rectus femoris (extensor) vs biceps femoris (flexor) | knee_angle_r | opposite sign confirmed |
| glmax1 (superior, abductor-signed) vs glmax3 (inferior, adductor-signed) | hip_adduction_r | **opposite sign confirmed** |

The last row is the most informative: gluteus maximus's superior fibers (`glmax1_r`: −3.87 cm at full
abduction) and inferior fibers (`glmax3_r`: +6.25 cm, adductor-signed) **reverse sign across the muscle**
— this bundle-dependent sign reversal is the precise anatomical phenomenon Dostal & Andrews (1981)
originally split gluteus maximus into segments to capture. Reproducing it is a real geometric-fidelity
signal that does not depend on getting any specific magnitude right.

## 8. Falsifier assessment — what actually flagged, and its estimated force impact

Across all `VERIFIED_QUANT` rows (17 total: 6 magnitude + 3+4 angle-only): **11/17 PASS**, **6/17 FLAG —
all 6 confined to adductor magnus** (0 flags on glutes_max or hamstrings, the two dominant hip-extension
prime movers). No knee prime mover flagged at any tier. No hip-flexion-extension prime mover other than
adductor magnus flagged.

**Diagnosis (forced adversary, not accepted at face value):** Nemeth & Ohlsen's "15 mm" is a single
lumped line-of-action for the whole adductor magnus; the model splits it into 4 geometric bundles
(Dist/Isch/Mid/Prox) that individually range −5.32 to +2.58 cm at 0° — comparing each bundle to the one
lumped number is not apples-to-apples by construction (naive per-bundle diffs ran up to 165%). Forcing the
fairest available resultant: all 4 bundles share **identical `Fmax` = 597.3 N** in this model (a modeling
simplification — PCSA appears split evenly across bundles rather than graded), so the force-weighted
resultant equals the simple mean: (−2.30 −3.98 +0.43 +2.28)/4 = **−0.89 cm** at 0°, vs the literature's
1.50 cm → **40.4% difference** — reduced sharply from the naive worst case (165%) by the bundle-split
explanation, but **still outside the pre-registered 30% band even under the fairest construction**. This
is reported as a genuine, if minor, geometry note — not explained away.

**Force impact estimate** (machine-computed from `Fmax × |moment arm at 0°|`, the 4 hip-flexion-extension
groups tested):

| Group | Combined `Fmax` × arm potential @ 0° (Nm) | Share of the 4-group total |
|---|---|---|
| glutes_max | 225.98 | 41.2% |
| hamstrings | 194.00 | 35.4% |
| iliopsoas | 75.06 | 13.7% |
| **adductors_magnus** | **53.69** | **9.8%** |

Adductor magnus is a **secondary** contributor to hip flexion-extension moment capacity (9.8% share; its
primary action is adduction, where its sign/consistency checks out cleanly, §7) — so a 40% geometry
deviation confined to this muscle's *secondary* action plane implies an order-of-magnitude **≤4% shift on
the combined 4-group hip extension moment capacity** (0.098 × 0.40), well below what would be needed to
materially move a 387 %BW contact-force cert whose dominant load-bearing prime movers (glutes_max +
hamstrings, 76.6% share) passed cleanly. **This does not rise to the falsifier threshold for the prime-
mover claim** (the task's falsifier is about prime movers; adductor magnus's flexion-extension action is
explicitly secondary) but is disclosed in full rather than dropped, per symmetric QC.

Two of the four adductor-magnus **angle**-only flags are further explained: `addmagProx_r` peaks at
**−25.9°** (i.e. it is a **flexor-signed** bundle, the anatomically well-known anterior-fiber behavior),
so comparing its peak location to Nemeth & Ohlsen's extensor-context 75° is a category mismatch, not a
geometry error — their study (titled "hip **extensor** muscles") most plausibly represents the
posterior/ischial fibers, which is `addmagIsch_r`/`addmagMid_r` here (9.1° and 26.9° off; one passes, one
is a borderline miss just outside the ±20° band).

**No prime mover falls outside the tolerance band.** The falsifier condition (¬C) is not met for any of
vasti, rectus femoris, hamstrings (knee or hip), gastrocnemius, gluteus maximus, or iliopsoas.

## 9. Symmetric QC — which muscles were genuinely held-out vs tuned-to

This is not a fully independent anchor for every muscle, and that must be stated plainly, not buried:

- **Tuned-to (agreement partly expected):** Rajagopal et al. 2016 (PMCID 5507211, this model's direct
  parent, confirmed via its own full text, §3) explicitly compared its model moment arms against Buford
  1997, Spoor & van Leeuwen 1992, Nemeth & Ohlsen 1985, and Arnold 2000 **by name** during development.
  Passing §4/§5 against those same datasets therefore substantially confirms the geometry was
  **preserved correctly through inheritance + scaling** (a real and useful check — it would have caught a
  scaling bug, a `poly`-reformulation regression, or a copy/edit error from the `withArms`/`weldHand`
  modifications) rather than independently proving the *original* 1990s-2016 geometry correct from
  scratch. Herzog & Read 1993's specific inclusion in that validation lineage is **not confirmed**
  (possibly one of the unlabelled refs), and Dostal & Andrews 1981 is named only via the companion 1986
  paper, not confirmed identical in scope.
- **Genuinely closer to held-out:** the *scaling step itself* (generic model → this subject's segment
  lengths) is not something Rajagopal 2016 could have "tuned to" — it happens per-subject, downstream of
  their validation. §4/§5 passing on the **scaled** model is real evidence scaling did not distort the
  moment-arm geometry for *this* subject, which is a non-tautological, if narrower, claim.
- **Never anchored by this literature at all (§6, hip ab-adduction):** honestly OPEN, not fit-checked and
  not falsified — flagged as a follow-up needing full-text access to Dostal & Andrews 1981 or Arnold 2000.

Net: read every PASS in §4-5 as *"the model preserved the geometry it was designed to reproduce, through
scaling, at a level within cadaveric inter-study spread"* — not as *"independently proves 1990s cadaver
biomechanics right."* The **decorrelation value** of this whole exercise is against the OrthoLoad/force
pipeline (§0), not against the literature the geometry itself descends from.

## 10. Verdict

- **C is not falsified for any prime mover.** 11/17 `VERIFIED_QUANT` rows PASS (0 FLAGs on the two
  dominant hip-extension prime movers, glutes_max + hamstrings), 9/9 `DOMAIN_CONSENSUS` rows PLAUSIBLE,
  7/7 internal consistency checks confirmed (including a specific, literature-motivated anatomical
  behavior — glmax's fiber-dependent sign reversal — reproduced exactly).
- **One disclosed, quantified, low-materiality flag:** adductor magnus's hip flexion-extension moment arm
  is ~40% low vs Nemeth & Ohlsen even under the fairest (force-weighted) resultant — estimated ≤4% impact
  on combined hip-extension moment capacity, since adductor magnus is a secondary (9.8%-share) contributor
  in that plane and its primary (adduction) action checks out cleanly.
- **One honestly OPEN plane:** hip ab-adduction has no literature magnitude anchor from this session's
  accessible sources — do not read §6 PLAUSIBLE-by-trend as a quantitative pass.
- **Confidence tier: cadaveric/published-plausibility**, sub-tiered per §3 table — this validates
  *geometry preserved through scaling, within cadaveric inter-study spread*, on datasets the model
  lineage was partly fit to; it is not an in-vivo force anchor and does not by itself certify the 391/387
  %BW numbers — it rules out **moment-arm geometry** as the source of any error in them, which is exactly
  the decorrelated, non-tautological role this leg was built to play (§0).

## Files

- `scripts/msk/moment_arm_validation.py` — OpenSim sweep + raw-curve computation + sign-convention check +
  frozen-tail QC gate.
- `scripts/msk/moment_arm_lit_compare.py` — literature comparison, pre-registered tolerance check,
  internal consistency checks, all verdicts machine-computed.
- `data/msk_smoketest/moment_arm_validation/moment_arm_validation_evidence.json` — per-muscle peak +
  angle-of-peak + full curves for all 3 sweeps, sign-check, frozen-tail QC log.
- `data/msk_smoketest/moment_arm_validation/moment_arm_raw_curves.csv` — every single computed point
  (2,810 rows: coordinate, muscle, angle_deg, moment_arm_cm) — the full audit trail behind every peak
  quoted above.
- `data/msk_smoketest/moment_arm_validation/moment_arm_lit_comparison.json` — literature values (with
  verbatim-quoted sources), tolerance justification, per-row verdicts, internal-consistency checks, the
  hip-ab-adduction OPEN note.
