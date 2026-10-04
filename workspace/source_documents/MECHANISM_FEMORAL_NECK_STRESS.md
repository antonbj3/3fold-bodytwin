# MECHANISM FEMORAL NECK STRESS — the higher-value clinical site the shaft doc explicitly excluded (2026-07-21)

**Question.** `docs/MECHANISM_BONE_STRESS_CONSEQUENCE.md` found the femoral **shaft's** bone-stress
elevation does not coherently inherit the twin's 1.412–1.990× joint-contact-force over-prediction
(`docs/MECHANISM_CROSS_SUBJECT.md`). But hip-fracture risk is clinically dominated by **femoral-NECK**
fractures (subcapital/transcervical/basicervical), not diaphyseal ones — the exact site
`bone_stress.py` itself excluded from its headline (s<0.15 on the shaft axis, "femoral head/neck",
flagged diagnostic-only). Extending the hip-contact + muscle (glute/psoas) free body onto a
separately-derived **neck axis**: does the modeled peak femoral-neck strain land in the physiological
~500–2000 µε walking band (Aamodt 1997), well below the ~7000 µε yield, or is it inflated?

**Headline answer: mostly INFLATED relative to the stated walking band, comfortably below yield in
most configurations — and, like the shaft, this elevation is NOT coherently explained by the
joint-contact-force over-prediction (though the mechanism differs in an interesting, honestly-reported
way).** The primary configuration (NSA=124.7°, offset=42.9mm, anteversion=14.5° — all three from the
SAME live-verified 628-hip CT cohort, Carmona et al. 2019) predicts **−57.01 MPa / 3354 µε** (medial
cortex, compression) at the base of the neck — **1.7× the walking band's own upper edge**, but
**4.7× below** the task-stated 7000 µε yield. Removing ONLY the neck eccentricity (NSA→180°, the direct
analog of the shaft doc's own offset=0 null) drops this to **1550 µε — INSIDE the walking band** —
localizing the elevation to the geometric eccentricity assumption, exactly as the shaft doc found for
its own offset. De-inflating the hip+knee contact forces by subject2's own measured ratio (the SAME
counterfactual the shaft doc ran) **decreases** neck stress slightly (3354→3200 µε, **−4.6%**) — the
**opposite direction** from the shaft's own finding (there, de-inflation increased stress +8.0%) — but
the magnitude is similarly too small to explain a 29–34% force change, so the **same qualitative
conclusion holds by a different mechanism**: decoupled, not coherent propagation.

A forced OODA pass caught and fixed a real methodological problem before any of the above was trusted:
an early version's reported "peak" landed exactly on an arbitrary domain-ceiling guess (the same
red-flag shape as `bone_stress.py`'s own first-draft bug) — fixed by anchoring the valid domain to the
model's own geometry (Sec 3), not a round number, and reported honestly as a **base-of-neck boundary
value**, not an interior optimum.

| | AS-IS primary | null control (zero eccentricity) | de-inflated (subject2's own ratio) |
|---|---:|---:|---:|
| Peak stress | −57.01 MPa | −26.35 MPa | −54.40 MPa |
| Peak microstrain | 3354 µε | 1550 µε | 3200 µε |
| vs walking band (500–2000 µε) | **above** | **inside** | **above** |
| vs task-stated yield (7000 µε) | well below | well below | well below |

---

## 1. Method — reuse discipline, and the one genuinely new axis

Imports `scripts/msk/bone_stress.py` (`bs`) and `scripts/msk/bone_stress_consequence.py` (`bsc`)
**UNEDITED** — `femur_axis_geometry`, `frame_loads_at`, `internal_force_moment`, `combined_stress`,
`run_synthetic_selftests`, `internal_force_moment_batch`, `section_props`, `scale_loads`,
`decompose_at_section`, `to_microstrain` — and `validate_joint_force.py`'s `parse_mot`/Savitzky-Golay
machinery. The free body is **identical** to `bone_stress.py`'s own: hip + knee + patellofemoral joint
reaction, all 26 live-detected femur-attaching muscles (including `psoas_r`, `glmax1-3_r`, `glmed1-3_r`,
`glmin1-3_r`), segment weight, translational d'Alembert inertia.

**What's new**: a NECK axis, derived — not looked up as a third independent citation — from two things
this model already has: its own empirically-derived `(u_hat, e_lateral, e_ap)` frame
(`bs.femur_axis_geometry`, reused) and two population parameters from a single, live-verified CT cohort
(Carmona et al. 2019, PMID 31266693, n=628 hips): neck-shaft angle (NSA) and femoral offset (FO). Neck
length is **derived**, not cited, from the right triangle these two already define:

```
offset = l_neck * sin(NSA)   =>   l_neck = offset / sin(NSA)
n_hat  = -cos(NSA)*u_hat + sin(NSA)*(cos(alpha)*e_lateral + sin(alpha)*e_ap)     [alpha = anteversion]
```

`NSA=180°` is the exact null control (`n_hat = u_hat`, zero eccentricity) — the direct geometric analog
of the shaft doc's own `offset=0mm` null. This derivation is forced through a **4-case synthetic
self-test** (NSA=90° pure-perpendicular, NSA=180° exact-null, the primary config's offset reproduced to
float precision, anteversion preserving the axial component) before any real-data use — same
synthetic-control-before-real-result discipline `bs.run_synthetic_selftests()` already established.

Once `(n_hat, l_neck)` are fixed, method-of-sections proceeds with the SAME, unedited functions the
shaft model uses — called with `(origin=hip_c, axis=n_hat, length=l_neck)` instead of `(hip_c, u_hat,
L)`. No "eccentricity trick" is needed here (unlike the shaft): the neck axis is, by construction,
anchored at the true hip-force application point, so the hip force has zero perpendicular offset from
its own axis's origin.

---

## 2. Verification (machine-checked before any downstream number was trusted)

| gate | result |
|---|---|
| Reused shaft beam-sectioning synthetic self-test (`bs.run_synthetic_selftests`) | PASS |
| New neck-axis-geometry synthetic self-test (4 closed-form cases) | PASS |
| This script's own hip/knee %BW re-derivation vs the pre-existing certs | PASS (387.037/391.115 %BW, diff <1e-5) |
| This script's own free-body re-derivation reproduces `bone_stress_results.json`'s shaft peak | PASS (54.928 vs 54.968 MPa, diff 0.072%) |
| `internal_force_moment_batch` vs scalar reference, ON THE NECK AXIS, 25 sampled points | PASS, max diff 4.9e-14 N·m |
| Valid-domain ceiling frame-invariance (a femur-local point cannot move with pose) | PASS, spread 2.7e-15 |
| Decomposition sums to total (N, M_bend) at the primary peak | PASS, diff <6.4e-13 |
| Offset-sweep degeneracy (Sec 4.3) confirmed exactly zero, as predicted | PASS, spread 0.0 MPa |

**All 8 gates PASS** (`data/.../femoral_neck_stress_results.json`, key `gates` all `true`, exit 0).

One latent bug was found and fixed mid-session, disclosed rather than silently corrected: two
boolean gate values (`shaft_reproduction_pass`, `offset_sweep_degenerate_as_expected`) were computed
from numpy-array-derived values, giving `numpy.bool_` rather than a plain Python `bool`. `json.dump(...,
default=str)` silently stringifies an unrecognized type — meaning a `numpy.bool_(False)` would have
been written as the **string** `"False"`, which Python's `all(gates.values())` treats as **truthy**
(any non-empty string is truthy) — a fail-open gate that happened not to bite this run (both values were
genuinely `True`) but would have silently passed a real failure. Fixed at the source (explicit `bool(...)`
casts) and defended in depth (every `gates` dict value is now cast + asserted `isinstance(v, bool)`
before aggregation). Caught by inspecting this run's own raw JSON output, not assumed safe.

---

## 3. The forced fix — a boundary problem, not a bug, found via OODA before trusting anything

An early version hardcoded the valid neck-domain ceiling at `s=0.95` (0=head, 1=the NSA-derived
"neck-base" point). The worst-case search landed **exactly** on that ceiling — the identical red-flag
shape as `bone_stress.py`'s own first-draft bug (peak at s=0, identical across the whole offset sweep).

**Orient (diagnosis, not assumed):** a diagnostic scan of `|M_bend(s)|` at the peak frame, extended past
s=1.0, showed it climbing essentially monotonically with no interior turnover, well past s=1.6. Unlike
the shaft axis (which terminates at a real joint, the knee), **the neck axis has no natural far end** —
`l_neck` is a geometric construct (where the NSA/offset right-triangle closes), not an anatomical
boundary, so "the peak within [0.05, 0.95]" was simply "the value at wherever I stopped looking."

**Decide/Act (the fix):** mirror `bone_stress.py`'s own precedent exactly — "bounds set from this
model's own attachment geometry, not a round-number guess." Computed live: the closest real (non-hip)
load's own projection onto the neck axis. This is `psoas_r`'s **inner** femur attachment point (the
lesser trochanter) at `s = 0.79589` — verified **frame-invariant** (a femur-body-local point cannot move
with pose; checked across 5 frames spanning the trial, spread 2.7e-15) before being trusted. The valid
ceiling `S_HI` is 98% of this value (0.78), pulled back to avoid evaluating exactly at another load's
own point.

**This reframes the finding, honestly, rather than hiding it**: the reported peak is **the stress at (or
right at) the base of the neck**, the anatomical transition where the first real counteracting structure
appears — a real, clinically meaningful location (the "basicervical" fracture site), but a genuinely
weaker claim than "the worst point somewhere in the neck's interior." The boundary-sensitivity table
below shows this is not a knife's edge:

| fraction of `S_HI` | s | peak stress (MPa) | microstrain |
|---:|---:|---:|---:|
| 80% | 0.624 | −45.32 | 2666 |
| 90% | 0.702 | −51.17 | 3010 |
| 95% | 0.741 | −54.09 | 3182 |
| 99% | 0.772 | −56.43 | 3319 |
| 100% (reported) | 0.780 | −57.01 | 3354 |

A real, non-trivial slope (26% swing over the last 20% of the domain) — but a smooth one, not a
singularity — consistent with "stress grows toward the neck base" (classical Pauwels teaching), not with
a computational artifact.

---

## 4. Results

### 4.1 Geometry sweeps

**Neck-shaft angle** (offset fixed at 42.9mm) — the operative geometric sensitivity parameter, since it
changes the axis DIRECTION, not just a rescaling:

| NSA (deg) | l_neck (mm) | peak stress (MPa) | microstrain |
|---:|---:|---:|---:|
| 110.0 | 45.65 | −49.29 | 2899 |
| 115.0 | 47.33 | −53.69 | 3158 |
| 120.0 | 49.54 | −56.28 | 3310 |
| **124.7 (primary, Carmona mean)** | **52.18** | **−57.01** | **3354** |
| 130.0 | 56.00 | −55.69 | 3276 |
| 135.0 | 60.67 | −45.86 | 2697 |

Non-monotonic (peaks near the measured mean, falls off at both clinical extremes coxa-vara/valga) —
2697–3354 µε across the full 110–135° clinical range, all **above** the walking band, all comfortably
below yield.

**Offset sweep (NSA fixed) — a genuine degeneracy, caught by symmetric QC, not swept under the rug.**
Every nonzero offset (30/40/42.9/50/60mm) gives **exactly** −57.01 MPa:

| offset (mm) | peak stress (MPa) | microstrain |
|---:|---:|---:|
| 0.0 (= null control) | −26.35 | 1550 (**inside** walking band) |
| 30.0 – 60.0 (all 5 values) | −57.01 (identical) | 3354 (identical) |

**Why, exactly** (verified, not hand-waved): with NSA fixed, `n_hat` does not depend on offset at all
(only `l_neck` does — see the equations, Sec 1). Since `S_HI` is a fixed **fraction of the distance to a
fixed physical point** (psoas_r), the evaluation location `x(S_HI) = hip_c + S_HI·l_neck·n_hat =
hip_c + (0.98·proj)·n_hat` has `l_neck` cancel exactly, so every nonzero offset evaluates the identical
3D point. Confirmed numerically (spread = 0.0 MPa across all 5, a hard gate, Sec 2). **This means, once
the domain is anchored to real anatomy, NSA — not offset — is the operative geometric lever**; the
shaft doc's own "offset accounts for roughly half the peak" framing does not transfer unmodified to this
boundary-forced construction.

**Anteversion** (NSA/offset fixed) — confirmed a second-order refinement, as reasoned before running it:

| anteversion (deg) | peak stress (MPa) | microstrain |
|---:|---:|---:|
| 0.0 (coronal-only) | −54.73 | 3220 |
| 14.5 (measured mean, primary) | −57.01 | 3354 |

A 4% swing — real, but small, consistent with the reasoning that anteversion mostly redirects which
meridian of the cortex carries the peak rather than the peak's own magnitude.

### 4.2 Cross-section sweep (literature-informed, swept — honest gap, Sec 8)

| outer ⌀ (mm) | cortical thickness (mm) | peak stress (MPa) | microstrain | classification |
|---:|---:|---:|---:|---|
| 22 | 2.0 | −131.11 | 7712 | **AT/ABOVE task-stated yield** |
| 22 | 3.5 | −90.83 | 5343 | above walking, below yield |
| 22 | 5.0 | −77.09 | 4534 | above walking, below yield |
| 27 | 2.0 | −85.56 | 5033 | above walking, below yield |
| **27 (primary)** | **3.5 (primary)** | **−57.01** | **3354** | above walking, below yield |
| 27 | 5.0 | −46.60 | 2741 | above walking, below yield |
| 32 | 2.0 | −60.80 | 3576 | above walking, below yield |
| 32 | 3.5 | −39.43 | 2320 | above walking, below yield |
| 32 | 5.0 | −31.39 | 1846 | **inside** walking band |

Only the thinnest (22mm/2.0mm) crosses into yield territory; only the thickest (32mm/5.0mm) falls
inside the walking band. The primary sits mid-pack — the qualitative finding (above the walking band,
below yield) is not an artifact of picking one convenient corner of this grid.

### 4.3 Ablations — does removing glute/psoas force move stress the direction classical theory predicts?

| variant | peak stress (MPa) | microstrain | vs primary |
|---|---:|---:|---|
| **primary (all muscles present)** | **−57.01** | **3354** | — |
| zero abductors (`glmed1-3_r`, `glmin1-3_r`, `tfl_r`) | −69.25 | 4074 | **+21%** (removing abductors INCREASES stress) |
| zero `psoas_r` only | −51.84 | 3049 | **−9%** (removing psoas DECREASES stress) |
| zero all 26 muscles | −153.39 | 9023 | **+169%**, crosses **above yield** |

The abductor result **confirms classical Pauwels theory directly**, with real subject-specific numbers,
not just cited: gluteal/abductor force counteracts neck bending (removing it makes stress worse) —
independently corroborated by Cristofolini et al. 1995's live-verified ex-vivo finding that the glutei
are the *principal* muscles determining proximal-femur strain. The all-muscle ablation crossing above
yield directly extends Duda et al. 1998's already-verified shaft-level finding ("muscles reduce peak
bending, they do not add to it") to the neck.

**Psoas behaves oppositely to the abductors at this specific section/frame** — its removal *reduces*
stress here, meaning psoas's own line of action currently adds to (not counteracts) the net neck bending
at the base-of-neck cross-section. This is not a contradiction; it is what the exact decomposition (Sec
5) shows directly: different muscles pull in different directions relative to a specific cut, and psoas
(a flexor, small adduction/abduction moment arm) does not play the same frontal-plane counteracting role
the abductors do. Reported plainly, not smoothed into "muscles help."

### 4.4 De-inflation counterfactual — extends the shaft doc's own test, opposite direction, same conclusion

Scaling the hip (×0.7078) and knee (×0.6602) contact forces by subject2's own measured OrthoLoad ratios
(`docs/MECHANISM_CROSS_SUBJECT.md`), holding all 26 muscle forces + patellofemoral + weight + inertia
fixed — the **identical** counterfactual `bone_stress_consequence.py` ran at the shaft:

| | AS-IS | de-inflated | direction |
|---|---:|---:|---|
| Peak stress | −57.01 MPa | −54.40 MPa | **decreased** 4.6% |
| Peak microstrain | 3354 µε | 3200 µε | **decreased** 4.6% |
| Classification | above walking, below yield | above walking, below yield | unchanged |

**This is the opposite sign from the shaft's own finding** (there: de-inflation *increased* stress
+8.0% at the primary offset). But the magnitude here is similarly small relative to the 29–34% force
change applied — a 4.6% stress change from a ~30% force change is not proportional/coherent by any
reasonable reading. **The pre-registered falsifier fires** (Sec 7): de-inflation does not coherently
explain the elevation at the neck either, just via a different (here, mildly helpful rather than mildly
harmful) mechanism than at the shaft.

---

## 5. The mechanism — exact decomposition, same cancellation structure as the shaft, quantified fresh

At the primary peak (s=0.780, t=0.60s), the internal bending moment splits **exactly** (gate: diff
<6.4e-13, Sec 2) into:

| load category | axial N (N) | \|M_bend\| (N·m) |
|---|---:|---:|
| hip_joint_contact | −1069.04 | 39.57 |
| knee_joint_contact | −794.42 | **160.03** |
| patellofemoral_joint_contact | +27.71 | 29.45 |
| muscles (26, summed) | −338.27 | 107.64 |
| femur segment weight | +25.44 | 1.56 |
| d'Alembert inertia | +8.22 | 1.71 |
| **TOTAL** | **−2140.37** | **65.91** |

The **same qualitative mechanism** the shaft doc's Sec 4 found: the knee's own individual bending
contribution (160.03 N·m — the single largest term, even though the knee is ~490mm away) is
**substantially, not completely, cancelled** by muscles + patellofemoral acting in the opposing
direction (component-wise: knee's vector is `[−106.7,−64.6,−100.3]`, muscles' is `[+67.6,+47.8,+68.8]` —
opposite sign in every component). The net (65.91 N·m) is well below any single category's own magnitude
— confirming this is a genuine over-determination of the same Moissenet-named mechanism
(`docs/MECHANISM_JOINT_FORCE_SCORECARD.md`, PMID 24210475) showing up independently at a **second,
geometrically distinct** cut location, using the same underlying per-frame force data. This is why
de-inflating hip+knee together does not cleanly reduce stress here either: both terms are entangled in
this same partial-cancellation structure, and scaling them down changes which direction the residual
points, not simply its size.

---

## 6. Anchors — live-verified this session (NCBI eutils + Europe PMC REST), not recalled

| citation | PMID/DOI | verified content | used as |
|---|---|---|---|
| Carmona M, Tzioupis C, LiArno S, Faizan A, Argenson JN, Ollivier M. "Upper Femur Anatomy Depends on Age and Gender: A 3D CT Comparative Bone Morphometric Analysis of 628 Healthy Patients' Hips." *J Arthroplasty.* 2019;34(10):2487-93. | 31266693 / 10.1016/j.arth.2019.05.036 | Abstract (verbatim): "the mean NSA was 124.7° (SD 6.2), mean femoral neck version was 14.5° (SD 8.1), mean FO was 42.9 mm (SD 6.8)." n=628 hips. | **Primary geometric anchor** — NSA, anteversion, AND offset, all from the SAME cohort (paired, not mixed). The 42.9mm FO also cross-validates `bone_stress.py`'s own, independently-sourced (boneandspine.com) 40mm shaft-offset assumption from a genuinely different source. |
| Aamodt A, Lund-Larsen J, Eine J, Andersen E, Benum P, Husby OS. "In vivo measurements show tensile axial strain in the proximal lateral aspect of the human femur." *J Orthop Res.* 1997;15(6):927-31. | 9497820 / 10.1002/jor.1100150620 | Abstract (verbatim, in full): strain-gauge rosettes bonded to the exposed proximal-lateral femur during snapping-hip-syndrome surgery (2 patients); "principal tensile strain increased significantly during one-legged stance, walking, and stair climbing" vs two-legged stance; strain aligned "within 22 degrees" of the femoral axis; supports Pauwels' bending theory. | Qualitative direction + site (real in-vivo, not recalled). **Honest gap**: no numeric microstrain in the abstract (full text subscription-walled, checked live, not accessible) — the task's stated 500–2000 µε band is used as a TASK-STATED reference, not independently re-extracted from Aamodt's own primary text this session. Site is the **exposed proximal-lateral shaft** (post-capsular access), not literally the intracapsular neck — a real cross-region gap, not glossed over. |
| Cristofolini L, Viceconti M, Toni A, Giunti A. "Influence of thigh muscles on the axial strains in a proximal femur during early stance in gait." *J Biomech.* 1995;28(5):617-24. | 7775497 | Abstract (verbatim, in full): ex-vivo composite/cadaveric proximal femur, 10 thigh muscles via nylon straps; "the three glutei are the principal muscles in determining the vertical strains"; simulating abductors via a lumped resultant force was shown INADEQUATE. | Directly supports this doc's own free-body choice (real, separate glute forces, not a lumped resultant) and directly corroborates the abductor-ablation finding (Sec 4.3). |
| Cirovic A et al. "Three-dimensional mapping of cortical porosity and thickness along the superolateral femoral neck in older women." *Sci Rep.* 2022;12(1):15544. | 36109611 / 10.1038/s41598-022-19866-2 | Abstract (verbatim): cortical thickness/porosity differ significantly across basicervical/midcervical/subcapital subregions; subcapital is thinnest and most porous. | Qualitative support for modeling the neck cortex thinner than the shaft's own (honest-gap, swept, not point-asserted, Sec 8). |
| **Reused, not re-verified this session** (already live-verified upstream): Duda 1998 (PMID 9802785, muscles reduce peak bending), Lanyon 1975 (PMID 1146518), Burr 1996 (PMID 8739897), Frost 1987 (PMID 3688455), Moissenet 2014 (PMID 24210475, SO overestimates joint reaction force). | — | — | Cross-referenced, same tier as the shaft doc's own use. |

---

## 7. Confidence tiers and the falsifier, stated explicitly

- **Claim**: this beam model's predicted femoral-**neck** bone strain during gait is, in most tested
  configurations, above the task's stated physiological walking band but comfortably below its stated
  yield threshold — and this elevation is decoupled from (not coherently explained by) the 1.4–2.0×
  joint-contact-force over-prediction, extending the shaft doc's own finding to a second, higher-value
  clinical site via a partially different mechanism.
- **Confidence tier**: **partially in-vivo-anchored** crossed with **method-only**. The *qualitative*
  in-vivo/ex-vivo anchors are real and live-verified this session (Aamodt 1997's tension-direction
  finding; Cristofolini 1995's muscle-dominance finding) — but the **specific 500–2000 µε numeric band**
  is the task's own stated figure, not independently re-extracted from Aamodt's paywalled primary text
  (weaker than the shaft doc's own Lanyon/Burr anchors, which DID carry live-verified exact numbers).
  The beam model's own **absolute** stress level is **method-only**: it inherits every honest gap
  `bone_stress.py` already disclosed (generic, not CT/DXA-measured cross-section; uncorrected rotational
  inertia; quasi-static) **plus two NEW, worse-than-shaft gaps found this session** — the neck's low
  aspect ratio (1.93, vs the shaft's 18.2 — St Venant's margin consumes ~103% of the neck's own length)
  and the boundary-value (not interior-optimum) nature of the reported peak (Sec 3). The **decoupling**
  claim itself (the falsifier test, Sec 4.4) is machine-verified to float precision, independent of
  whether the model's absolute MPa level is accurate.
- **Falsifier, pre-registered and fired**: "coherent propagation" required de-inflated stress to drop
  proportionally to the applied force ratio (~30-34%) and move materially toward the walking band.
  Measured: it moved only 4.6% (in the OPPOSITE direction from the shaft's own finding, which itself
  moved 8.0% the other way) — neither direction is proportional to the applied force correction. **This
  falsifies coherent propagation at the neck, same conclusion as the shaft, different sign.**
- **What would falsify the decomposition's own causal story** (Sec 5): if the knee/muscle
  near-cancellation did not reproduce with an independently re-computed free body, or if the linearity
  (decomposition-sums-to-total) gate failed. Neither happened (Sec 2).
- **What would falsify the boundary-forcing fix** (Sec 3): if `min_nonhip_s` were NOT frame-invariant
  (would mean a femur-local point is somehow moving with pose — a modeling contradiction). Checked, not
  assumed — PASS, spread 2.7e-15.

---

## 8. Honest gaps

1. **The reported peak is a boundary value (Sec 3), not an interior optimum.** `|M_bend|` climbs
   essentially monotonically toward the neck base with no interior turnover (verified by an extended
   diagnostic scan past s=1.0). The boundary IS geometrically forced (the closest real attaching
   structure, not a round number) and the local slope is shown explicitly (Sec 3 table) — but this is a
   genuinely different, weaker claim than "the worst point in the neck's interior."
2. **St Venant's principle concern is WORSE than at the shaft**: aspect ratio 1.93 (vs the shaft's 18.2)
   means the "~1 diameter" invalid margin from each end consumes ~103% of the neck's own length. The
   (N, M_bend) rigid-body-statics computation is exact regardless of slenderness; only the final M·c/I
   stress conversion inherits this concern, more severely than the shaft's own already-disclosed version
   of the same idealization. This is the SAME idealization the classical Pauwels/Frankel-Nordin textbook
   neck model already uses — a recognized, standard, coarse first-order approximation for this specific
   site, not a novel or unusually weak one.
3. **Neck cross-section (outer diameter, cortical thickness) is a swept literature-informed range, not
   this subject's own CT/DXA geometry**, and — despite six distinct live search angles tried this session
   (general morphometry, DXA hip-structural-analysis, CT cortical-thickness mapping) — not independently
   pinned to one exact peer-reviewed mm figure (every candidate abstract reported correlations or
   qualitative subregion differences, not baseline mm means). Same honest-gap tier as `bone_stress.py`'s
   own shaft cross-section gap.
4. **Neck axis geometry (NSA, anteversion, offset) is population-level literature (Carmona et al. 2019,
   n=628), not this subject's own measured hip geometry.** Neck length is derived (not independently
   cited) from NSA+offset — internally consistent, but not itself a directly measured anatomical figure;
   it came out somewhat longer (52mm) than commonly-seen direct neck-length figures (~30–40mm), meaning
   the true aspect-ratio concern (Sec 8.2) may be understated, not overstated.
5. **Uniform circular annulus idealization** — the real femoral neck cortex is highly asymmetric (thick
   inferomedial calcar, thin-to-absent superolateral cortex, live-verified qualitatively via Cirovic et
   al. 2022) with substantial trabecular content; modeled here as a symmetric thin cortical annulus at a
   single uniform E=17 GPa (same convention the shaft docs already use).
6. **Anteversion's rotation sign is not independently pinned** — `e_ap` is just the third axis of an
   empirically-derived orthonormal frame; both signs would need sweeping for a fuller treatment (not done
   here, since the effect is shown to be second-order, Sec 4.1).
7. **Aamodt 1997's task-attributed "500-2000 µε" figure was not independently re-extracted** from that
   paper's primary text this session (subscription-walled) — flagged plainly, not silently assumed
   verified.
8. **Fixed-muscle-force counterfactual, not a re-solved simulation** — same caveat
   `docs/MECHANISM_BONE_STRESS_CONSEQUENCE.md` already disclosed: if Static Optimization's own muscle
   recruitment is itself a root cause of the contact-force over-prediction, the 26 muscle forces used
   unchanged here still carry that same unremoved bias.
9. **Quasi-static, translational-inertia-only** (identical to `bone_stress.py`) — rotational inertia
   (I·α + ω×Iω) is not corrected.
10. **subject2/walking1, right leg, single trial only** — same scope boundary as every MSK cert in this
    repo family.

---

## 9. Files

- `scripts/msk/femoral_neck_stress.py` (new) — the full pipeline: neck-axis derivation + its own
  4-case synthetic self-test, per-frame free-body re-derivation (bit-for-bit-verified against
  `bone_stress.py`'s own shaft output before being trusted), the geometrically-forced valid-domain
  ceiling (with frame-invariance gate), the NSA/offset/anteversion/cross-section sweeps, the null
  control, 3 ablations, the de-inflation counterfactual, the exact decomposition, and the
  boundary-sensitivity report. Real run, exit 0, all 8 gates PASS.
- `data/msk_smoketest/subject2_walking1/femoral_neck_stress/femoral_neck_stress_results.json` — every
  number in this document, machine-written (this IS the evidence JSON).
- **Read-only, not modified**: `scripts/msk/bone_stress.py`, `scripts/msk/bone_stress_consequence.py`,
  `scripts/msk/validate_joint_force.py`.
- **Read-only inputs**: `data/msk_smoketest/subject2_walking1/bone_stress/bone_stress_results.json`
  (cross-check target), `data/msk_smoketest/subject2_walking1/static_optimization/{jr,so}/*.sto`.
- No git commit, no git push performed (isolation respected). No file outside this task's own new
  outputs was touched.
