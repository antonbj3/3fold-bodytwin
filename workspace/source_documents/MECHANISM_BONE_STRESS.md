# MECHANISM BONE STRESS — cortical femur stress below joint contact force (2026-07-21)

Adds the layer BELOW the twin's already-validated joint contact forces (hip 387.04 %BW, knee 391.11
%BW — both re-verified live below from the exact same `.sto` the prior certs used) and Static
Optimization muscle tensions: what CORTICAL BONE STRESS do those forces induce in the femur? Every
number below is machine-measured this session (`scripts/msk/bone_stress.py`, exit 0), not recalled.
Isolation respected: `.venv-msk` only, all outputs new files under `data/msk_smoketest/subject2_walking1/
bone_stress/`, no git commit/push.

## Headline result

| | value |
|---|---:|
| **Peak combined cortical stress (primary: 40mm offset, ⌀27mm/6.5mm-thick section)** | **−54.97 MPa (compression)** |
| Location | s=0.156 (78mm distal to the hip center) — **subtrochanteric/proximal shaft** |
| Side | **medial cortex** |
| Instant | t=0.56s (36% of this trial) |
| Axial component | −3.13 MPa (compression) |
| Bending amplitude component | 51.84 MPa |
| Task-stated 40-60 MPa walking band | **PASS** |
| Sane-order-of-magnitude gate (10-150 MPa, pre-registered) | PASS |

Both the **location** (subtrochanteric) and the **side** (medial) match the task's own stated
expectation — an independent corroboration from actual subject-specific forces, not just repeating
the expectation. **But read the honesty section (§6) before trusting the single number above** — the
quasi-static residual leaves a real, quantified uncertainty band of roughly 19-92 MPa on the bending
component alone, which this document does not hide.

## 1. Method — Euler-Bernoulli beam, method of sections, femur free body

**Boundary loads** (the task's own ask): the corrected hip contact force (387.04 %BW) + corrected knee
contact force (391.11 %BW) + every muscle that actually **attaches to the femur** — found LIVE from
the model's own current path geometry every frame (`getCurrentPath`), never a hand-picked "major
muscles" list. Measured, n=26 distinct muscles with ≥1 femur attachment across the trial:
`addbrev_r, addlong_r, addmagDist_r, addmagIsch_r, addmagMid_r, addmagProx_r, bfsh_r, gaslat_r,
gasmed_r, glmax1_r, glmax2_r, glmax3_r, glmed1_r, glmed2_r, glmed3_r, glmin1_r, glmin2_r, glmin3_r,
iliacus_r, piri_r, psoas_r, sart_r, tfl_r, vasint_r, vaslat_r, vasmed_r` — includes the abductors
(glmed/glmin, greater trochanter), iliopsoas (lesser trochanter), adductors (linea aspera), vasti +
biceps-femoris-short-head (shaft to patella/fibula), and gastrocnemius (condyles). This is a
**superset** of "abductor forces" as asked — using the full live-detected set is strictly more
complete and less cherry-picked than hand-selecting "the major ones," and the sensitivity of the
result to the abductor/hip mechanism specifically is quantified directly (§3, the offset=0 null).

**Beam axis**: the straight line from the hip joint center to the knee joint center, in the femur's
OWN body frame (rigid, pose-independent) — **length 491.08mm, this model's own scaled geometry**, not
a literature number.

**Femoral-neck offset** (the dominant lever, per the task's own framing): a rigid-body joint model has
no medullary-canal geometry, so the ~30-60mm literature offset between the true head center and the
shaft's own centroidal axis (boneandspine.com clinical reference; consistent with THA-literature
"native offset" targets ~40-44mm) is applied as an eccentricity on the hip force's application point,
in a direction **empirically derived from this model's own abductor attachment geometry** (not
assumed) — this model's own implied trochanteric offset from the same axis is **66.48mm**, the same
order of magnitude as the literature head-offset band, a soft consistency cross-check.

**Method of sections**: cut the beam at parametric position s∈[0,1] (0=hip, 1=knee); the internal
axial force N (tension-positive) and bending-moment vector at each cut come from the AVERAGE of the
distal-sum and proximal-sum conventions (identical in a perfectly balanced free body — see §4 for why
averaging was necessary, not optional). For a circular cortical annulus, peak bending stress magnitude
at any section is exactly `c·|M_bend|/I` regardless of bending direction (closed form, no eyeballed
figure) — this also gives which side (medial/lateral) carries the tensile vs compressive peak, via the
model's own empirical lateral-direction vector.

**Valid-domain restriction** (forced by an OODA pass, §5): the headline result is restricted to
s∈[0.15, 0.85] — the true diaphysis. The excluded region (femoral head/neck, s<0.15; condylar flare,
s>0.85) is where (a) St. Venant's principle says beam theory itself is invalid within about one
section-diameter of a concentrated load, and (b) a hollow circular-tube cross-section does not
anatomically apply. The unrestricted whole-domain scan is still computed and reported (peak −57.77 MPa
at s=0.044, "femoral head/neck") — close in magnitude, clearly labeled DIAGNOSTIC, not the headline.

## 2. Cortical cross-section — literature, swept (explicit honest gap)

No subject-specific CT/DXA femoral cross-section geometry exists for this subject. A live search this
session for a peer-reviewed exact midshaft mm²/mm⁴ figure hit paywalls (Springer, a 62-femur
histomorphology study PMID 28548343 covering ages 51-95 — wrong population anyway) without a citable
number. Rather than assert a single fabricated-precision figure, cross-section is modeled as a
circular annulus and **swept**, not point-asserted:

| outer ⌀ (mm) | cortical thickness (mm) | peak \|σ\| (MPa) |
|---:|---:|---:|
| 24 | 5.0 / 6.5 / 8.0 | −81.83 / −75.30 / −72.59 |
| 27 (primary) | 5.0 / **6.5** / 8.0 | −60.85 / **−54.97** / −52.20 |
| 30 | 5.0 / 6.5 / 8.0 | −47.02 / −41.82 / −39.17 |

Range across this plausible adult geometry grid: **−39 to −82 MPa** — same order of magnitude
throughout, moving the headline number but not the qualitative finding.

## 3. Femoral-offset sensitivity — the task's own flagged mechanism, quantified

| offset (mm) | peak \|σ\| (MPa) | location |
|---:|---:|---|
| 0 (null control) | −26.75 | s=0.248, diaphysis, t=0.09s |
| 30 | −47.19 | s=0.156, subtrochanteric, t=0.56s |
| **40 (primary)** | **−54.97** | s=0.156, subtrochanteric, t=0.56s |
| 50 | −62.82 | s=0.156, subtrochanteric, t=0.56s |
| 60 | −70.72 | s=0.156, subtrochanteric, t=0.56s |

The eccentricity accounts for **28.22 MPa of the 54.97 MPa primary peak — just over half**. This is
the direct, measured version of the task's own qualitative warning ("the femoral-neck moment arm
dominates — flag it"): confirmed, and quantified rather than merely asserted. It also shifts WHERE the
peak occurs (diaphysis at offset=0 → subtrochanteric once the eccentricity is present), matching the
clinical/biomechanical expectation that the neck offset is what makes the subtrochanteric region the
proximal femur's highest-bending-moment location.

## 4. A forced fix, not a first-draft result — two real bugs found via OODA, not narrated away

Two earlier drafts of this analysis returned nonsense, and each was traced to a concrete, fixed cause
rather than accepted or hand-waved:

1. **First draft**: peak stress landed at s=0.000 (the hip joint center itself), IDENTICAL across every
   offset from 0-60mm. Diagnosis: summing only the "distal" side of each cut (proximal part acts on
   distal part) means the hip force — always on the proximal side for any cut beyond it — can
   structurally never enter the calculation for s>0, regardless of how good or bad the underlying
   physics is. **Fix**: compute the internal moment via BOTH the distal-sum and proximal-sum
   conventions (identical in a perfectly balanced system — proved via total-moment-about-the-cut=0)
   and average them; the classical Pauwels/Frankel-Nordin proximal-sum convention is what actually lets
   the hip eccentricity propagate. Verified against three closed-form synthetic beam cases (pure axial,
   a balanced simply-supported beam with a known triangular moment diagram, and a deliberately
   IMBALANCED case whose discrepancy diagnostic was checked against its own closed form) before
   trusting it on real data.
2. **Second draft** (after fix #1): whole-femur free-body residual was 565N (73.75 %BW) — implausibly
   large. Diagnosis: the model has a `patellofemoral_r` joint (femur↔patella) whose reaction force
   (the patella pressing into the trochlear groove — a real contact force the vasti's tension routes
   through) was missing from the free body. **Fix**: added it (mirrors the knee's Newton's-3rd-law +
   frame-transform pattern). Residual dropped to 45.0N (5.87 %BW) — confirms the diagnosis rather than
   just asserting it.

A third, smaller fix: a **translational d'Alembert correction** (`−m·a_com`, from the femur segment's
own live COM-acceleration trajectory, same Savitzky-Golay differentiator `validate_joint_force.py`
already uses) drops the residual further to 23.1N (3.02 %BW) in the trial's interior. The remaining
42.6 N·m moment residual is the **uncorrected rotational-inertia term** (I·α + ω×Iω) — an explicit,
quantified, NOT-forced-to-zero honest gap (§6).

## 5. Anchors — live-verified PMID/DOI this session, not recalled

- **Duda GN, Heller M, Albinger J, Schulz O, Schneider E, Claes L. "Influence of muscle forces on
  femoral strain distribution." J Biomech. 1998 Sep;31(9):841-6.** PMID
  [9802785](https://pubmed.ncbi.nlm.nih.gov/9802785/), DOI
  [10.1016/s0021-9290(98)00080-3](https://doi.org/10.1016/s0021-9290(98)00080-3). Finite-element
  femur, real muscle+joint forces over a gait cycle. Finding used here: peak surface strain **<2000
  microstrain with all thigh muscles** (~45% gait cycle) vs **~3000 microstrain under simplified
  loading** — muscles REDUCE peak bending, they do not add to it. Converted to a stress band (E=17
  GPa, a widely-cited textbook figure attributed to Reilly & Burstein 1975, J Biomech 8(6):393-405,
  PMID [1206042](https://pubmed.ncbi.nlm.nih.gov/1206042/) — that citation's existence is verified live,
  but PubMed carries no abstract for this 1975 paper, so the exact GPa figure itself was NOT
  independently re-extracted from the primary source this session, flagged as such): **~34-51 MPa**.
  This session's primary result (54.97 MPa) sits just above this band.
- **Taylor ME, Tanner KE, Freeman MA, Yettram AL. "Stress and strain distribution within the intact
  femur: compression or bending?" Med Eng Phys. 1996 Mar;18(2):122-31.** PMID
  [8673318](https://pubmed.ncbi.nlm.nih.gov/8673318/), DOI
  [10.1016/1350-4533(95)00031-3](https://doi.org/10.1016/1350-4533(95)00031-3). Finding used here: with
  **physiological muscle+joint contact forces, the femoral diaphysis is loaded PREDOMINANTLY IN
  COMPRESSION**, not bending. **This model's own result is bending-dominated at the primary peak
  (axial −3.13 MPa vs bend amplitude 51.84 MPa) — IN TENSION with Taylor 1996's finding, not
  agreement.** Flagged honestly, not smoothed over (§6 discusses candidate reasons: different subject,
  different gait phase/instant than Taylor's own reference posture, a coarse circular-tube section vs
  Taylor's real 3D geometry, and the still-uncorrected rotational-inertia residual).

## 6. Honest gaps (full list, symmetric QC — read before trusting the headline number)

1. **The distal-sum/proximal-sum discrepancy at the primary peak (148.08 N·m) is LARGER than the
   averaged bending moment itself (92.93 N·m)** — ratio 1.40. The two pure-sided bounds are
   **18.81 MPa (distal-only) to 91.56 MPa (proximal-only)** bending-stress-equivalent, bracketing the
   reported 51.84 MPa average. This is the single largest uncertainty in this document and is NOT a
   rounding-level correction — a different (equally defensible) choice of which side to trust could
   move the answer by roughly ±35 MPa. Root cause: the uncorrected rotational-inertia term (§4).
2. **Only translational segment inertia is corrected** (d'Alembert, live COM acceleration);
   **rotational inertia (I·α + ω×Iω) is not** — would need the femur's angular acceleration and
   inertia tensor (both obtainable via the OpenSim API, not attempted this session — the concrete next
   step, §7).
3. **Cortical cross-section is a generic adult literature range, not this subject's own CT/DXA
   geometry** (§2) — a live search this session could not pin an exact peer-reviewed midshaft mm²/mm⁴
   figure behind open access; swept, not asserted, rather than fabricating a citation.
4. **Femoral neck-shaft offset — the single largest lever on the result (§3)** — a literature range
   (30-60mm), not measured on this subject's own femur; direction is empirically derived from this
   model's own abductor geometry, magnitude is not.
5. **Circular-annulus idealization** — the true femoral shaft cross-section is mildly elliptical; this
   also means the peak-stress location around the cortex perimeter is only as accurate as that
   simplification.
6. **Bending-dominated vs Taylor 1996's compression-dominated finding (§5)** — an honest disagreement,
   not resolved this session. Candidate explanations (not adjudicated): different subject/geometry,
   different loading instant (this trial's t=0.56s vs whatever single reference posture/phase Taylor
   used), the residual/rotational-inertia gap (§6.1-2) could plausibly be biasing the bending term up,
   and a coarse beam model cannot capture Taylor's real 3D FE geometry effects.
7. **Coarse Euler-Bernoulli beam vs subject-specific CT-FE** (e.g. Duda et al. 1998's own finite-element
   mesh) — explicitly a first, cheap, fully-auditable layer below the joint-contact-force twin, not a
   replacement for one.
8. **Quasi-static**: even with the translational correction, this is a per-instant snapshot treatment
   (matching Duda 1998's own "four phases of a gait cycle" convention), not a continuous dynamic
   simulation.
9. **This trial's t=0.56s is NOT the same instant as the twin's own reported peak joint contact forces**
   (hip 387.04 %BW @ t=0.55s, knee 391.11 %BW @ t=0.51s, both in the trial's FIRST stance phase) — the
   peak COMBINED bending stress happens to occur at a comparably-loaded SECOND stance phase in this
   two-stride trial instead, because bending depends on the full moment balance (transverse force
   components × lever arms across every muscle/joint load), not on any single joint's own peak contact
   magnitude. This is scanned for and reported honestly, not assumed to coincide with the joint-force
   peaks.
10. **Single trial, right leg, subject2** — same scope boundary as every other MSK cert in this repo.

## 7. Next step

If a tighter number is later needed: correct the rotational-inertia term (angular acceleration +
inertia tensor, both available via the OpenSim API) to close the remaining distal/proximal
discrepancy properly, and re-run the same offset/geometry sweeps to see how much the ~19-92 MPa
bending-stress bound narrows. Until then, the headline number should be read as −55 MPa **±** roughly
35 MPa from the quasi-static residual alone, before even the cross-section/offset literature ranges
(§2-3) are folded in — a first, honestly-bounded layer, not a precision estimate.

## Files

- `scripts/msk/bone_stress.py` — the full pipeline (self-contained, re-runnable; imports
  `validate_joint_force.py` for proven `parse_mot`/`MODEL_FILE`/`IK_MOT`/`G`; reuses the existing
  `static_optimization/{so,jr}` outputs in place, never re-runs the expensive SO/JR solve).
- `data/msk_smoketest/subject2_walking1/bone_stress/bone_stress_results.json` — every number in this
  document, machine-written.
- Reused, not written by this script: `data/msk_smoketest/subject2_walking1/static_optimization/so/
  walking1_StaticOptimization_force.sto` and `.../jr/walking1_JointReaction_ReactionLoads.sto` (both
  pre-existing from the knee/hip contact-force certs; re-verified — §3 of those certs' own headline
  numbers reproduced bit-for-bit here before trusting anything downstream — not re-run).
