# DAY1_REVIEW — BodyTwin, dygn 1 (2026-09-22 ~21:10 → 2026-09-23 ~21:10)

**FINAL VERSION 23/9 22:15 — reviews AU6–AU8 of lane_runner results are ongoing and will be added when they are ready.**

The numbers are copied from `notes/RESULTS_INDEX.md` (line number A.. / B..) or from the specified source file. Where a review has corrected or reformulated, the review's wording is used: AU1 for demo (`results/AU1/AUDIT.md`, `results/D5/DEMO_CLAIMS_v2.md`), AU2 for P3, J1, H7 and AB1b (`results/AU2/AUDIT.md`), AU3 for J2 (`results/AU3/AUDIT.md`), AU4 for JS1 (`results/AU4/AUDIT.md`), AU5 for AT1 (`results/AU5/AUDIT.md`), BT-CERT-AUDIT for LV1/LV2 (`results/BT-CERT-AUDIT/AUDIT.md`).

Status: **A** = the agent who drove has reported; **R** = root has read the source file; **G** = independent review (AU1–AU5, BT-CERT-AUDIT, or field session audit A225 for FV1). **lane_runner, not reviewed** = run or completed by a lane_runner lane; waiting for AU6 (IM1, IM2, FV3), AU7 (LV1, LV2) or AU8 (FV2, H3b, C4, X3c, V7, P2b, C1b, R2).

Version and data for the whole day:

| What | Value |
|---|---|
| Published BodyTwin | `b95a8dc` (untouched; patches in R1/R2 clone) |
| Demo D5 | `results/D5/demo_results.json`, sha256 `3aac5c84…` (AU1 reproducerade bitidentiskt) |
| File index with sha256 and PREREG order | `tasks/index/index.json` (`python3 tasks/index/build_index.py`) |
| Dataset | `tasks/index/DATASETS.json`; tabell sist i `notes/RESULTS_INDEX.md` |

---

## 1. Det viktigaste

**For the collaborator (meeting 2026-09-25)**

Message after reviews AU2 and AU3: the profit lies in task-oriented measures (H7 T5) and in averaged landmarks, not in the population model as such. Against what the reference model can do itself (landmark morph with the hip joint in the caput sphere), P1's five X-ray measurements are not significantly better. RM1 confirms: neither regularized TPS nor SSM priors add anything to the landmark-morph baseline.

1. **Show:** daemon D5 runs the entire chain with one command; AU1 got bit-identical result and 106 of 106 numbers match the source files. Secure according to AU1: B1, B3, B5, B6, B8, B12 (as counted), B13, C3, C6 and hinged nets.
2. **Do not display as stated in v1:** "10 mm anterior → +22 %". Right: 10 mm **posteriorly** +22,0 %, 10 mm anteriorly −2,5 %. B2's "1,3·10⁻¹³ mm" is deleted (0 per construction).
4. **Master model (P3, AU2):** with five X-ray measurements, the population model gives moment arm error 2,5 mm against a TLEM morphed proxy, half the error when scaling with length and epicondyle width (4,4–4,7 mm). A TLEM morph to the individual's 9 landmarks, which the reference model can do, gives the same level (2,3 mm). With sphere-HJC in the landmark-based arm: 4,38 mm, −43,9 %. Band covers 0,895 (requirement 0,90), dropped. Truth is a proxy. **RM1:** affine 2,27/2,26 mm, regularized TPS 2,21/2,19 (equally good), SSM condition 3,71/3,86 (worse) → pre-registration failed.
5. **J1 (AU2):** against the collaborator's scaling, TPS the peak of the hip reaction moves −2,3…+18,7 %. TPS and affine differ systematically in 3 of 5 individuals (same direction in all raters). With a single rater, the spread of the TPS result is as large as the difference, because TPS tracks each landmark exactly.
   **J2 (AU3):** in the collaborator's box lift, the master model is closer to a TLEM-based proxy truth than his ScalingLengthMass. The hip reaction peak has RMSE 9,1 → 4,5–5,4 % of R_J, gmed error 12 → 5–6 % and damage ratio error 0,66–0,77 → 0,31–0,53. The win is borne by three individuals where the collaborator's scaling underestimates the peak by 18–22 %. Per individual, X5 is better with 12 of 19 and T5 with 15 of 19. Against what the reference model can do itself, affine morph to the same 9 landmarks with the hip joint in the caput sphere of the landmarks, X5 is not significantly better (7,0 vs 5,3 %, p 0,10). T5 is (4,5 %, 16/19, p 0,003 vs T1 and 0,045 vs T2). TPS morph on rater means is as close to the proxy as the master model. The hypothesis of discriminability against affine per individual failed (2/5). The proxies are based on TLEM mounts, so the atlas bug is not tested.
6. **H7 measurement proposal (AU2):** if the trochanter top–caput distance and a second length measure replace epicondyle width and caput diameter, the abductor arm error decreases from 2,0 to 1,7 mm, and hip center's ML position relative to the knee from 12,6 to 8,0 mm. Length and inertia become slightly worse. In the collaborator's box lift, these metrics (T5) are the only ones significantly better than the affine landmark baseline morph with landmark-HJC (16/19, p 0,003/0,045; AU3).
7. **Measuring strategy (JS1, reviewed AU4):** With atlas mount with isotropic error σ 8 mm per axis, no measurement set reaches H2's requirements ±3,6 mm in the gmed arm. The floor is ±12,7 mm for CT and up to ±17,8 mm for the others. CT with MR pictured bracket hits the mark. The floor is almost entirely determined by the atlas assumption: 4 mm gives 6,5–14 mm and 12 mm gives 19–23 mm. The number 8 mm is H2's recalculation of Pellikaan 2014 (69 % of the attachment surfaces < 15 mm mean contour distance), not a measured value per axis. AT1 (3 cadavers, thin substrate) found the error anisotropic: the insertion on the trochanter 2 mm on surface registration, origin 12–14 mm mainly anteroposteriorly. Then the floor becomes 5–7 mm for EOS and CT, and the EOS error in the hip reaction becomes ~5 % instead of 7 %. The effect of the moment arm on the hip reaction is in the collaborator's box lift −0,49, not H2's −0,67; it further lowers the EOS/CT errors by 1–2 percentage points. Palpation, video and mocap with functional hip center give greater reaction errors than height/weight even when height/weight is palpated pelvic frame (14,7 vs 17,5 / 18,0 / 22,7 %), but the margin against palpation is small, and mocap with the literature 9,4 mm is no worse. EOS gives the most per cost if EOS costs less than approximately 45–90 % of CT (the breakpoint depends on the scenario); costs are assumed (ordinal scale without source). Absolute damage risk is determined by the ε-N exponent, not by the measurement.
   **H2b (anisotropic atlas error):** full CT ±12,8 → ±5,3 mm (4,2–7,7); EOS 13,3–13,5 → 6,6–7,1 mm; ±3,6 mm can only be accessed with CT + both mounts depicted on 2 mm (3,56). The hypothesis that only the origin needs to be imaged failed (4,0 mm; the arm 3,3× more sensitive to the lateral position of the insertion). EOS never reaches ±3,6.
   **AU5 (A105):** AT1's insertion error (~2 mm) holds; the origin error is definition sensitive (8–19 mm, most comparable 16–19 mm) and is based on a single pure donor comparison; Δr_ab 2,9 mm is not covered (KI 1,5–18). H2b's halved floor therefore rests on a thin base. LHDL may not be used commercially or be in a public repo.

**Anton's Track of Innovation: Motion Determines Geometry** (IM1, IM2; lane_runner, not reviewed — AU6 in progress)

9. **IM1 (in silico, 5 VSD × 64):** regular movement gives almost no information about the mount. Reduction of posterior-SD for gmed r_ab against only landmarks: gait/squat-HJC 5,5 %, isometric moment 4,5 %, GRF+ID 1,3 %. Directed measurements give: stellar motion-HJC (5 mm/axis) 48,5 %, passing 30 % only at 5 mm (9,4 mm → 26,5 %); ultrasound excursion (assumed 1 mm, optimistic) 86 %. The mount needs to be imaged or excursion measured.
10. **IM2 (OpenCap):** common shape and function adaptation against IK+ID residual did not give improvement on held movements and people (pelvic residual +0,05–0,09 %; 0/13, 0/6 0/13) → fell.

**Certified calculations** (LV1, LV2; scope reviewed in BT-CERT-AUDIT, AU7 in progress)

11. **LV1 (w/COM/inertia, H5):** Lean proves real interval composition from already valid term intervals (`MomentChain.lean` compiles without `sorry`); fixed femur has reported full Fraction inclusion and 10 000 accidental without fracture. The general Python/binary64 certificate is conditional (fsum limit, tetrahedron formulas, finite intermediate results and mesh validity pending). **LV1-SUM:** exact rational sum replaces the fsum assumption; 0 fracture; +40,1 % time.
12. **LV2 (case selection for muscle trajectories):** contact case = X1b in 6 989 observations; 0 sign error in 10 000 sample; with forced side 0 torque arm jump > 2 mm/0,5°. At geometry uncertainty u = 12 mm 98,43 % is ambiguous, at u = 19 mm 100 % → arithmetic cannot determine an uncertain mount. Side arch decisions and anatomical accuracy are not certified.


> **AU7 (22:3x, A113):** LV1 holds for interval composition (Lean without sorry); 0 breaks in 10 200 own adversarial cases), but tetrahedron formula, floating-point limits, and coupling Python↔Lean are tested/accepted, and "the error bound governs the work" was not done. LV2: tangent condition for straight path to winding is correct; the choice right/left is an assumed page condition with no source, and "0 jump" is X1b's forced page, not new evidence.


> **AU8 (22:4x, A114):** all eight lane_runner lanes reproduced exactly. Confirmed: FV2, H3b, C4, V7 (reproduction of G1b), P2b (negative). Subject to: X3c (the damage increase requires parameters outside of PREREG and reverses at day 365), C1b (code changed after the counter samples; overridden 9/9), R2 (at least three cells misclassified as B).

**For BodyTwin wide**

13. **Contract C1/C1b:** 9 of 9 counter samples are captured. +10 mm femur length → 2,43 N·m, of which 2,13 from the position of the knee. C1b adds covariance and provenance key: double counting, undeclared correlation, wrong covariance unit, and H6 refusal is stopped (lane_runner, not reviewed). **The fusion cell H6** lowers the percentage of sure wrong decisions from 0,515 to 0,025 (with CT) and from 0,775 to 0,039 (without CT), but refuses wrong person, 15° anteversion and +40 mm video bias only in 0,037–0,636 of the cases (requirement 0,90).
14. **Error budget (C2, C3, C4):** knee position dominates over leg inertia 6,5–10,7×. With measured σ, the hip center term falls to the reference floor (2,85 N·m). In the whole chain shape parameter → damage (C4) is the order relation > slope (ε-N) > geometry in all — the geometry never largest (lane_runner, not reviewed).
15. **Video:** joint center from video is ~26 mm after best correction (V5); F1's requirement 3,8 mm is not reached. Physics constrained COM gives walk-GRF horizontally −60 % (V3). Our own video chain is 3,1× worse than OpenCap in the lead center (V4). Body mass out of visual hole fell: median +44 % (V6). Bias-corrected mirror-TO holds for departure (median 20,7 ms) and clock (V7, lane_runner, not reviewed).
16. **Tissue:** at 2 mm 41–48 % of cortex < 1 mm disappears silently (K2). The laminate dropped where the radius of curvature < cell (H3), and curvature adaptive refinement also dropped (H3b 2/4; lane_runner, not reviewed). On real CT, the relation selection stands for density→modulus of 71–81 % of the variance, phantom-free calibration is rejected as insensitive, and no method measures cortex < 2 mm right (K3).
17. **Mechanobiology:** the unloading × anti-RANKL interaction is +12,42 pp BV/TV with the same sign in 1280 of 1280 runs, but the size is not validated (X3). The FREEDOM case c2 passes via mineralization: BMD +9,81 % vs BV/TV +2,24 % (X3b). Antiresorptive + load deals more damage than load alone in the basic setting, but opposite sign in 99/640 susceptibility tests (X3c, not robust; lane_runner, not reviewed).
18. **Foot contact:** anatomical contact points give GRF peak 87,2 % within ±15 % but HS 37,0 % within ±20 ms (M4). Foot–floor conditions in IK fell (FV2: foot height within ±5 mm 75,12 %, HS 0/50). Shoe model dropped (FV3: HS 0/81, TO 31/78, GRF top 0/39). A new hypothesis is needed.
19. **Published Code:** R1 (5 patches; 777/3/13 → 855/2/13) and R2 (82 fake gates classed A/B/C, only class A changed, commit ee07a39 in the clone). No pushing. **Incident A104:** H2b overwrote `results/H2/results_run.json`; the file is recreated with the same key but not byte-identical, and the graph proposal artifact for the file is out of date.


---

## 2. Track A — the collaborator

Figures: [D5 morph before/after](../results/D5/fig1_morph_before_after.png) · [D5 TPS/affine versus noise](../results/D5/fig2_tps_affine_noise.png) · [D5 error budget with UNKNOWN](../results/D5/fig3_error_budget_unknown.png) · [D5 P1 few measurements](../results/D5/fig4_p1_few_measures.png) · [P3](../results/P3/fig_p3_summary.png) · [J1 z001](../results/J1/fig1_z001_hip_reaction.png) · [J1 individuals](../results/J1/fig2_subjects.png) · [H7 LOSO](../results/H7/fig_loso.png) · [H7 task versus variance](../results/H7/fig_task_vs_variance.png) · [K1b sweep](../results/K1b/fig_k1b_sweep.png) · [X1 TPS/affine](../results/X1/fig_t2_tps_vs_affine.png) · [X1 hip centre rules](../results/X1/fig_t3_hjc_rules.png) · [P1 error versus measurement count](../results/P1/fig_error_vs_k.png) · [J2](../results/J2/fig1_j2.png) · [JS1](../results/JS1/fig_js1.png) · [RM1](../results/RM1/fig1_rm1.png) · [H2b](../results/H2b/fig_h2b.png) · [H2b → JS1](../results/H2b/fig_js1_AN.png) · [H6 refusal](../results/H6/fig_h6_refusal.png)

Note: the figures for P3 and J1 show the agents' original comparisons (P3 against two-measure scaling with A6-HJC; J1 with the A6 source). AU2's corrected numbers are only in `results/AU2/audit.json`.

### 2.1 What can be shown

| Result | Number | Source | Status | Review |
|---|---|---|---|---|
| the collaborator's eget problem | individuella modeller "too cumbersome" (NATO STO 2024) | A1; `notes/JOHN_REQUIREMENTS_AND_SOFTWARE.md` | R | — |
| Demo D5, a command | 8/8 tests; 127 tracked numbers in v2 | `results/D5/README.md`, `DEMO_CLAIMS_v2.md` | A, G | AU1: reproduced bit identically, 106/106 (v1) |
| Certificate of femur and pelvis | both approved; landmarks ≤ 1,20 mm from the surface | D5 B1 | A, G | AU1: confirmed |
| Mounts as triangle + weights | follow the surface through each morph with the same ID; 209 of 232 anchored | D5 B2 | A, G | AU1: the number 1,3·10⁻¹³ mm crossed out (tautology); ID confirmed |
| Rotational invariant folding control | 0 folded under rigid rotation in 80/80; 38/38 VSD morphs; local inversion is found after 126° | A78; D5 B3 | A, R, G | AU1: confirmed |
| Affine mapping misses hip center | 6,80 mm vs 2·SD noise 0,98 mm (z001); 19 of 19 over the noise | A55; D5 B4 | A, R, G | AU1: “Affine least-squares mapping misses the individual's hip center with 6,8 mm (19/19 over the noise). TPS goes by construction through it." |
| TPS vs affine, GMED arm | median difference 2,85 mm vs 2·SD 2,84 mm; 10 by 19 | D5 B5 | A, G | AU1: confirmed |
| Hip center out of caput sphere against regressions | Harrington 12,1±4,2 mm; Bell 15,3±6,7; 6 Point Sphere 0,6 | A56; D5 B8 | A, G | AU1: confirmed |
| the collaborator's hip reaction without runtime | median 0,012 % (p=2) | A36; D5 B10 | A, R, G | AU1: "With the collaborator's geometry and other muscles fixed, our quadratic recruitment reproduces his hip reaction to 0,01 % (the pre-registered cubic gave 3,3 %)." |
| Joint center → hip reaction (box lift) | 10 mm posterior +22,0 %; 10 mm anteriorly −2,5 %; 10 mm lateral +17,6 %; ~1,2–1,5 %/mm | A37; D5 §D | A, R, G | AU1: direction corrected in v2 |
| Femur from few measurements (P1) | 5 X-ray measurements: 22 % (KI 16–28 %) lower surface error than longitudinal scaling; against affine 12,6 % (KI 8,5–16,6 %) | A72; D5 B11 | A, R, G | AU1: "The requirement 20 % is passed in the point estimate." |
| Master model → torque arm (P3) | r_ab RMSE 2,46 against 4,68 mm (−47,5 %), KI −60…−26 %; carried by 5 individuals (without them −11 %) | A89, A95 | A, G | AU2: "With five X-ray measurements, the population model gives moment arm error 2,5 mm against a TLEM morphed proxy, half the error when scaling with length and epicondyle width (4,4–4,7 mm). A TLEM morph to the individual's 9 landmarks, which the landmark baseline can do, gives the same level (2,3 mm). |
| P3 with sphere-HJC (what the landmark baseline uses) | B 4,38 mm; −43,9 % (T1), −40,9 % (T2, p 0,06 two-sided) | A95; AU2 | G | AU2: P3 of the text "HJC = TLEM A6 (what the landmark baseline uses)" is wrong |
| Master model → hip reaction (P3) | RMSE 8,1 %BW vs. 17,2 %BW; TLEM to 9 landmarks 8,5 (affine) / 8,9 (TPS) %BW | A90; AU2 | A, G | AU2: "Enter along with the line above, not as a win against the reference model in general." |
| Where P3's profit comes from | — | AU2 | G | AU2: "The gain comes from using the individual's hip center and trocars, either via measurements or landmarks. Two-dimensional scaling can't do that." |
| Master model in the collaborator's boxlift (J2), error against proxy truth | RMSE R_topp (% R_J) T1: the collaborator 9,1 / affine 7,8 / TPS 6,3 / X5 5,3 / T5 4,5; T2: 8,7 / 7,5 / 7,0 / 5,4 / 5,0; gmed error 12 → 5–6 %; damage ratio error 0,66 → 0,31–0,36 (T1), 0,77 → 0,52–0,53 (T2) | A97, A99 | A, G | AU3: reproduced (≤ 5e-5 relative). Against the collaborator's scaling: X5 better 12/19 (T1, p 0,06), 9/19 (T2); T5 15/19 (p 0,002); the profit is borne by z035/z036/z061 (the collaborator −18…−22 %). "Systematic −4,4 %" crossed out: bias 24 % of the squared error |
| J2 to affine landmark baseline (AU3's missing arm) | affine with hip joint in landmarks caput sphere 7,0/6,5 % → X5 5,3/5,4 not significantly better (p 0,10/0,15); T5 4,5/5,0 significantly better (16/19, p 0,003/0,045) | A99 | G | AU3 |
| TPS on averaged landmarks | TPS on mean of 20 assessments 3,2 % vs. T2 (best with same input; X5 4,9, T5 4,4) | A99; `results/AU3/AUDIT.md` | G | AU3: the master model beats TPS only by a single assessment |
| Measurement suggestion (H7) | abductor arm 2,0 → 1,7 mm; hc_ML 12,6 → 8,0 mm; without hc_ML −1,5 % (p 0,46); i the collaborator's box lift T5 16/19 better than affine landmark baseline (p 0,003/0,045) | A88, A95, A99 | A, G | AU2: "Avoid the balanced −15,6 %." Wording in §1 point 6 |
| Measurement strategy for the collaborator (JS1) | gmed arm ±95 % with atlas attachment: height/weight/palpation/mocap/video 16,6–17,8 mm; EOS 13,3–15,3; CT 12,7; CT + MR-imaged attachment 3,6 (exactly at requirement ±3,6); hip reaction median/p95: first group 11–23 / 33–123 %, EOS 7–8 / 21–23, CT 6,4 / 19, CT+MR 2,1 / 6 %; injury factor p95 6–165 / 3,4–3,9 / 2,8 / 1,4 | A98, A101 | A, G | reviewed AU4; the floor depends on atlas-σ (AT1: anisotropic, limited supporting material); e box lift −0,49 → EOS 6,0 %, CT 4,7 %; height/weight with palpated frame 14,7 % (palpation 17,5, video 18,0, mocap 22,7; mocap with the literature's 9,4 mm 13,9 % no worse); costs without a source |
| Anisotropic atlas error in H2 and JS1 (H2b) | atlas contribution to r_ab 6,47 → 2,65 mm; full CT ±12,8 → ±5,3 (4,2–7,7); ±3,6 mm only with CT + both mounts depicted on 2 mm (3,56); JS1: EOS 13,3–13,5 → 6,6–7,1 mm, hip reaction median 7 → 5 %; CT 12,7 → 5,4 mm, 6,4 → 2,9 % | A103 | A | is based on AT1 — origin can be larger, Δr_ab validation not covered (A105) |
| Walk-ID through the contract | peak moment 48,03 N·m; Δτ 2,43 = knee 2,13 + leg 0,30 | A51; D5 B12 | A, R, G | AU1: status "counted, not compared" |
| UNKNOWN stoppar kedjan | `ChainStop`, inget tal | D5 B13 | A, G | AU1: "kontrollerat beteende" |

### 2.2 What is not suitable to show, or only with reservations

| Result | Number | Source | Status | Review |
|---|---|---|---|---|
| Hamstrings/Quadriceps TPS vs affine | distinguishable 7/19 respectively 0/19 | A55; D5 B6 | A, G | AU1: confirmed (fell) |
| Surface fit thigh | TPS better 12 of 19, below the requirement 2/3 | A55; D5 B7 | A, G | AU1: "population first" |
| "The master model halves the error against the landmark baseline" | TLEM fitted to 9 landmarks 2,27 (affine) / 2,33–2,60 mm (TPS) against A 2,42–2,46; not significant | A95 | G | AU2: the gain only applies to two-dimensional scaling |
| Regularized morph / SSM prior (RM1) | r_ab RMSE T1/T2: affine 2,27/2,26 mm; regularized TPS 2,21/2,19 (difference −0,06 [−0,58; 0,41], equally well, SD 1,56 > affine 1,41); SSM condition 3,71/3,86 (+1,44 [0,90; 1,93], worse) → the pre-registration FAILED | A102 | A | box lift against T2: regularized TPS 4,8 % but noise 3,6 % against affine 1,6; T5 4,5–5,0 % with 2,0 % noise best error/noise ratio |
| J1 "raw TPS indistinguishable from affine (0/5)" | formally dropped; with rectified HJC and paired spread pass the 3/5 gate; corrected span −2,3…+18,7 % | A92, A95 | A, G | AU2: "TPS and affine differ systematically in 3 of 5 individuals (same direction in all raters). With a single rater the spread of the TPS result is as large as the difference, because TPS follows each landmark exactly." The comparison is affine on the same 9 landmarks, not the collaborator's ScalingLengthMass |
| J2 distinguishable vs. affine | different vs. affine > 2·SD 2/5 (X5, T5), TPS 1/5; outside gate 11/19 (X5), 6/19 (T5); over lift 0/5 → H FALL | A96 | A, R, G | AU3: confirmed (gate follows PREREG) |
| H2b: only origin needs to be imaged | origin imaged + insertion from atlas 4,0 mm → FALL | A103 | A | arm 3,3× more sensitive to the lateral position of the insertion |
| P3 band | coverage 0,895 (170/190), requirements 0,90 → dropped | A89 | A, G | AU2: confirmed |
| Tibia from few measurements (P2, P2b) | P2 both scenarios failed; length takes 96 %. P2b joint surface measurements: plateau −0,6 % against B0 / −5,3 % against B1 (p 0,76); ankle +0,3 / −4,4 (p 0,60); both −4,9 / −9,9 (p 0,88) — all FALLS | A83, A107 | A; P2b lane_runner, not reviewed | — |
| Knee contact force against implant (X1b Test 4) | e_DM 465 N > d_gen 90 N; underestimate 346 N | A80 | A, G | AU1: "The test does not decide the geometry question. The simplified The SO pipeline underestimates the power by 346 N." |
| Movable knee + wrapping (X1b) | 10/11 within 10 mm (hinge 7/11) | A79 | A, G | AU1: "With OpenSim's knee kinematics transferred to TLEM: 10/11 within 10 mm of the same model." (partially circular) |
| A measurement determines the GMED arm (H2) | requirement ±3,6 mm; best individual ±7,5 mm; mount must be imaged | A69; D5 C4 | A, G | AU1: based on assumed atlas error σ 8 mm/axis; The D5–H2 correspondence is not independent; AU4: 8 mm is H2's recalculation of Pellikaan, not measured; AT1 found error anisotropic (A100) |
| Damage (F1) | exponent-SE → factor 3,5 in D; joint center 2 mm → 1,6× | A46; D5 C5 | A, G | AU1: "The exponent dominates. The damage/reaction ranking is almost the same (ρ 0,93)." |
| Straight muscle path as error source (K1) | +2,2 % (prediction ≥ 20 % fell) | A39 | A | — |

PREREG order for P3, H7, J1 (AU2): mtime for PREREG is older than first output, but there is no sha256 before run and the directory is not in git. J2 (AU3): mtime + sha correct, no git.

### 2.3 Left in track A

| ID | Question | Status |
|---|---|---|
| AB1b-runtime | the reference model calculates the predicted thigh (≤ 1e-6 m) | waiting for license/Wine (Anton) |

---

## 3. Track B — BodyTwin wide

### 3.1 Geometry, resolution and contract

Figures: [B1 error against pitch](../results/B1/fig/fig1_error_vs_pitch.png) · [B1 refinement](../results/B1/fig/fig3_K3_refinement.png) · [H3](../results/H3/fig_dq3_dq1.png) · [H3b cost/error](../results/H3b/fig_cost_error.png) · [H6 refusal](../results/H6/fig_h6_refusal.png)

| Result | Number | Source | Status |
|---|---|---|---|
| 2 mm is enough for mass/COM/axle | 15 of 19 combinations at 2 mm; surface defined moment arm requires 1 mm | A30 | A, R |
| Required pitch predictable | 19 of 19 within one step (15 exact) | A32 | A, R |
| Voxel 2 mm vs exact inertia | In 0,06–0,17 %, mass ≤ 0,22 % | A10 | A |
| Exact second moment | error ≤ 5,4e-13 | A9 | A, R |
| Thin cortex at 2 mm | 41–48 % of cortex < 1 mm disappears silently, volume error ≤ 0,74 % | B41 | A, R |
| Laminate: cortex volume and tunnels | cortex volume 0,0 %; tunnels 95–100 % discovered | B65 | A |
| H3b curvature adaptive refinement | A: stiffness+strain ±5 % on 4 test tibia FALL (2/4; +10,93 % and +6,27 %); B: ≤ 1/3 element FALL (112802 34,22 %); C random counter sample HOLDS; validation against H3 HOLDS | B92 | A (lane_runner, not reviewed) |
| Composition contract C1 | 6 rules; 9 of 9 counter sample | A49, A52 | A, R |
| C1b: covariance + provenance key + H6 link | C1's 9/9 counter sample left (zero parameter 3,55e-15; hip flexion exact 2,4296 N·m); duplicate count → DuplicateProvenanceError; undeclared correlation → UNKNOWN; wrong covariance unit/shape is captured; H6 refusal stops next link — all HOLD | B96 | A (lane_runner, not reviewed) |
| Osim → pinocchio adapter (M2) | segment positions 1,8e-13 mm RMS; ID ≤ 3,1e-9 | B37, B38 | A, R |
| H6 fusion cell: wrong side / m↔mm | refusal 1,00 (holding); false alarm 0,053/0,005 (holding) | B90 | A, R |
| H6 sure wrong number | naive → cell 0,515 → 0,025 (with CT), 0,775 → 0,039 (without) | B90 | A, R |
| H6 wrong person, anteversion, video bias | refusal 0,636/0,634, 0,068/0,037, 0,421/0,526 (requirement 0,90) FAILED; coverage with CT 1,00, without CT 0,847 FAILED | B90 | A, R |

### 3.2 Tissue, material and damage chain

Figures: [H5b](../results/H5b/fig_h5b_summary.png) · [K3 FE Grid](../results/K3/out/fig_fe_grid.png) · [K3 thickness 006](../results/K3/out/fig_thick_006.png) · [V6 summary](../results/V6/fig_v6_summary.png) · [C4 error budget](../results/C4/fig_c4_budget.png)

| Result | Number | Source | Status |
|---|---|---|---|
| Density→modulus determines mechanostat regime (synthetic density) | subchondral strain 5,9× relationships between; Morgan 443–616 µε vs. Keller/Carter–Hayes 1806–2621 µε | B42 | A |
| No isotropic mixing rule for thin layer | Voigt strain −30,9 %; Reuss stiffness −20,1 % | B43 | A |
| K3 HU → density without phantom | spread max/min trab 1,49/1,67, short 1,21/1,27; without cortex reference cortex 2,25–2,51 (17–31 % over ICRU 1,92); calibration FAILS as insensitive | B83 | A, R |
| K3 T1 on real density | relation selection spread 15,4/13,7×; Λ 0,47/0,71 (intermediate position); relation 71–81 % of the variance; Morgan always below 1500 µε, Carter–Hayes always above, Keller determined by calibration | B84 | A |
| K3 cortical thickness | median \|Δ\| 0,28/0,12 mm but 0,77–1,18 mm in metaphyses; no method measures cortex < 2 mm right | B85 | A |
| C4 shape parameter → hip reaction → neck strain → damage | order relation > slope (ε-N) > geometry of all — geometry never greatest; the hypothesis (slope > relation > geometry) intermediate position; M1 bit-identical null change and M2 UNKNOWN stops FE (beam continues) HOLDS | B93 | A (lane_runner, not reviewed) |
| Contour from video is sufficient if | contour error ≤ ~1,4 mm (mass out of contour) or ≤ ~2,9 mm (known mass) | A86 | A |
| V6 mass from visual hull | 0/9 within 5 % (+28…+50 %, median +44 %) FAILED; 3/4/5 cameras +101/+62/+44 % | B88 | A, R |
| V6 thigh mass and contour | thigh mass versus de Leva −1,2 %BW median; contour error 6,5–11,2 mm → 0/9 ≤ 5 mm FALL | B89 | A |

### 3.3 Mekanobiologi

Figurer: [X3 armar](../results/X3/fig_arms.png) · [X3 Sobol](../results/X3/fig_sobol.png) · [X3b](../results/X3b/fig_x3b.png) · [X3c](../results/X3c/fig_x3c.png)

| Result | Number | Source | Status |
|---|---|---|---|
| Interaction unloading × anti-RANKL | +12,42 pp BV/TV dag 390; same sign 1280/1280 | B6 | A, R |
| Validation X3 | c2 +2,40 % against FREEDOM +9,2 % (fell); c4 29 moon against ~9 (fell); c3 passed | B7 | A |
| X3b secondary mineralization (c2) | BMD +9,81 % vs. BV/TV +2,24 %; 20/20 grid points | B17 | A, R |
| X3b interaction | positive sign 28/28 variants, +1,1–12,4 pp | B21 | A, R |
| X3c damage → targeted remodeling | antiresorptive + load deals more damage than load alone: r1 3,722 dag 120 in default, but opposite sign in 99/640 sensitivity test (not robust); counter sample M1–M3 HOLDS | B94 | A (lane_runner, not reviewed) |
| Decisive measurements | CTX 1–2 weeks after simultaneous start (SNR ≈ 14); P1NP after sclerostin neutralisation 1 and 3 months | B9, B22 | A |
| Synth flag in published code (V1) | is never read; 5 of 10 gates False → exit 0 | B10 | A, R |

### 3.4 Kinematics, dynamics, motor and foot contact

Figures: [C2](../results/C2/fig_c2_summary.png) · [C3](../results/C3/fig_c3_summary.png) · [M1 ID against OpenSim](../results/M1/fig_A1_id_vs_opensim.png) · [M3 HS breakdown](../results/M3/fig_M3_hs_decomp.png) · [M4 breakdown](../results/M4/fig_M4_decomp.png) · [M4 error](../results/M4/fig_M4_errors.png) · [FV2 foot height and HS](../results/FV2/fig_FV2_height_and_hs.png) · [FV3 comparison](../results/BT-FV3/fig_BTFV3_comparison.png) · [FV3 static displacements](../results/BT-FV3/fig_BTFV3_static_offsets.png)

| Result | Number | Source | Status |
|---|---|---|---|
| Engine = OpenSim-ID | RMS ≤ 3e-9 N·m | B11 | A, R |
| Segment masses from a gait cycle | 0 of 12 identifiable; total mass 72,50±0,18 kg | B14 | A |
| Kneeling position against inertia (C2, 19 VSD) | ratio 6,51–10,67 (median 7,91) | A62 | A, R |
| Donor bone | femur length scaling 0,34 N·m; no scaling 5,15 | A63 | A |
| Error budget with same filter (C3) | HJC assumed 4,08; soft part 4,07; STA 3,70; marker 3,70; HJC SCORE 3,25; reference floor 2,85 N·m | B52 | A |
| Knee as hinge | pelvic residual 9,74 % | B38; D5 §D | A, R, G |
| Foot switch, frozen model (M3) | HS +75,7/+81,2/+76,6 ms, 0 % within ±20 | B67 | A, R |
| Zeni marker method for events | HS −1,8 ± 11,6 ms (90 % within ±20) | B68 | A |
| M4 anatomical contact points | HS 37,0 % within ±20 ms FALL; TO 10,3 % FALL; GRF peak 87,2 % within ±15 % HOLDS | B86 | A, R |
| M4 mechanism | heel point 18,1 mm above floor at HS > cushion 10,3 mm; stiff toes in IK (mtp 3°); video-Zeni HS 81,5 / TO 42,3 % | B87 | A |
| FV2 foot-floor conditions in IK | H1 foot height within ±5 mm 75,12 % (requirement 90) FALL; H2 marker residual C/U 1,043 (≤ 1,20) HOLD; H3 contact model's HS within ±20 ms 0/50 FALL; counter sample N1 wrong phase 1,060 (requirement > 1,5) FAIL; HS bias +77,6 → +57,5 ms, foot +13,9 mm above the floor at HS | B91 | A (PREREG) hash OK; lane_runner completed, not reviewed) |
| FV3 shoe model (sole, midsole, markers on shoe, mtp spring) | HS within ±20 ms 0/81, TO 31/78 (39,7 %), GRF top 0/39 → all FALLS (requirement 80/80/70 %) | B98 | A (lane_runner, not reviewed) |

### 3.5 Video

Figures: [E2 hip bias](../results/E2/fig_hip_bias.png) · [V3](../results/V3/fig_v3_examples.png) · [V4 chain](../results/V4/fig_chain.png) · [V5 reduction](../results/V5/fig_v5_reduction.png) · [G2 ρ](../results/G2/fig_rho.png) · [FV1](../results/FV1/fig_reduction.png) · [E1 coverage](../results/E1/fig_coverage.png) · [V6 hull subject2](../results/V6/fig_hull_subject2.png)

| Result | Number | Source | Status |
|---|---|---|---|
| Keypoint bias at hip | median 30,7 mm; correction −21 % (requirement 50 %) | B32, B34 | A, R, G (AU1 C3) |
| Posture Conditional Correction (V5) | 42,0 → 26,5 mm (−36,5 %, requirement 40 %) → dropped; new activity 21,5 % | B81 | A, R |
| Joint center from video, conclusion | model residual ~26 mm + person term ~14 mm vertical; F1's 3,8 mm far away | B82 | A |
| GRF out of physics limited COM (V3) | time horizontal 3,7 %BW (−60 %) holds; DJ 19,7 %BW (claim < 15) fell | B56 | A, R |
| V3 → ID | walk hip 14,1 N·m (V2 77) | B58 | A |
| Gait events from video | HS median 9,0 ms | B47 | A |
| Our chain against OpenCap (V4) | chain center 72,4 against 26,6 mm (3,1×); angle 18,5° to 4,4° | B72 | A, R |
| Camera geometry predicts errors per joint (G2) | dropped: hip ρ −0,02/0,06; uncertainty = empirical constant ~39–57 mm | B70, B71 | A |
| Self Occlusion (FV1) | ρ hip 0,23 (requirement 0,30); weighting 1,0 % (requirement 25 %) | B53, B54 | A, G (A225) |
| Gravity as ruler/clock (G1, G1b) | ruler 4,1 %; watch 100/100/98,3 % holds | B61, B76 | A |
| Bias-corrected mirror-TO, LOSO (V7) | (a) TO error median 20,7 / p95 72,5 ms HOLD; (b) ruler 8/10 ≤ 5 % (median 2,62 %) FALL (requirement 9/10); (c) clock 100/100/98,3 % HOLDS; (d) V3 DJ-MAE 19,51 %BW (KI 15,05–23,41) FALL (requirement < 18,3) | B95 | A (lane_runner, not reviewed) |
| Camera uncertainty → moments (E1) | geometry covered 0,91–0,98; moments fell (corners 0,544–0,676) | B23 | A, R |
| Visual hull (V6) | control gate markers in silhouette 0,94–1,00; mass +44 % median; 1° calibration error shrinks hole → closer to scale (volume criterion rewards error) | B88 | A, R |
| Video inventory (V0, V0b) | 8 candidates; 22 missed clusters | B28, B59 | A, R |

### 3.6 The movement determines the geometry (Anton's innovation track)

Background and plan: "Innovation track from Anton (2026-09-23 ~11:25)".

| Result | Number | Source | Status |
|---|---|---|---|
| IM1 ensemble matching 9 landmarks + OED (5 VSD × 64) | reduction of posterior-SD against only landmarks (gmed r_ab / HJC / hip reaction): ultrasound excursion (assumed 1 mm, optimistic) 86 / 24 / 62 %; star motion-HJC (5 mm/axis) 48,5 / 57 / 46 % — passes 30 % only at 5 mm (9,4 mm → 26,5 %, 15,5 → 12,3 %); walk/squat-HJC 5,5 / 16 / 3,4 %; isometric moment 4,5 / 1,3 / 3,2 %; GRF+ID 1,3 / 9,5 / 0,8 % | A106 | A (lane_runner; PREREG hashed; not reviewed) |
| IM2 common form+functional alignment against IK+ID residual (OpenCap) | kept movements, held persons and LOSO: pelvic residual +0,05–0,09 % (no improvement), 0/13, 0/6, 0/13 → FALL | A108 | A (lane_runner, not reviewed) |

### 3.7 Certified calculations

Background and plan: "Innovation track from Anton (2026-09-23 ~12:40)". Arithmetic guarantee ≠ correct geometry, bracket or page conditions.

| Result | Number | Source | Status |
|---|---|---|---|
| LV1 interval for m/COM/inertia (H5) | Lean proves real interval composition from already valid term intervals; fixed femur has reported full Fraction inclusion and 10 000 accidental without fracture; general Python/binary64 certificate VILLKORLIGT (fsum limit, tetrahedron formulas, finite intermediate results and mesh validity pending); inertial width is not total biomechanical error margin | A109 | G (range/sources reviewed; no independent large overtaking) |
| LV1-SUM exact final sum | fsum assumption replaced by exact rational endpoint sum; 1 964 accepted/2 000, 36 overflow, 0 violation; 43 324 femoral faces with all 13 moments included; full range path 5,61 to 4,00 s (+40,1 %) | A111 | G (root code/statement review and five hash checks; numeric full sample agent run) |
| LV2 certified case selection for muscle courses | contact case = X1b in 6 989 segment–angle observations (6 956 scertain in first step, 33 need second); 0 sign error in 10 000 test against rational oracle; with forced side 0 moment arm jump > 2 mm/0,5° (raw X1b: VAS 12, SM 2, ST 1); at u = 12 mm 98,43 % ambiguous, u = 19 mm 100 % | A110 (lane_runner line) | A (lane_runner, not reviewed) |
| LV2 range | exact rational contact classification conditional on derived local coordinates and specified error bounds; lateral arch trigonometric decision and anatomical accuracy EJ certified; no vested speed gain | A110 (review line) | G (conditional contact result, other claims limited) |

### 3.8 Hands

Figur: [HD1 grepp mot friktionskon](../results/HD1/fig_grip_vs_cone.png)

| Result | Number | Source | Status |
|---|---|---|---|
| Friction cone as lower limit | 99,8 % | B79 | A |
| Grip force from formula | 72 % error; k median 2,1–2,5 (spread 1,4–5,0) | B79 | A |
| OpenSim hand (private) | 45–134 N vs 304–500 N norm | B80 | A |

### 3.9 Published code, graph and incidents

| Result | Number | Source | Status |
|---|---|---|---|
| R1 patch series | 5 patches; 777/3/13 → 855/2/13 (+77 tests); 2 remaining = time limits at load | B75 | A, R |
| R2 false gates with exits 0 | 82 entries graded A/B/C (file:row); only class A changed; mutation test exit 0 → 2; full run same 2 timeout error as R1, no new ones; commit ee07a39 only in clone | B97 | A, R (%B: 0 attribution); lane_runner, Unreviewed |
| Graph proposal GR1/GR2 | GR2 covers A1–A100, B1–B90 (65 folds, 1319 artifacts, 1765 links, 32 gaps); built before A101–A111, B91–B98 | `tasks/graph_proposals/README.md` | A, R |
| Incident A104 | H2b overwrote `results/H2/results_run.json` via `from h2_core import *`; recreated, key figures match H2's README, but not byte-identical (hash f721c62d… → 9fbbbd68…); graph_proposals artifact for file is out of date | A104; `results/H2b/README.md` §Deviations | R |

---

## 4. Negative and tied results

| Result | What fell | Source | What it sharpens |
|---|---|---|---|
| P3 versus the landmark-morph baseline | TLEM to 9 landmarks 2,27–2,60 mm versus A 2,42–2,46; not significant | A95 | the contribution is information on proximal geometry (hip centre, trochanters), not the atlas; demo must compare against landmark morph, not just two-measurement scaling |
| RM1 prior-morph | regularized TPS as good as affine; SSM condition worse (+1,44 mm) | A102 | no geometric operation with measurable gain in r_ab which the landmark baseline lacks; the robust is noise weighting and task-oriented measures (T5) |
| P3 band | coverage 0,895 < 0,90 | A89 | band needs ×1,18; misses at the edge of the population |
| JS1 atlas floor (AU4) | floor is determined by atlas-σ: 4 mm → 6,5–14 mm, 12 mm → 19–23 mm; 8 mm is recalculation, not measurement | A101, A100 | AT1 (thin substrate): insertion ~2 mm, origin 12–14 mm → floor 5,3–13,5 mm; more cadavers with pictured mounts needed (AT1 §5) |
| AT1 (AU5) | origin definition sensitive 8–19 mm; a pure comparison; Δr_ab 2,9 not paved | A105 | H2b's floor rests on thin substrate |
| H2b only origin | 4,0 mm > 3,6 | A103 | both mounts must be imaged; EOS never reaches ±3,6 |
| J2 gate | 2/5 over the judge noise (11/19 outside gate) | A96, A99 | against affine landmark baseline with landmark-HJC is X5 not significantly better (p 0,10/0,15), T5 is (p 0,003/0,045); TPS on assessor means 3,2 % → profit sits in task-oriented measures and averaged landmarks; atlas error untested → depicted mounts |
| J1 gate | formal 0/5; with corrected HJC and paired spread 3/5 | A92, A95 | one assessor not enough: mean of 2–4 assessor or regularized/surface-fit morph required |
| H7 combined measure | without hc_ML −1,5 % (p 0,46) | A95 | the measurement proposal must be specified per quantity |
| IM1 ordinary movement | gait/squat-HJC 5,5 %, GRF+ID 1,3 % reduction in r_ab-SD | A106 | attachment from movement requires targeted measurements: star movement (≤ ~9 mm HJC error) or tendon excursion |
| IM2 shape+function | 0/13, 0/6, 0/13 | A108 | dynamic consistency adds no geometry information in this setup (cf. IM1: GRF/ID ~1 %) |
| LV2 uncertain attachment | 98,43 % ambiguous at u = 12 mm, 100 % at 19 mm | A110 | the certificate determines case selection only when the geometry is known; the bracket must be measured first |
| P2/P2b tibia | all scenarios | A83, A107 | neither angles nor joint surface shape carry surface information beyond affine on tibia |
| W1 brackets/wrap | difference within the noise 0/19 | A58 | whichever way is correct requires mounts pictured |
| X1b Test 4 | 465 N > 90 N | A80 | the test does not determine the geometry; recruitment must be resolved first |
| H2 individual measurement | best ±7,5 mm against requirements ±3,6 | A69 | bracket must be imaged |
| F1 ranking | ρ 0,93 | A45 | the exponent dominates → ε-N at Δε 0,002–0,003 |
| C4 hypothesis | slope > relation > geometry: intermediate position; relationship first | B93 | the density–modulus relationship, not the geometry, is the biggest damage uncertainty |
| K1 straight path | +2,2 % vs predicted ≥ 20 % | A39 | wrapping is not K1's biggest error source |
| B1 first order / directed refinement | 49 % (requirement 90); 0 of 24 | A31, A33 | second order (curvature) needed |
| H3/H3b laminate and curvature adaptive refinement | H3 0/5; H3b 2/4 and 34,22 % elements | B64, B66, B92 | curvature refinement is not enough for ±5 % under 1/3 | cost.
| K3 phantom-free calibration | insensitive to wrong reference tissue (−31/−30 % < method spread) | B83 | the air point controls the slope; tight anchor point required |
| K3 T1 | Λ 0,47/0,71, hypothesis Λ ≥ 1 hnot valid | B84 | the ratio selection, not the calibration, determines the mechanostat side |
| K3 thickness | cortex < 2 mm not measured correctly | B85 | requires defolding or HR-pQCT |
| C2/C3 filter | 4 Hz moves top 2,24 > gain 1,78; same filter −0,2 % | A67, B49 | filter limit is own unknown term |
| C3 functional hip center | lateral indeterminate (47,2 mm) | B51 | requires star movements (cf. IM1) |
| H6 refusal | wrong person, anteversion, video bias < 0,90; coverage without CT 0,847 | B90 | the contract needs covariance and provenance key → C1b |
| X3/X3b validation | c4 not validated; PK fell | B7, B18, B19 | decisive measurement: P1NP after sclerostin neutralization |
| X3c robustness | opposite sign in 99/640 | B94 | antiresorptive × load is not sign-stable in sensitivity space |
| M2/M3/M4/FV2/FV3 foot switch | M3 0 %; M4 HS 37,0 %; FV2 HS 0/50, foot height 75,12 %; FV3 HS 0/81, GRF top 0/39 | B39, B67, B86, B91, B98 | neither anatomical points, foot–floor condition or shoe model solves HS; new hypothesis needed |
| E2/V5 video correction | 21 % resp. 36,5 % | B34, B81 | joint center must be imaged |
| V6 mass out of hole | 0/9 within 5 %; contour error 0/9 ≤ 5 mm | B88, B89 | frontal cameras overestimate volume; the volume criterion rewards calibration errors |
| G1/G1b/G2/FV1/V7 | ruler, cert, mirror-TO, per-joint prediction, occlusion; V7 ruler 8/10, DJ-MAE 19,51 | B61–B63, B76, B70, B53, B95 | detector bias dominates and is not predicted. out of geometry; the clock keeps |
| V4 own video chain | 3,1× worse than OpenCap | B72 | replacement plan per component (B73–B74) |
| HD1 grip force | 72 % error | B79 | k (safety margin) dominate, not pose |

---

## 5. Unknown — requires the collaborator or new data

| Question | Why | Source |
|---|---|---|
| Target organs, left. (AnyBody-version, licensserver), input, output format, acceptance dimensions, API | determines what demon should optimize | `notes/JOHN_REQUIREMENTS_AND_SOFTWARE.md` §4, questions 1–10 |
| How does his workflow place the hip center (frames, regression, functional, image)? Does he use landmark morph or scaling with few dimensions? | determines the relevance of B4 and of P3 (gain only applies to two-dimensional scaling) | D5 v2 §F question 8; A95 |
| Can the femur population from the osteotomy work be shared? | P3 based on Imperial 70 | ask 6 |
| external solver runtime under Wine; node locked trial? | AB1b cannot be run without | AB1b |
| Illustrated muscle mounts (both gmed mounts) | only independent reference for W1/C2/H2; H2b: ±3,6 mm requires both | A69, A58, A103 |
| Senexkursion (ultraljud) med verklig noggrannhet | IM1 antog 1 mm (optimistiskt) | A106 |
| Joint center from image (MR/CT/EOS) or star movements | mocap-HJC separate 22 mm between definitions; star motion helps at ≤ ~9 mm | B35, A106 |
| Multiple raters per individual | J1: one rater gives TPS spread equal to the difference | A95 |
| ε-N-kurva vid Δε 0,002–0,003 | exponenten dominerar skadan | A47 |
| Dense density reference in CT; validated density–modulus per site | calibration lacks anchor at dense end; the relationship determines the mechanostat and the damage chain (C4) | B83, B84, B93 |
| Shoe information (sole thickness, midsole) and why HS is missed | FV3 fell with shoe model | B98 |
| Recruitment Criterion | dominates damage and knee contact force | A38, A80 |
| AT1 query 1: does the group have multiple MR- or CT based individual models where gmed/gmin mounts are imaged or segmented, and not morphed from TLEM? | the atlas error is tested on 3 cadavers | `results/AT1/README.md` §5 |
| AT1 question 2: why were they moved 12 `GluteusMedius*Node` 5,9–31 mm (centroid 9,4 mm) between TLEM2.1 and 2.2? Which source applies? | determines which attachment version is reference | AT1 §5; A28 |
| AT1 question 4: is the Pellikaan's other carcass available? | more independent comparisons of the atlas error | AT1 §5 |

---

## 6. Recommended sequel

| # | Alternativ | Evidens | Kostnad |
|---|---|---|---|
| 1 | **Mount imaged or measured via excursion.** Data with imaged gmed mounts (both) from the collaborator's group, TLEMsafe, or the LHDL entity; examination of tendon excursion with true ultrasound accuracy; star movement-HJC | JS1/AU4: only CT + MR image mount reaches ±3,6 mm; H2b: both mounts required (3,56 mm), only origin 4,0; AT1/AU5: a pure comparison; IM1: excursion 86 %, star movement 48,5 % at 5 mm, passage 5,5 % | data and agreement (AT1's questions 1–4 to the collaborator); the calculation is available |
| 2 | **Task-oriented measures and averaged landmarks** as message and demo for the collaborator; always compare against affine landmark baseline with landmark-HJC | H7 T5 16/19 (p 0,003/0,045, AU3); TPS on means of 20 assessments 3,2 % against T2; RM1: prior-morph adds nothing; AU1-safe lines in D5 | rewording + an extra arm in existing code; multiple raters require data |
| 3 | **Certified calculations** further: bridge code ↔ Lean for LV1, C1 chain certificate, and LV2 only where geometry is known | LV1 Lean without `sorry`, 10 000 case without break, general certificate conditionally; LV1-SUM 0 offense (+40,1 % time); LV2 0 sign error but 98,43–100 % ambiguous at u 12–19 mm | lane_runner lanes; AU7 reviews first |
| 4 | **Foot contact and shoe: new hypothesis needed.** Anatomical points, foot–floor conditions in IK and shoe model have all fallen on HS | M4 HS 37,0 %; FV2 HS 0/50, foot +13,9 mm at HS; FV3 HS 0/81, TO 39,7 %, GRF top 0/39; M4: stiff toes in IK (mtp 3°) | new PREREG before driving; AU6/AU8 reviews FV2/FV3 first |
| 5 | external solver runtime according to AB1b; R1/R2 into staging after Anton's decision; tissue chain via validated density module and tight density anchor | AB1b 7 includes static ok, HJC fix confirmed (AU2); B75, B97; K3 relationship 71–81 %, C4 relationship first | Anton's license/Wine; decision; density anchor requires data |

---

## 7. How to drive

**Demo (D5)**

```bash
cd 
python3 results/D5/demo/run_demo.py --subject z001          # hela kedjan → demo_results.json, run_log.json, graph_node_z001.json
python3 results/D5/demo/figures.py --subject z001           # fig1–fig4
python3 results/D5/collect_claims_v2.py                     # claims_v2.json
cd results/D5 && python3 -m pytest -q -p no:cacheprovider tests   # 8 tester
```

Claims: `results/D5/DEMO_CLAIMS_v2.md`. The dynamics step runs in `../.venv-motion/bin/python`.

**The the reference model Attempt (AB1b)**

```bash
cd results/AB1b
nice -n 19 python3 -m pytest -q -p no:cacheprovider tests          # 34 tester
for t in identity z001 p1_z001; do nice -n 19 python3 validate_static.py --subject $t; done
WINEPREFIX=$HOME/.wine nice -n 19 python3 run_anybody_dump.py --subject off --dataset TLEM2.2       # [LICENS][RUNTIME]
WINEPREFIX=$HOME/.wine nice -n 19 python3 run_anybody_dump.py --subject identity --dataset TLEM2.2
```

The rest of the steps and the comparison (`compare_runtime.py`): `results/AB1b/TOMORROW_v2.md`.

**lane runner-lanes**

```bash
cd 
tasks/lanes/run_lane_runner.sh BT-IM1                      # one lane: reads tasks/lanes/BT-IM1.md, writes results/BT-IM1/RESULTS.md
tasks/lanes/run_codex.sh BT-FV2 results/FV2          # lane med annan resultatmapp
MODEL=lane-model EFFORT=medium MAX_ROUNDS=3 tasks/lanes/run_codex.sh BT-V7   # default values, may be changed
PAR=3 tasks/lanes/run_lane_runner_queue.sh                 # queue in priority order, at most PAR at once
```

Lane description: `tasks/lanes/<LANE>.md`. Log: `tasks/lanes/lane_runner_<LANE>.log` and `tasks/lanes/lane_runner_queue.log`. The script runs `lane_runner exec` with workspace-write-sandbox (write-roots sdc1/bodytwin, bodytwin_work, /tmp), 2 threads, `nice 10`, max 3 rounds, and is finished when `results/<LANE>/RESULTS.md` starts with `# <LANE>`. A lock per result folder prevents double-start.

**Granskningarna**

| Granskning | Vad | Fil | Skript |
|---|---|---|---|
| AU1 | D5-demon | `results/AU1/AUDIT.md` | `results/AU1/scripts/` |
| AU2 | P3, J1, H7, AB1b | `results/AU2/AUDIT.md`, `audit.json` | `results/AU2/scripts/` (t.ex. `p3_arms.py main`, `h7_audit.py agg`, `j1_gate.py`, `ab1b_hjc_check.py`) |
| AU3 | J2 | `results/AU3/AUDIT.md`, `audit.json` | `results/AU3/scripts/` (`au3_j2.py repro`, `au3_j2.py arms`, `au3_stats.py`) |
| AU4 | JS1 | `results/AU4/AUDIT.md` | `results/AU4/scripts/` (`au4_k1.py spot/elast/row1p`, `au4_atlas.py`, `au4_combine.py`, `au4_mcnoise.py`) |
| AU5 | AT1 | `results/AU5/AUDIT.md`, `audit.json` | `results/AU5/scripts/` |
| BT-CERT-AUDIT | LV1, LV2 (reachability) | `results/BT-CERT-AUDIT/AUDIT.md`, `audit.json` | — |
| AU6 | IM1, IM2, FV3 | in progress | — |
| AU7 | LV1, LV2 | ongoing | — |
| AU8 | FV2, H3b, C4, X3c, V7, P2b, C1b, R2 | ongoing | — |

The probes run with nice 19, 2 threads and 1 process, reading the source folders without changing them and only writing in their own `results/AU<n>/`.

**R1-/R2-patchar**

R1: clone `external_mount`, branch `fix/day1-findings`; patches `results/R1/patches/0001–0005*.patch`, applied with `git am` on top of `b95a8dc`; test logs before/after in `results/R1/`. R2: `results/BT-R2/RESULTS.md` and `results/BT-R2/patches/` (commit ee07a39 only in the clone).

**Graph Suggestion**

```bash
nice -n 19 python3 tasks/graph_proposals/build_proposals.py
nice -n 19 python3 tasks/graph_proposals/validate.py
```

Loading: `tasks/graph_proposals/README.md`. No generated graphs are changed. The artifact for `results/H2/results_run.json` is out of date (A104).

**Other results:** run command is in `results/<ID>/README.md` or `RESULTS.md` respectively.
