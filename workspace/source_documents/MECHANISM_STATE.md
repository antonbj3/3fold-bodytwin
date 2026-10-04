> **STANDING DIRECTIVE 2026-07-25 — READ `docs/MECHANISM_COLLISION_COURSE.md` FIRST.**
> The audit menu in the tick prompt is maintenance, not the mission. Aim the stress-point method at
> paralysis, Alzheimer's and cancer, and state findings as claims a clinician could argue with —
> not as internal audit results. Same rigour, bolder aim.

# MECHANISM — STATE OF THE CERTIFIED TWIN (2026-07-18)

## ⚠ STALE / SUPERSEDED (dated 2026-07-18 — the oldest doc in this family; flagged 2026-07-21)

This doc's "~92 designed cells are EXECUTED / 20 re-verified" claim (line below) **cannot be
reconciled against the current live graph** and should be read as stale, not current. Re-verified
directly this session: `data/MECHANISM_ANCHOR_GRAPH.json` currently has **998 total nodes**, status
breakdown `{OPEN: 923, ASSUMED: 66, REFUTED: 5, DEFERRED: 3, PROVEN: 1}` — its own `_meta.note` says
plainly **"status OPEN = designed, NOT measured."** Neither the "~92"/"20" figures nor this doc's own
tier vocabulary (PROVEN/MEASURED-B/MEASURED-B-partial/NEGATIVE-CERTIFIED) line up with that graph's
schema, and `bt_memory/LEDGER.jsonl` is an event stream, not a per-cell tally, so it cannot arbitrate
either. Note also (per `docs/MECHANISM_TRUST_LEDGER.md` §11): `data/MECHANISM_ANCHOR_GRAPH.json` /
`docs/MECHANISM_CELL_INVENTORY.md` are a **different, broader whole-body design catalog** (not
MSK-specific) from whatever this doc's own "~92 cells" figure was originally counting — so this is
most likely two different tracking schemes that were never the same graph, not a regression in
either one. This is flagged **stale/likely-superseded, not asserted false** — no attempt was made
this session to reconstruct what the original "~92/20" count referred to. Full disposition:
`docs/MECHANISM_TRUST_LEDGER_REMEDIATION.md`. The rest of this document is preserved verbatim below
(historical record) — do not delete.

---

What the twin can *certifiably* do, organized by the value proposition the session's own self-tests established.
Source of truth = `data/MECHANISM_ANCHOR_GRAPH.json` + `bt_memory/LEDGER.jsonl`; this file is the legible map.
Status tiers: PROVEN > MEASURED-B > MEASURED-B-partial; NEGATIVE-CERTIFIED / MEASURED-NEGATIVE = honest negatives
(a PASS in this project's discipline). **All ~92 designed cells are EXECUTED.** A 2026-07-19 ADVERSARIAL-HARDENING
phase re-verified 20 of the strongest-headline cells (reproduce from raw + force the specific risk): ~half HELD under a
forced adversary (incl. both mission flagships), ~half were honestly downgraded to MEASURED-B-partial (over-claimed
*evidence strength* — feature-sensitivity, same-modality-≠-decorrelated, discreteness-floor p, pooled-≠-within-subject,
trivial-by-construction — none a fabricated effect), and the identifiability thesis was FIRMED under pre-registration.
The status FIELD is stale on ~36 executed cells; the `verify` TEXT is the record. See §7 + §8.

## 0. THE PREMISE — measured, not asserted
Two falsifiable self-tests, each with the decision-relevant metric (error/discrimination + CI, never just a p-value):
- **Identifiability value — CONFIRMED (firmed under pre-registration)** (`MSK-FAT-FUSION-IDENTIFIABILITY`): fusing two
  decorrelated OUTER modalities (surface geometry + bioimpedance) resolves regional fat **6.3–14.4%** better (held-out
  RMSE) than the best single modality — 6/6 region×cycle CIs exclude 0 under a *pre-registered* 50kHz-BIS spec,
  seed-robust; common-mode falsified (BIS reactance is 75–77% anthropometry-orthogonal). (The informal "8–18%" was
  feature-sensitive; 6.3–14.4% is the pre-registered number.)
- **Prediction value — MARGINAL** (`ORG-BIOLOGICAL-AGE-MORTALITY` discrimination test): a 9-biomarker composite does
  NOT beat the single best lab (albumin) at out-predicting a well-predicted outcome (mortality) on C-index
  (ΔC=+0.0024, CI includes 0), though it is overwhelmingly *significant* (HR p=6.6e-58).
→ **The twin bets on IDENTIFIABILITY** — inferring what no single channel resolves (internal force, fat distribution,
neural drive) — not on out-predicting strong single biomarkers. See
`bt_memory/significance-is-not-discrimination-and-over-determination-buys-identifiability-not-necessarily-better-prediction.md`.

## 1. IDENTIFIABILITY SPINE — inferring the unmeasurable interior from the exterior (the mission core)
- **★ Internal knee contact force from motion** — `KNEE-CELL` (MEASURED-B, flagship): OpenSim-JAM COMAK on
  Grand-Challenge gait, predicted 2.503 BW vs measured 2.584 BW = ~3% peak err (RMS 0.566, r=0.802).
  ⚠️ **READ THAT NUMBER WITH ITS FIVE CAVEATS (2026-07-24, extended 2026-07-26 — the child node carries them, this line used to omit them):**
  (1) **BEST-OF-3, not typical** — og1 3.04%/r0.802 > og4 3.76%/r0.741 > og3 10.35%/r0.658; the honest range is
  **3–10% peak / r 0.66–0.82**, single subject, single limb, no cross-subject replication.
  (2) **n_eff=2, not 4** — the 4 "decorrelated legs" are a design, not an over-determination (`HOLE-KNEE-GC-NEFF2-STRUCTURAL-CAP`).
  (3) ★ **VOID-FLOOR: r and peak-error are WEAK DISCRIMINATORS** (`VOID-FLOOR-KNEE-GENERIC-WAVEFORM-CLEARS-BAR`) — a generic
  8-subject OrthoLoad ensemble waveform with ZERO subject-specific content, given only the measured stance-window timing,
  **matches or beats** this pipeline on r, RMS *and* peak error in all 3 trials. What the model DOES demonstrably own: it
  beats a cross-subject-scrambled-activation null decisively (RMS 5–7× worse, peak err +149/+208%, r goes negative).
  **Quote that margin as the evidence of mechanism — not r, not 3%.**
  (4) ★★ **THE DECISIVE STRENGTH-CONVENTION FALSIFIER IS NOW COMPLETE ON ALL THREE TRIALS** (2026-07-25). Real COMAK
  re-solves with every Fmax divided by 2.276 (toward the measured 268 kN/m² specific tension):

  | trial | peak err base→corr | max activation (ceiling 1.0) | margin | frames converged |
  |---|---|---|---|---|
  | **og1** (this headline's own trial) | **−2.79% → −5.00%** | 0.9941 | 0.59 pp | 82/83 (residual 50.0) |
  | og3 | −10.03% → −12.07% | 0.9505 | 4.95 pp | 91/92 (residual 1.78) |
  | og4 | −3.43% → −7.62% | 0.9992 | 0.076 pp | 88/88 |

  **The correction is absorbed everywhere — nothing saturates, nothing crashes — but never cleanly:** peak error
  degrades in all three, and og1/og3 develop a non-converged frame that did NOT exist at baseline (all three baselines
  converge 100%). ★ **And risk is multi-dimensional: og3 is worst by peak error, og4 by activation margin, og1 by
  convergence residual — no ranking predicts the others.** The 6 untested trials cannot be assumed safe.
  Note the distinction: −5.00% may be judged acceptable; the claim is about MARGIN and ROBUSTNESS, not that the number
  itself is unacceptable.
  (5) ★★ **DOMAIN-DEFAULT AUDIT REVERSES THE HEADLINE ON MEAN ERROR (2026-07-26 — this doc, dated Jul25, predates it)**:
  recomputed over all 3 available trials, OrthoLoad's independent population-median baseline (zero fit to this
  subject) beats the model on mean absolute error — domain-default **3.06%** vs model **5.72%** — winning **2 of
  3 held-out trials (og1, og3)**; the model wins only og4, by a shrunken 2.1–2.3pp margin
  (`AUDIT-VOIDFLOOR-DOMAINDEFAULT-KNEE-SARCOMERE-2026-07-26`, `AUDIT-KNEE-BODYWEIGHT-WINDOW-CORRECTION-OG1-SELECTION-2026-07-26`).
  Same caveat as (1)/(3): n=3 trials/one subject, and OrthoLoad's own between-subject spread (sd≈44 %BW) keeps
  this an *existence* result, not a significance-tested one. This reverses the "single strongest cert in the
  repo" framing that `MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md`, `MECHANISM_UNIFIED_MODEL.md` and
  `MECHANISM_LAYER_COMPLETENESS.md` had carried unqualified since 07-21 (now corrected alongside this line) —
  read the one-line summary of this cert as "beats a scrambled-activation null by 5–7×, loses to a zero-fit
  population median on 2 of 3 trials," not as a bare "~3% error."
  The mission's first predicted-vs-measured internal-force cert. Calibration to the ~0.33 BW bar characterized as needing an MVC
  trial (naive EMG-informing proven to WORSEN it — recipe in `bt_memory/knee-comak-flagship-close-and-fork-build.md`).
- **Foot→hip force-transfer chain** — `MSK-KINETIC-CHAIN-RUNNING` (MEASURED-B): multi-stance n=223 inverse dynamics;
  ankle in-band, knee a characterized systematic ~0.7 Nm/kg low bias (correctable model property), hip heterogeneous.
- **Fat-vs-muscle under the surface** — cornered by 4 decorrelated legs (CT geometry `MSKEL-GEOM-OVERDET`, optical
  appearance `MSK-SURFACE-BODYFAT-DEGENBREAK`, NIRS `MSK-PERFUSION-NIRS`, bioimpedance `MSK-BIOIMPEDANCE-BODYWATER`),
  fusion gain confirmed above. Fat is the shared bottleneck layer across geometry/appearance/perfusion/composition.
- **Neural drive** — `MSK-NEURALDRIVE-MOTORUNIT` (both sub-claims: rate-coding via 2 decomposition pipelines +
  size-principle LIFO on voluntary contractions) + `MSK-NERVE-CONDUCTION` (raw NCS, upper>lower, age-slowing).
- **Molecular→functional** — `MOL-SAMESAMPLE-PATCHSEQ-CALIBRATION`.
  ⚠️ **RETRACTED HEADLINE (banner added 2026-07-25, `HOLE-DOCVSCHILD-CAPSTONE-RESTATEMENT-STALENESS-CENSUS`):** this
  line previously quoted a held-out R²=0.72–0.74 with no hedge. **The node's own later re-verify BREAKS that number as a
  cell-type-identity confound**. Verified against the node directly: its verdict reads *"WEAKENED/BROKEN-headline …
  the R²=0.72–0.74 same-cell transcriptome→ephys COUPLING is BROKEN"* — a leakage-free within-type test (train-fold-only
  cell-type residualisation) collapses it to **R² = 0.014 / 0.041**, i.e. **94–98% of the signal was cell-type
  separation, not within-type prediction**. Do not quote 0.72–0.74 as a coupling result.

## 2. THE VIDEO→CERTIFIED-FORCE CHAIN — ★ CLOSED END-TO-END on real walking (2026-07-19, PARTIAL)
LabValidation (OpenCap, 10 subjects, calibrated multi-view + mocap ground truth + force-plate GRF) landed and the whole
chain ran on real walking video, ground-truth-validated at every link:
- **Link-1 (video→kinematics)** `VIDEO-LINK1-KINEMATICS` (MEASURED-B): ~4.7° lower-limb RMSE / knee 4.2–5.3° over 540
  comparisons (matches the ~4.5° anchor). 2 cameras ≈ 5 (accuracy is pose-pipeline-limited, not camera-count-limited).
- **Link-2 (kinematics→knee joint-reaction force, SO→JR)** — pipeline verified: reproduced mocap-GT JR to +0.09%.
- **End-to-end** `VIDEO-TO-KNEEFORCE-ENDTOEND` (MEASURED-B-partial): video-driven knee JR within **mean 10.75% peak**
  of mocap-GT-driven (median 10.83%, range 0.1–22.2%, 83% <15%, **n=9**).
  ⚠️ The older **8.78% / r=0.94 at n=3** figure that stood here until 2026-07-24 was **superseded by the child node itself**,
  which states the n=3 value "sat on the favorable tail; 10.75% is the honest powered number." Do not re-quote 8.78%. ~4.7° kinematic error propagates *nonlinearly*
  + subject/side-dependently through SO→JR (no fixed multiplier).
- **REMAINING for full autonomy:** the GRF above is the MEASURED force plate. A phone-only twin needs *video-estimated*
  GRF (`MSK-COM-GRF-VIDEO-ONLY`, ~8–17% RMSE, stacks on top) — being chained now; plus running my own Pose2Sim on the raw
  `.avi` (vs the provided OpenCap video-IK). The core loop (action video → certified internal load) is DEMONSTRATED.

## 3. HEALTH-STATE / ORGAN CELLS (prediction & clustering; anchor-passed, adversary-forced)
Cardiometabolic substrate: `ORG-CARDIAC-MECHANICS` (EF/HFrEF; EF-preserved-disease shown), `ORG-GLYCEMIC-METABOLIC`
(HbA1c/glucose 39% discordant, age-drift), `ORG-LIPID-CARDIOMETABOLIC` (Friedewald-derived LDL flagged),
`ORG-METABOLIC-SYNDROME-CLUSTER` (real latent cluster beyond adiposity, 3 methods), `ORG-BIOLOGICAL-AGE-MORTALITY`
(gap→mortality HR 1.02/yr, non-circular, hardest anchor). Also `ORG-RESPIRATORY-SPIROMETRY`, `MSK-BONE-DENSITY-AGE`,
`ORG-HRV-AUTONOMIC`, `MSK-CONNECTIVE` (tendon elasticity), `MSK-FLEXIBILITY-PASSIVE`, `MSK-MUSCLE-ARCHITECTURE`,
`ORG-BARRIER-BBBGUT`, `ORG-TRANSPLANT-GRADIENT`. Energy 1st-law `NODE-0-ENERGY` (~2.3-of-3 terms measured).
PROVEN: `ORG-WAVEFORM-ANCHOR`.
⚠️ **`MOL-GENOME-VARIANT-GROUNDED-ANCHOR` REMOVED FROM THIS LINE (2026-07-25):** it was listed as PROVEN here, but the
node's own `verify` text records an explicit **DOWNGRADE PROVEN → MEASURED-B**. The doc had not been updated. Cite it as
MEASURED-B.

## 4. HONEST NEGATIVES (correctly refuted — a PASS, do NOT rescue)
- `MSK-EMG-FORCE-XSUBJ` — calibration-free cross-subject EMG→force does NOT generalize (shared-protocol-timing leak).
- `MOL-LIFESTAGE-PUBERTY-MENOPAUSE` — power-limited measured-null.
- `MEMBRANE-NERNST-CELL` (OPEN) — no same-well-paired TEER+Papp disruption data exists open; the anchor also
  decouples by permeation sub-state (pore vs leak). Correctly not rescued, not re-refuted.
- Genuine refutations (MITO-CELL same-sample disagreement, FACIAL confirmed on 2nd dataset, KIDNEY/SOD1 measured-null):
  keep refuted — that is the correct science.

## 5. OPERATOR-ACTIONABLES (human owns identity; tool cannot/should not self-pull)
1. **LabValidation_withVideos** (simtk.org/projects/opencap, file 6689 group 2385, 19GB, your SimTK login) →
   closes Link-1 accuracy + the first video→internal-force run.
2. **Dryad 7947gd0** (150 athletes vs 194 HCM, CMR) — AWS-WAF bot-gated; a browser/Playwright pull → athlete's-heart
   cardiac-adaptation contrast.
3. An **MVC/reference trial** for the Grand-Challenge subject → the calibrated (0.33 BW) flagship via the recorded
   Lloyd-Besier EMG-to-activation spec.
Incoming athlete action-videos + known-load lifting clips + training logs = operator cartography of the highest-value
data the tool cannot reach.

## 6. METHOD LEARNINGS (booked in bt_memory — the disciplines that kept this honest)
Aging effects are often composition/load confounds (force the confound; bone-density is the counter-case where it
survives) · a conservation-law closure is a tautology if all terms share one instrument · open-access is a file-API
property not a license badge · an imputed/derived reference can encode your predictor · two-pipeline artifact control ·
tail-leveraged & threshold-blind statistics · significance ≠ discrimination · decorrelated legs must measure the SAME
hidden state · citation-graph data-hunt vector · recompute-to-prove-raw. Full index: `bt_memory/MEMORY.md`.

## 7. CELLS ADDED SINCE FIRST DRAFT (2026-07-18/19 — 18 cells, mostly athlete-performance + integrative)
**Athlete performance (the mission-central layer for incoming action-sequence data), all byte-verified + adversary-forced:**
`MSK-JUMP-POWER-SSC` (CMJ>SJ elastic energy return +12.9%, arm-swing ruled out externally) · `MSK-SPRINT-FV-PROFILE`
(Samozino F0/V0/Pmax, model-fit recovery in trained bands) · `MSK-RUNNING-IMPACT-LOADING` (VALR↑ with speed p=1e-11,
rearfoot>forefoot; MD5 caught a data defect) · `MSK-RUNNING-ECONOMY` (cost-of-transport ~flat, corrected; economy
spread 20-30% over-determined) · `MSK-MUSCLE-FATIGUE-EMG` (MDF-down/RMS-up over-determination, window-invariant) ·
`MSK-KNOWN-LOAD-LIFTING` (pipeline-prep: open data exists but Visual3D-container-blocked; operator-actionable).

**Identifiability — confirmed AND bounded (the mission's core test):** `MSK-FAT-FUSION-IDENTIFIABILITY` (fusing
decorrelated outer modalities cuts regional-fat error **6.3–14.4%** below best single, pre-registered 50kHz-BIS
spec — CONFIRMED, replicated; §0 above — the older informal "8-18%" is superseded, feature-choice-sensitive) ·
`MSK-VISCERAL-FAT-IDENTIFIABILITY` (the exterior corners VAT only R²≈0.68, ~32% irreducible — the identifiability
CEILING for a surface-invisible state; fusion untestable in NHANES, honest).

**Organ / integrative / medical:** `ORG-CARDIAC-MECHANICS` · `ORG-RESPIRATORY-SPIROMETRY` · `ORG-GLYCEMIC-METABOLIC` ·
`ORG-LIPID-CARDIOMETABOLIC` · `ORG-METABOLIC-SYNDROME-CLUSTER` (real latent cluster beyond adiposity) ·
`ORG-BIOLOGICAL-AGE-MORTALITY` (PhenoAge gap → all-cause mortality, hardest anchor) ·
`ORG-FUNCTIONAL-MORTALITY-DISCRIMINATION` (grip adds discrimination where labs don't — modality-decorrelation) ·
`MSK-BONE-DENSITY-AGE` · `MSK-SARCOPENIA-COMPOSITE` (dynapenia dissociation) · `ORG-POSTURAL-BALANCE` (fall-risk).

**The prediction-vs-identifiability premise is now fully measured (§0 refined):** same-modality stacking → ~null
discrimination; decorrelated MODALITY → small-but-real (+0.5 C-index, CI excludes 0); IDENTIFIABILITY → large
(6.3-14.4% error reduction, pre-registered spec — see §0). The twin bets on identifiability. Booked:
`significance-is-not-discrimination-and-over-determination-buys-identifiability-not-necessarily-better-prediction.md`.

**Backlog re-designs resolved (§4):** `AGING-ETA-DRIFT-CELL` → certified POWERED MEASURED-NEGATIVE (η uncoupled from
epigenetic age, 92.7% power) · `ORG-ENDOCRINE-AXES` → anchor-block resolved, concordance method-robust but n-limited ·
`MEMBRANE-NERNST-CELL` → honest OPEN (no same-well data + a method caveat: decorrelated legs must measure the SAME state).

New method learnings booked this stretch: decorrelated-legs-must-measure-the-same-hidden-state · significance≠discrimination
· imputed-reference-tautology · aging=composition/load-confound (+ bone-density counter-case) · parseability axis ·
low-df-R²-doesnt-validate-form · citation-graph data-hunt vector. Full index: `bt_memory/MEMORY.md`.
