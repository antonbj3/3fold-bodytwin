# MECHANISM REPO-WIDE TRUST LEDGER — every headline claim across the twin, machine-cross-checked (2026-07-21)

**Purpose.** `docs/MECHANISM_JOINT_FORCE_SCORECARD.md` did this for the 6-joint family. The operator
cannot personally judge biomechanics and needs to see, AT A GLANCE and without domain expertise, which
of this project's headline claims are backed by real machine-checked evidence vs. prose-only assertion
vs. no external anchor at all. This document extends that scorecard's discipline to **every** doc in
the family: **66 `docs/MECHANISM_*.md` files scanned as of this session, 72 distinct headline-claim rows**
(a handful of docs make two genuinely independent claims; 13 are process/registry docs with no empirical
claim to score). The doc count moved from 58 to 66 mid-session — this is a **live, shared repo** (two
mechanism instances write into it concurrently per `COORDINATOR.md` §0) — the 8 that appeared while this ledger
was being built are included, verified to the same standard as everything else (including one,
`MECHANISM_ANATOMICAL_HAND.md`, that appeared after this ledger's first draft was already written — added
in a follow-up pass, §3).

**Method (what "verified" means here, so the tiers are trustworthy and not just relabeled prose).** For
every doc below, the headline number was re-derived directly from on-disk JSON/CSV/JSONL/`.sto` evidence
files via `python3`/`grep`/direct re-execution this session — never copied from doc prose, never eyeballed
from a table. ~40 of the 71 rows were verified directly by the process that produced this ledger; the
remaining ~31 were extracted by two parallel research passes, each of which was then itself spot-checked
by re-running 4 of its most decisive/surprising findings independently (§10) — all 4 held up exactly,
including one case where a sub-process had to **live re-execute a script** because no evidence file for
that claim existed on disk at all (`MECHANISM_CLIMBING_SCENE.md`, confirmed by reproducing it a second
time, independently, here). Two of my own initial evidence-file guesses to the research passes were
themselves **wrong** (stale checkpoint files) — caught and corrected, not silently trusted; this is
reported explicitly in §10, not hidden, because the operator asked for exactly this kind of finding.

Isolation respected throughout: bodytwin repo only, read-only (no code/model edits), no git commit, no
git push. This document is the only write.

## Tier legend (pick ONE per claim; do not inflate — a rigorous-sounding INTERNAL process with no external
anchor is still "method-only," never upgraded just because the process itself was careful)

- **in-vivo-anchored** — compared against a real in-vivo instrumented/telemetry measurement in a living
  human (OrthoLoad implant telemetry is the source for every row in this tier).
- **cadaveric-or-published-plausibility** — compared against a real, genuinely external/decorrelated
  published dataset or literature figure (cadaveric mechanics, published EMG-timing patterns, published
  metabolic-cost studies, published endurance-time curves, published cartilage-pressure figures, an
  external world-record) but NOT a direct in-vivo force measurement of this twin's own quantity.
- **method-only-no-external-anchor** — the internal process may be rigorous (real citations, self-
  consistency checks, toy-test proofs, pre-registered gates) but there is no external dataset the number
  is actually compared against. **This is the correct tier even when a real PMID is cited for a donor
  model's geometry** — a provenance citation is not a functional-accuracy anchor. This is the single
  largest tier in this repo (31/71 rows, 44%) — expected for an ambitious from-scratch anatomical build;
  not itself a bad sign, but not to be read as more validated than it is.
- **diagnosed-gap** — the doc itself reports a structural limitation, honest negative, or generalization
  failure. This is a PASS in this project's own discipline (an honest negative is not a failure to hide),
  not a lesser tier of embarrassment — but it means there is no positive claim to trust yet.
- **N/A-process** — the doc is a process/mandate/onboarding/design-registry document with no empirical
  headline claim to score.

## Cross-cutting caveat that applies to MANY rows below (read once, not repeated everywhere)

Many docs in §2 explicitly reuse the SAME subject2/`walking1` Static-Optimization solve ("zero re-solve",
stated honestly in each). Good practice — but it means these docs are **not mutually decorrelated
confirmations of each other**; they share a common-mode dependency. If that one SO solve has a systematic
bias, it propagates identically into every "independent-sounding" downstream number (EMG timing, tendon
energy, metabolic cost, fatigue, cartilage kinematics). Agreement between them is internal consistency,
not external validation. Only a row's own named EXTERNAL anchor (OrthoLoad, Koelewijn, Lichtwark & Wilson,
Frey Law & Avin, Arokoski, a world record, published cartilage-pressure figures) actually tests whether
the shared SO solve is right, not just self-consistent.

The joint-force family (§1) is by far the most densely cross-checked cluster in this repo — the same
headline numbers independently reproduce across 3-4 structurally different computation paths (self-computed
subtraction, `opensim.JointReaction`, the sign-bug audit, the diff-reconcile re-derivation) to <0.1%
agreement. Most of the rest of the repo (anatomical grafts, video-corpus scenes, foot-fidelity prototypes)
has not been cross-checked to that density — visible directly in the TIER column below, not asserted here.

---

## §1. Joint-force validation family (the most rigorously anchored cluster — 10 docs, 10 rows)

| DOC | CLAIM | EXTERNAL ANCHOR | VERIFICATION STATUS | TIER | FALSIFIER |
|---|---|---|---|---|---|
| JOINT_FORCE_VALIDATION.md (knee) | Muscle-driven twin over-predicts in-vivo knee contact force 1.515x (391.10%BW self-computed / 391.11%BW JointReaction vs OrthoLoad median 258.22%BW, n=72 trials/9 subj) | OrthoLoad knee instrumented-implant telemetry, in-vivo, n=72 | cross-checked-vs-live-JSON (3 independent files agree to full float precision) | in-vivo-anchored | Self-computed vs JointReaction disagreement growing beyond ~0.1% (currently 0.0039%) |
| HIP_FORCE.md | Twin over-predicts in-vivo hip force 1.412x (386.77%BW self-computed / 387.04%BW JointReaction vs OrthoLoad median 273.93%BW, n=162/18subj) | OrthoLoad hip instrumented-implant telemetry, in-vivo, n=162 | cross-checked-vs-live-JSON | in-vivo-anchored | Self-computed vs JointReaction disagreement growing beyond ~0.1% (currently 0.070%) |
| ANKLE_FORCE.md | Twin lands almost exactly on anchor, 1.016x (484.49/484.48%BW vs cadaveric/FE literature band 390-540, mid 476.7) | 3 cadaveric/quasi-static/FE studies 1977-2000 — no in-vivo ankle program exists anywhere | cross-checked-vs-live-JSON | cadaveric-or-published-plausibility | Any in-vivo ankle telemetry program surfacing and disagreeing with the 445-540%BW band |
| SPINE_FORCE.md | Twin over-predicts in-vivo (partial-load-sharing) lumbar compression 1.74x(no_box)-2.28x(with_box) (296.59/389.17%BW vs OrthoLoad VBR median 170.94%BW, n=13/4subj) | OrthoLoad spine VBR implant telemetry — BUT VBR+posterior-fixation shunts load away, measures a fraction of true compression | cross-checked-vs-live-JSON | in-vivo-anchored (anchor structurally partial — flagged) | Quantifying the VBR load-sharing fraction directly (not yet done) could shrink/grow the ratio |
| SHOULDER_FORCE.md (native) | 0/80 native muscles cross the shoulder at any tested pose -> static optimization structurally impossible; pure-kinematic estimate 5.306%BW vs OrthoLoad median 72.2%BW | OrthoLoad shoulder telemetry exists, but nothing on the twin side to compare | cross-checked-vs-live-JSON | diagnosed-gap | Any native-model muscle actually crossing the glenohumeral joint in a future revision |
| ARM_MUSCLES.md | Donor-grafting real shoulder muscles produces a genuine 3.59x increase over the zero-muscle baseline (5.306%BW -> ~19%BW) | Donor model Seth/Dong/Matias/Delp 2019, PMID 31780916 — geometry provenance, NOT a functional-accuracy anchor | cross-checked-vs-live-JSON; doc itself caught+fixed its OWN earlier wrong PMID (26816372) via live NCBI lookup before publishing 31780916 | method-only-no-external-anchor | A held-out cadaveric shoulder moment-arm dataset comparison (not run) |
| SHOULDER_MUSCLES_CORRECTION.md | Sign-bug-corrected grafted shoulder: 19.06%BW self-computed vs 19.36%BW JointReaction (1.54% apart), ratio 0.264x vs anchor — closes 20.6% of the pre-graft gap | OrthoLoad shoulder median 72.2%BW, in-vivo | cross-checked-vs-live-JSON | in-vivo-anchored (anchor strong; twin-side explicitly fidelity-limited — missing scapulothoracic rhythm+wraps, only 47/70=67% prime-mover sign checks pass) | Adding scapula/clavicle+wraps and re-running; if the ratio doesn't move toward 1x, missing anatomy wasn't the limiting factor |
| ELBOW_FORCE.md | 0/80 muscles cross the elbow at any of 3 poses; 10 N*m ideal actuator is 1.81x under capacity at 5kg load; NO in-vivo elbow program exists anywhere in the literature | NONE — both twin estimate and anchor absent | cross-checked-vs-live-JSON | diagnosed-gap | Any elbow muscle crossing the joint in a future revision, or an in-vivo elbow telemetry program appearing |
| CROSS_ACTIVITY.md | Extends OrthoLoad anchor beyond walking: knee/hip stairs-up 329.6/300.2%BW, stairs-down 348.1/309.4%BW, STS 268.0/184.2%BW, squat 249.8/202.0%BW (all >= walking's 255-302%BW); twin's OWN squat/STS runs over-predict FAR worse (5.04x/3.93x knee/hip) than walking's 1.4-1.5x | OrthoLoad knee/hip stairs+STS+squat buckets, in-vivo, n=9 subj across activities | cross-checked-vs-live-JSON | in-vivo-anchored — flags the over-prediction is activity-DEPENDENT, not constant | A corrected model bringing squat/STS ratio down toward walking's 1.4-1.5x would change "activity-dependent" to "uniform" |
| JOINT_FORCE_SCORECARD.md (synthesis) | At knee+hip, twin over-predicts in-vivo ~1.4-1.5x, a property of the two-step SO-then-subtract method architecture (Moissenet 2014, PMID 24210475), not noise; does NOT generalize uniformly across all 6 joints | Moissenet et al. 2014, PMID 24210475 (quoted verbatim, live-fetched) + all 6 joints' own anchors | cross-checked-vs-live-JSON — independently re-verified against all 8 source JSONs it cites, every number matched to full float precision incl. its own §6 QC finding (identity residuals 3737.5N/4338.0N) | in-vivo-anchored (for the core knee/hip claim; correctly downgrades ankle/shoulder/elbow per-joint) | Any of the 8 source JSONs changing on script re-run |

## §2. Method/mechanism validation layer — does the pipeline itself hold up? (14 docs, 14 rows)

| DOC | CLAIM | EXTERNAL ANCHOR | VERIFICATION STATUS | TIER | FALSIFIER |
|---|---|---|---|---|---|
| STATIC_OPT.md | Original knee cert; corrected headline 391.10%BW (was 233.20%BW pre-sign-fix) | OrthoLoad knee (same as §1 row 1) | **PARTIAL** — the doc's OWN cited `static_opt_knee_results.json` still holds the STALE pre-fix 233.20%BW on disk; 391.10 verified only via SIBLING files | in-vivo-anchored (via siblings) | A naive checker reading only that one file gets the wrong (superseded) number — permanent flag |
| SIGN_BUG_AUDIT.md | Confirmed real sign bug (index ternary backwards both branches), masked by vector-norm output; buggy vs corrected peak vectors point 157.7 deg(knee)/162.6 deg(hip) apart | NONE — internal code-correctness diagnostic | cross-checked-vs-live-JSON (identity residuals 3737.53N knee / 4338.04N hip, independently reproduced) | method-only-no-external-anchor | Re-running post-fix and getting pass=True would falsify the doc's own prediction it must read False (it does) |
| SIGN_BUG_REMEDIATION.md | Bug fixed at source; a SEPARATE copy of the same bug independently found+fixed in the shoulder-with-muscles script | NONE — code-correctness | cross-checked-vs-live-JSON | method-only-no-external-anchor | Any other unaudited caller of the old buggy pattern surfacing elsewhere |
| DIFF_RECONCILE.md | The ~1.68-1.69x self-computed-vs-JointReaction gap was the sign bug, not a diff-scheme difference; collapses to 0.001-0.07% at knee/hip/ankle, re-derived via a 3rd route | NONE of its own (JointReaction is a strong INTERNAL cross-validation — same model/software, different code path — not an external empirical anchor) | cross-checked-vs-live-JSON | method-only-no-external-anchor | Its own secondary/exploratory differentiation apparatus FAILED its own sanity gate (accel_sanity max 7882 vs 2000 ceiling) — disclosed, doesn't affect the primary finding |
| RRA.md | HONEST NEGATIVE: Residual Reduction Algorithm did NOT shrink the pelvis residual — force peak 175.32N (~unchanged) and moment residual got WORSE (83.68 vs pre-RRA 45.11 N*m); both post-RRA gates FAIL | NONE — internal engineering QA gate | cross-checked-vs-live-JSON | diagnosed-gap (honest negative) | RRA succeeding on a different trial/subject or tuning |
| ANKLE_RESERVE_FIX.md | Boosting ankle/subtalar/mtp reserve optimal_force 2.5->150-300N*m PARTIALLY fixes RRA: force-residual FAIL->PASS (0.0508->0.0491) but moment-residual stays FAIL and gets WORSE (0.0868->0.0957) | Hicks et al. 2015 residual-magnitude QC convention — a published methodological threshold, not a force/motion measurement | cross-checked-vs-live-JSON (all ratios reproduced to full precision) | method-only-no-external-anchor | A boost level clearing BOTH gates simultaneously would change "partial" to "full" |
| EMG_TIMING.md | Pre-registered gate: 0 PASS / 3 PARTIAL / 5 FAIL of 8 muscle groups vs published EMG timing, null-calibrated (circular-shift). Mean precision 0.809, recall 0.631, Jaccard 0.514, null-rank 77.0% | Published surface-EMG gait-timing literature (Perry/Rajagopal 2016), decorrelated from this twin's SO | cross-checked-vs-live-JSON | cadaveric-or-published-plausibility | A muscle group's null-rank falling below the pre-registered floor (none currently do) |
| CONTRACTION_DYNAMICS.md | Forward Millard2012 activation-dynamics ODE confirms both pre-registered directional predictions, 6/6 groups (naive excitation LATER 5.9-13.4ms; deconvolved EARLIER 2.5-14.8ms), round-trip falsifier max 4.4ms vs 5.9-13.4ms effect | NONE — internal forward-vs-SO-assumption directional test | cross-checked-vs-live-JSON | method-only-no-external-anchor | Round-trip residual growing to exceed the claimed onset shift |
| TENDON_ELASTIC.md | Achilles peak stored elastic energy 7.76J during walking = 20% of Lichtwark & Wilson 2005's 38J measured in one-legged HOPPING (different, more elastic-demanding activity) | Lichtwark & Wilson 2005 in-vivo ultrasound hopping study — genuinely external, but cross-activity not same-task | cross-checked-vs-live-JSON | cadaveric-or-published-plausibility | A same-activity (walking) published number disagreeing materially with 7.76J |
| METABOLIC_COST.md | Umberger2010 COT_net 7.683 J/kg/m, Bhargava2004 6.198 at 1.065 m/s; twin reads 1.57-1.95x higher than Koelewijn's published level-walking COT | Koelewijn et al. published metabolic-cost study — different subjects/lab | cross-checked-vs-live-JSON | cadaveric-or-published-plausibility | The ratio narrowing to ~1x would need independent confirmation the gap is real physiology not model artifact |
| MUSCLE_FATIGUE.md | 3-compartment fatigue model validated against a REAL external anchor: Frey Law & Avin 2010's directly-fit joint-specific endurance-time power law, 7 joints x 5 intensities, ratio range 0.767-1.572 (median 1.146), ALL 35 points pass gate | Frey Law & Avin 2010, PMID 20069487 (PMC OA) — genuinely external meta-analysis of published endurance-time data | cross-checked-vs-live-JSON | cadaveric-or-published-plausibility — **one of the strongest-anchored claims in this repo** | Any of the 35 points falling outside the pre-registered gate on re-run |
| CARTILAGE_CONTACT.md | TF peak pressure 2.11-4.78 MPa, PF 1.43-6.13 MPa via elastic-foundation contact; 16/16 machine cross-checks PASS incl. mesh-pair agreement <5% median | Literature cartilage-pressure figures, live-verified via NCBI eutils; doc itself: "anchored against a literature MPa figure, not implant telemetry" | cross-checked-vs-live-JSON | cadaveric-or-published-plausibility | Subject-specific cartilage thickness/modulus producing pressures outside the 1-6 MPa literature band |
| KNEE_LIGAMENTS.md | Blankevoort1991Ligament added for ACL/PCL/MCL/LCL; structural self-consistency PASSES (slack-length, force-law, 0/2982 tension-only violations, MCL/LCL opposite-sign engagement) — but ABSOLUTE force/strain magnitudes explicitly NOT validated (stiffness left unscaled) | NONE executed here (Grand-Challenge in-vivo knee-load dataset name-checked as available elsewhere, not run as a comparison); doc DID catch+fix 2 wrong task-given PMIDs via live NCBI lookup | cross-checked-vs-live-JSON | method-only-no-external-anchor | Scaling stiffness by cross-sectional area and comparing resulting forces against a cadaveric/Grand-Challenge dataset (not yet done) |
| UNIFIED_MODEL.md | Merges 6 fidelity forks into ONE model: 232 muscles/84 ligaments/27 bodies/329 total force elements all match expected counts exactly; re-confirms the 1-DOF-vs-6-DOF knee-ligament mismatch persists post-merge (disclosed, not hidden) | NONE — internal software-composability/integration test | cross-checked-vs-live-JSON (expected==measured counts, exact) | method-only-no-external-anchor | Any component-count mismatch on merge-script re-run |

## §3. Anatomical model-building / muscle grafts (9 docs, 11 rows)

| DOC | CLAIM | EXTERNAL ANCHOR | VERIFICATION STATUS | TIER | FALSIFIER |
|---|---|---|---|---|---|
| ANATOMICAL_HAND.md | First muscle-driven hand dial-turn (index finger only): 2 real forearm flexors (FDP2_r/FDS2_r) + 1 extensor (ED2_r) grafted; moment-arm cross-check (OpenSim `computeMomentArm` vs independent finite-difference) 12/12 PASS (rel. error <=6.3e-5); sign-contrast (extensor opposite flexors) PASS at all 4 joints; pre-registered literature order-of-magnitude moment-arm band 11/12 PASS, 1 disclosed miss (FDP2_r@DIP 10.56mm vs [1,10]mm band, 5.6% over) | Real, live-verified anatomy literature (An 1983 PMID 6619158, Lee 2008 PMID 18387615, Buchholz 1992 PMID 1572336, etc.) for the order-of-magnitude bands — BUT several primary numeric tables (An 1983 especially) are paywalled/could not be retrieved this session, disclosed rather than fabricated; the correct full-hand donor (McFarland 2023, PMID 36301780, 43-muscle) is identified but NOT used (fetch not attempted, disclosed, not silently substituted) | cross-checked-vs-live-JSON (`anatomical_hand_evidence.json`: 12/12 cross_check_pass=True, 11/12 band_pass matching the doc's disclosed single miss, exact) | cadaveric-or-published-plausibility (real literature bands, honestly caveated including a disclosed miss and disclosed paywall gaps) | A finite-difference cross-check failing on any of the 12 pairs, or the literature band miss growing beyond a single joint |
| MUSCLE_AUDIT.md | Are the 80 native muscle-tendon actuators "correctly modeled"? 0/80 flagged outside pre-registered physiological bounds; 16/16 posterior-chain muscles present | NONE — cross-anchor uses 2 OTHER already-scaled `.osim` files from the same lineage, not literature/in-vivo data | cross-checked-vs-live-JSON (`n_flagged_bounds=0`, `n_posterior_chain_present=16/16`) | method-only-no-external-anchor | Any of the 80 muscles found outside pre-registered bounds, or <16/16 posterior-chain present |
| ERECTOR_SPINAE.md (build) | 88/126 candidate lumbar-extensor fascicles kept (38 excluded); 88/88 verified positive (extensor) moment arm at 2 poses | Beaucage-Gauvreau 2019, PMID 30714401 — donor identity only, not a moment-arm-accuracy anchor | cross-checked-vs-live-JSON | method-only-no-external-anchor | Any of the 88 muscles showing negative/zero moment arm at either pose |
| ERECTOR_SPINAE.md (capacity) | Peak extension capacity 228.0 N*m (neutral) / 148.6 N*m (30 deg flexed) vs Arokoski 2004 in-vivo extension torque (147.3-170.1 N*m), ratio ~1.3-1.5x, "physiologically sane" | Arokoski JP et al. 2004, PMID 15129408 (in-vivo dynamometry, n=9) — genuine independent cohort, BUT it is a **chronic low-back-pain cohort, not healthy** (this caveat is disclosed only later, in TRUNK_FLEXORS.md, not in this doc itself) | cross-checked-vs-live-JSON | cadaveric-or-published-plausibility | Capacity falling far outside a plausible multiple of Arokoski's measured torque |
| TRUNK_FLEXORS.md (required force) | Hip-hinge/band pose requires -41.749 N*m at `lumbar_extension`, doc's own words: **"required by gravity alone"** | NONE — internal virtual-work computation | **cross-checked-vs-live-JSON for the single-pose value, BUT with two confirmed issues** (below) | method-only-no-external-anchor | See notes — this row is one of the two most important flags in this ledger (§10) |
| TRUNK_FLEXORS.md (capacity) | Peak flexion capacity 226.0/282.7 N*m vs Arokoski measured flexion torque (72.0-93.5 N*m) -> 2.4-3.9x, explicitly flagged "NOT flattering" | Arokoski JP et al. 2004, PMID 15129408 (same CLBP cohort as above; caveat disclosed explicitly and honestly HERE, unlike the erector-spinae doc) | cross-checked-vs-live-JSON | cadaveric-or-published-plausibility | A healthy-cohort re-comparison (2 candidate papers located, not resolved this session) landing the ratio near 1.3-1.5x would soften but not eliminate the discrepancy |
| SCAPULA_CLAVICLE.md | Scapulohumeral rhythm "measured" 0.489 (mean) vs target 0.500 (classic 2:1); independent-FK cross-check diff 1.06 deg vs 0.5 deg gate | Inman 1944/de Groot & Brand 2001/Ludewig 2009 (classic 2:1 ratio) + reused donor PMID 31780916 | cross-checked-vs-live-JSON — **BUT largely circular**: `RHYTHM_RATIO_SCAPULAR_TO_GH = 0.5` is HARD-CODED as the `CoordinateCouplerConstraint`'s slope in the script itself (independently confirmed by reading `add_scapula_clavicle.py` line 113 and the evidence JSON), so "0.489 vs 0.500" mostly checks whether 3D-composition math reproduces its OWN imposed constant, not real scapular kinematics — doc's own §7.5 admits no real scapular mocap trial exists in this corpus | method-only-no-external-anchor | A real scapular mocap trial (none exists) giving a materially different ratio |
| CROSSBRIDGE_MODEL.md | History-dependence discriminator: RFE = -0.085% (approx 0) — RFE/FDE (the #1 reason to prefer cross-bridge over Hill) provably ABSENT in this uniform-sarcomere 2-state construction; isometric force-length plateau reproduces Gordon-Huxley-Julian 1966 exactly | Huxley 1957 PMID 13485191 + Gordon-Huxley-Julian 1966 PMID 5921536 (frog plateau 2.05-2.2um) — both real, both live-verified | cross-checked-vs-live-JSON — GHJ-1966 reproduction is itself doc-labeled "a pipeline-correctness check, not a novel claim" (solved-to-match); the RFE=0 headline has no independent quantitative anchor, only internal hold-duration convergence | method-only-no-external-anchor | A hold-duration sweep NOT shrinking the RFE residual toward zero |
| MOTOR_UNIT.md | Fuglevand-style recruitment: size-principle order holds under stochastic noise (all pre-registered gates PASS); soleus-specific "anti-onion-skin" firing-rate direction (Spearman -0.965) matches Oya 2009's real finding | Oya 2009 — anchors only the qualitative DIRECTION of a design choice; magnitudes are explicitly generic/illustrative | cross-checked-vs-live-JSON | method-only-no-external-anchor | A real soleus motor-unit decomposition dataset showing classic (not anti-onion-skin) ordering |
| MYOFASCIAL.md | Quantifies epimuscular myofascial force transmission from primary literature: proximal-distal EDL force differences 0-22.7%/0-14% (cat/rat dissection); 40.5+-5.9% plantaris-to-calcaneus (rat); caught+fixed a real citation author-order error live | Huijing/Baan/Maas 2001-2010 series + Sandercock/Maas studies — real PMIDs, live-verified | cross-checked-vs-live-JSON | cadaveric-or-published-plausibility — **but every quantitative study cited is cat/rat ex-vivo dissection**; the one human study (Bojsen-Moller 2010) gives no % | Any cited PMID being mis-quoted on independent re-verification (none found this session) |
| FIDELITY_ARCHITECTURE_AUDIT.md (self-audit synthesis) | Method/discipline CONFIRMED strong; cross-bridge-emerges-Hill claim CONFIRM-NARROW (concentric matches Hill R^2=0.9987, but eccentric over-predicts 2.57x vs Hill's 1.4x cap, and RFE is provably absent); composability verdict = NOT YET COMPOSABLE; 4 substrate ceilings confirmed | Inherits its children's anchors (Moissenet 2014, GHJ 1966, etc.) | cross-checked-vs-live-JSON — independently re-verified both the joint-force numbers AND the crossbridge citation metadata it cites; appropriately self-critical, no inflation found | method-only-no-external-anchor (a synthesis doc — inherits the tier of whichever child claim is weakest for any given sub-point) | Any child doc's number changing on re-run |

## §4. Video-corpus / motion-capture pipeline (13 docs, 16 rows)

| DOC | CLAIM | EXTERNAL ANCHOR | VERIFICATION STATUS | TIER | FALSIFIER |
|---|---|---|---|---|---|
| POSE_PIPELINE.md | 264/264 corpus clips processed; 228/264 (86.4%) pass the full gate suite | NONE | cross-checked-vs-live-JSON (`data/msk_ik/_overall_pass_manifest.json`: n_total=264, n_pass=228, exact) — **note: the doc's own originally-implied log (`data/msk_pose/_batch_log.jsonl`) only has 258/264, the true authoritative file is the manifest** | method-only-no-external-anchor | Any clip missing an OVERALL_PASS key, or a recount != 228 |
| CORPUS_IK.md | On the priority clip: marker RMS 0.0857m PASS; knee/hip_flexion range gates PASS; ankle range gate FAILS | NONE for this monocular-IK claim itself (cites a different doc's mocap-based validation as background context only, not a direct anchor) | cross-checked-vs-live-JSON, exact on every number, raw and smoothed | method-only-no-external-anchor | A re-run producing an ankle range inside [5,60] deg would flip FAIL to PASS |
| CORPUS_IK_AGGREGATE.md | 228/228 clips aggregated; CLEAN(lenient)=181/DEGENERATE=47; CLEAN(strict)=67/DEGENERATE=161 | NONE (internal regression vs an earlier pilot clip from the same project) | cross-checked-vs-live-JSON — **note: the initially-suggested evidence file (`data/msk_ik/_aggregate_summary.json`) is a STALE mid-batch checkpoint frozen at 19/228; the correct, current file is `scripts/msk/corpus_ik_aggregate_evidence.json` (228/228) — independently re-confirmed both facts here** | method-only-no-external-anchor | A nonzero mismatch in the pilot-clip regression self-check (script exits 1 on failure; none occurred) |
| BAND_NOBAND_CONTROL.md | Semimembranosus's near-saturated activation under the band is a POSE effect, not a band effect (detaching the band leaves it pinned at ~1.0) | NONE — internal cross-check against 2 other models from the same pipeline | cross-checked-vs-live-JSON — **one small confirmed rounding slip**: doc states semimem_l no-band activation as "1.00000," raw value is 0.99995988 (rounds to 0.99996) — negligible, but a real, confirmed mismatch | method-only-no-external-anchor | Detaching the band and finding activation drops meaningfully below ~0.999 |
| BAND_STATIC_OPT.md | Gastroc/hamstrings/glutes carry the load; 88 erector-spinae muscles stay at the activation floor (wrong-sign requirement, not a solver bug); moment-balance cross-check approx 0 | NONE — internal virtual-work vs SO+moment-arm cross-check | cross-checked-vs-live-JSON, all 4 headline numbers match to stated precision | method-only-no-external-anchor | Any erector found with a genuinely negative (flexion-capable) moment arm at this pose |
| BAND_POSTERIOR_CHAIN.md | Only the "deadlift feet-to-hands" band config loads the full ankle-to-lumbar chain directly (100% direct at all 9 coordinates); a pelvis-anchored config gives exactly zero coupling | NONE — pure geometric/moment-arm derivation | cross-checked-vs-live-JSON, exact | method-only-no-external-anchor | Any nonzero indirect-GRF component at one of the 9 coordinates for this config |
| MSK_ELASTIC_BAND.md (ID-tool bug) | OpenSim 4.6 `InverseDynamicsTool` returns physically-absurd knee/hip moments (-1065 to +5741 N*m) at a hip-hinge pose; workaround = virtual-work method | NONE — a tool-defect finding, not a biomechanics claim | cross-checked-vs-live-JSON for the hip-hinge-pose numbers (reproduced from raw `.sto`) — **the doc's SEPARATE neutral-pose sub-table (-1252.62 N*m, "1815x too large") has no on-disk artifact anywhere; doc-prose-only for that specific number** | diagnosed-gap | A re-run at the file-default (all-zero) pose with gravity only giving a plausible (<10 N*m) moment would falsify the bug |
| MSK_ELASTIC_BAND.md (propagation table) | Band's direct effect is largest at lumbar (+4.58 N*m); hip/knee/ankle effect is small (0.6-1.7 N*m) and mostly indirect via GRF | Uchida et al. 2016 Thera-Band tension table — anchors the INPUT band-tension calibration only, not the output joint-moment numbers (doc itself calls this a "soft plausibility check," not a tight anchor) | cross-checked-vs-live-JSON, exact to stated rounding | method-only-no-external-anchor | A re-run showing indirect (not direct) coupling at lumbar |
| PLYO_JUMP_CASCADE.md | A real ballistic jump still does not resolve an ankle-knee-hip-back cascade at ~30fps; corpus screen 37->21->6; freefall misread inflates peak 1.29-1.37x | NONE | cross-checked-vs-live-JSON, exact on every figure | diagnosed-gap | Any propulsion-phase peak with a nonzero inter-joint gap surviving all 3 margin values |
| CORPUS_SLOWMO_CASCADE.md | Of 12 db-flagged ballistic candidates, only 2/12 physics-confirm real slow-motion; negative controls cluster near speed-ratio 1; cascade "mostly still tied" across 5 re-timed clips | NONE — internal physics (gravity/projectile) cross-check of the corpus DB's own metadata flag | cross-checked-vs-live-JSON — "4/5 tied" is a Bottom-Line simplification of a messier 3-clean-tie/1-partial(111ms lead)/1-clean-20ms-stagger split; the doc's own detail table discloses this nuance, so it's rounding, not a hidden error | method-only-no-external-anchor | The one confirmed stagger clip failing to replicate under an angular-velocity-peak cross-check (doc already reports it does not replicate) |
| FORCE_SCENES_BATCH.md (Newton method) | Jump clip peaks 23-30% higher than squat/other clip at ankle/knee/hip; squat vs other within 0.4-3.5%; regression check PASS at 0.0000% | NONE for these 3 clips (self-regressed against this project's own earlier pilot clip — not decorrelated) | cross-checked-vs-live-JSON, computed ranges match exactly | method-only-no-external-anchor | A 4th clip reversing the squat approx other / jump-higher pattern |
| FORCE_SCENES_BATCH.md (SO implausibility) | Static-Optimization tier numerically converges but is NOT dynamically plausible on any of 3 clips (179-1044%BW after the sign-bug fix) | NONE | cross-checked-vs-live-JSON, exact extremes | diagnosed-gap | A 4th independent clip landing inside a plausible band (<~300%BW) |
| FORCE_TRANSMISSION_SCENE.md (force magnitudes) | 5-rep squat: ankle 119-178%BW, knee 113-170%BW, hip 86-131%BW, back 116-176%BW; mass-shedding cross-check within ~2-3% of segment mass | Segment mass from the subject's OWN scaled OpenSim model — Newton's-3rd-law self-consistency, NOT an independent cadaver/population anthropometry table | cross-checked-vs-live-JSON, every rep's values fall inside stated ranges, mass-shedding ratios match to 3 decimals | method-only-no-external-anchor | A mass-shedding ratio outside ~0.95-1.05 on a held-out rep |
| FORCE_TRANSMISSION_SCENE.md (no cascade resolves) | All 4 chain levels (ankle/knee/hip/back) peak in the SAME video frame, in all 5 reps — no resolvable sequence at ~30fps | NONE | cross-checked-vs-live-JSON, exact 5/5 | diagnosed-gap | Any single rep with a >1-frame (33ms) gap between ankle and hip peak |
| CLIMBING_SCENE.md | 2/3 solo-run clips cleanly pass every pre-registered gate; speed estimates match the IFSC 15m-wall world record to within ~3% | Wikipedia/IFSC world record (Zhao Yicheng 4.54s, 15m) — genuinely independent, fetched live | **cross-checked by INDEPENDENTLY RE-RUNNING the script live a second time this session** (no persisted gate-table JSON exists anywhere — confirmed absent both times) — v_climb 3.387/3.375 m/s vs WR pace 3.304 m/s = +2.5%/+2.1%, reproduced exactly | cadaveric-or-published-plausibility | Re-running the 3rd (partial-pass) clip and getting all 3 gates to pass, or >5% WR-pace deviation on either full-pass clip |
| YOUTUBE_CASCADE_POC.md | Network fetch (yt-dlp) is bot-blocked in this environment; 23/23 direct-ID fetches failed; zero clips downloaded | NONE | "zero clips downloaded" confirmed directly by filesystem check (no files under `data/youtube/` except an unrelated note); the specific "23/23" count has **no persisted artifact anywhere in this repo** | diagnosed-gap; the 23-count specifically is **unverified** | Finding even one successfully-downloaded video file under `data/youtube/` |

## §5. Foot-fidelity sub-project (4 docs, 4 rows)

| DOC | CLAIM | EXTERNAL ANCHOR | VERIFICATION STATUS | TIER | FALSIFIER |
|---|---|---|---|---|---|
| FOOT_MULTISEGMENT.md | Held-out generalization FAIL: DJ2 max-abs 5.76 deg (`subtalar_angle_l`) vs pre-registered 5.0 deg ceiling (walking1 in-sample 0.25/3.69 deg PASS) | NONE — the 0.5/5.0 deg ceiling is an internal project convention (traced to an earlier smoke-test baseline), not a literature/clinical tolerance. Genuinely pre-registered BEFORE this specific test (confirmed via file mtimes: the ceiling first appears hours earlier in FOOT_FIDELITY_PLAN.md), not tuned post-hoc | cross-checked-vs-live-JSON, exact (`regression_max_abs_diff_deg: 5.7603526`) | diagnosed-gap (doc's own honest generalization-failure finding) | Additional held-out trials consistently <5.0 deg would overturn "does not generalize" |
| FOOT_FIDELITY_PLAN.md | Hallux-split prototype: regression PASS (median 0.057deg/max 1.34deg vs 0.5/5.0deg ceiling); new DOF stays near-null (range 2.18deg, max 7.20deg) exactly as predicted by the 3-marker sensing-ceiling argument | NONE — self-regression + an internally-predicted null result (a sound geometric observability argument, not a tautology gate, but not checked against outside data either) | cross-checked-vs-live-JSON, exact | method-only-no-external-anchor | New DOF showing a LARGE (not near-null) range despite zero informing markers |
| FOOT_DRILL_CORPUS.md | 0/893 corpus captions name a midtarsal-family term — no clip can anchor the failed midtarsal DOF; 76/893 reconciled as trustworthy foot/ankle content (vs naive over-counts of 132/360) | NONE (corpus-mining against the project's own video corpus — a legitimate direct caption-text check, not circular) | cross-checked-vs-live-JSON, exact (28+48=76 STRICT+MODERATE) | diagnosed-gap (doc's central point: the corpus cannot anchor the placeholder axis) | A fresh grep of the 893 captions finding >=1 genuine midtarsal/subtalar/chopart mention |
| FOOT_CLOSEUP_POSE.md | A toe/palm-width-ratio "lead" separates 2 true-barefoot detections (0.070-0.140) from 1 shod false-positive (1.558); raw hand-confidence score does NOT separate them | NONE (n=3 internal frames; doc itself explicitly flags this as unvalidated, "nowhere near the diverse instance-space required") | cross-checked-vs-live-JSON, exact | diagnosed-gap (an honestly-labeled exploratory lead, not a validated result) | Testing more shod/barefoot examples and finding overlapping ratio ranges |

## §6. Inherited components (1 doc, 1 row)

| DOC | CLAIM | EXTERNAL ANCHOR | VERIFICATION STATUS | TIER | FALSIFIER |
|---|---|---|---|---|---|
| HAND_INHERIT.md | CS's inherited hand-retarget RMS 0.0749m (bodytwin climbing clip) vs pre-registered <0.05m band -> FAILS (reproduces CS's own clip's 0.0689m miss); the raw MediaPipe hand-span (0.084m) falls OUTSIDE a real published adult anthropometric band (0.15-0.21m) | Published adult wrist-to-mid-fingertip anthropometric band (0.15-0.21m) — genuinely independent; pipeline's raw span falls outside it, a real disclosed miss, not circular | cross-checked-vs-live-JSON, exact | diagnosed-gap (doc's own bottom line: this component "cannot be the donor" for the muscle-driven-hand goal) | A different real clip reproducing RMS <0.05m would falsify "generalizes as a genuine miss" |

## §7. Data registry / OrthoLoad infrastructure (2 docs, 3 rows)

| DOC | CLAIM | EXTERNAL ANCHOR | VERIFICATION STATUS | TIER | FALSIFIER |
|---|---|---|---|---|---|
| ORTHOLOAD_FETCH.md | All 6 implant IDs fully crawled/downloaded: 14,897 attempted -> 14,893 succeeded (4 genuine 404s); AKF counts per joint sum to 3,942 | NONE — a data-provenance/download-integrity claim, self-verified via HTTP status | cross-checked-vs-live-JSON, exact on every sub-number (independently recomputed from the raw CSV: 14897 rows, 14893 True/4 False) | method-only-no-external-anchor | A 7th distinct implant ID, or a per-joint AKF recount disagreeing with the doc |
| ORTHOLOAD_INDEX.md (parse integrity) | 3,942/3,942 AKF files parsed, 0 parse errors | NONE | cross-checked-vs-live-JSON, exact | method-only-no-external-anchor | Any file producing a parse exception on a from-scratch re-run |
| ORTHOLOAD_INDEX.md (literature sanity check) | Knee level-walking (n=169) median 254.8%BW inside a pre-registered 250-300%BW band; hip walking (n=382) median 302.4%BW inside 210-330%BW band | "Pre-registered bands from the task brief" — **no DOI/PMID actually quoted in the doc**; plausibly traces to the SAME OrthoLoad investigators' (Bergmann/Kutzner) own published summary of this SAME cohort, which would make this same-cohort corroboration, not independent literature | cross-checked-vs-live-JSON, exact | in-vivo-anchored, **flagged**: decorrelation of the sanity band itself is unconfirmed | Tracing the band's source: if it has zero subject/instrumentation overlap with this corpus, decorrelation is confirmed; if it traces to Bergmann/Kutzner's own OrthoLoad analyses, downgrade to same-cohort corroboration |

## §8. Process / meta / registry documents — no empirical headline claim (13 docs, 13 rows)

| DOC | CLAIM | VERIFICATION STATUS | TIER | FALSIFIER |
|---|---|---|---|---|
| CELL_INVENTORY.md | 49-cell wave-1 design catalog (REAL 6/WEAKENED 6/RE-SCOPED 10/REFUTED 7/FENCED 5/un-audited 15); explicitly "none MEASURED yet" per its own header | doc-prose-only (by its own admission) | N/A-process | Any cell's verdict changing once executed |
| HARDENED_CONVENTIONS.md | Process/convention doc (isolation, fold path, node schema) | N/A | N/A-process | N/A |
| LAUNCH_MESSAGE.md | Onboarding prompt block | N/A | N/A-process | N/A |
| MISSION.md | Vision/architecture statement | N/A | N/A-process | N/A |
| MSK_BUILD_CHARTER.md | Mandate/pivot announcement | N/A | N/A-process | N/A |
| MSK_BUILD_PLAN.md | Capability inventory (paths/versions/imports) "checked live" — an infra claim, not biomechanics | unverified this session (did not re-run its specific checks; the shared foundational one — OpenSim version — was verified live separately, see next row) | method-only-no-external-anchor | Re-running its cited commands and getting a different answer |
| MSK_ENV.md | `.venv-msk`: OpenSim `4.6-2026-06-22-85aaf64`, Python 3.12.8, built fully offline | **cross-checked LIVE by me this session** — re-ran `opensim.GetVersion()` directly, got the exact string every downstream doc in this ledger depends on | method-only-no-external-anchor (infra fact, but foundational — everything else assumes this is true) | `opensim.GetVersion()` returning a different string (it did not) |
| PREAMBLE.md | Subagent discipline block | N/A | N/A-process | N/A |
| SCALING_DOCTRINE.md | Forward-looking architecture, explicitly flagged "NOT current implementation" | N/A (self-flagged) | N/A-process | N/A |
| SEED_CANCER_IMMUNE_WARGAME.md | "a partial seed to be expanded (not a spec)" per its own header | N/A (self-flagged) | N/A-process | N/A |
| STARTUP.md | Boot-sequence doc | N/A | N/A-process | N/A |
| STATE.md | "All ~92 designed cells are EXECUTED"; 20 strongest-headline cells adversarially re-verified, "~half HELD... ~half honestly downgraded" | **UNVERIFIED — attempted, could not reconcile.** Live graph query: 998 total nodes, status {OPEN:923, ASSUMED:66, REFUTED:5, DEFERRED:3, PROVEN:1} — does not obviously match "~92/20"; `bt_memory/LEDGER.jsonl` is an event stream (214 lines), not a per-cell tally; `data/mechanism_catalog/cell_designs_wave*.json` only goes up to wave1. Doc dated 2026-07-18, older than every other doc in this ledger — likely stale | method-only-no-external-anchor, **flagged STALE/UNRECONCILED** | A live per-cell tally (by whatever tag distinguishes "cell" from the other 998 nodes) matched against 92 and 20 |
| SUBAGENT_PLAYBOOK.md | Process rules for subagent fan-out | N/A | N/A-process | N/A |

---

## §9. Tier counts (72 rows across 66 docs)

| TIER | count | % |
|---|---:|---:|
| in-vivo-anchored | 8 | 11% |
| cadaveric-or-published-plausibility | 11 | 15% |
| method-only-no-external-anchor | 31 | 43% |
| diagnosed-gap | 12 | 17% |
| N/A-process | 10 | 14% |

**Reading this honestly:** only 11% of headline claims reach a real in-vivo comparison (concentrated almost
entirely in the original joint-force family + its cross-activity extension — the part of this project that
has had the most rounds of adversarial re-verification). 44% — the plurality — are internally rigorous but
have never been compared to anything outside this twin's own pipeline; this is the expected shape for an
ambitious from-scratch anatomical-modeling effort at this stage, not itself a red flag, but it means "passed
its own checks" and "true" are NOT the same claim for most of this repo. 17% are honest, explicitly-disclosed
negatives — a good sign for the project's discipline, not a bad sign for the ledger.

## §10. FLAGGED — highest scrutiny (prose-only sub-claims, unresolved decorrelation questions, and confirmed
mismatches; these are the rows the operator most needs surfaced, not buried)

1. **TRUNK_FLEXORS.md's "-41.749 N*m, required by gravity alone" is confirmed MISLABELED.** Its own
   pipeline's JSON (`band_static_opt_results.json`) shows `gravity_only = -50.100 N*m`; -41.749 N*m is
   `gravity_plus_GRF_plus_band` — i.e. gravity **net of the band's own +8.351 N*m contribution**, not
   gravity alone. Independently re-derived and confirmed directly from the JSON this session (not just
   trusted from the research pass that first caught it).
2. **The "8-point robustness sweep" that both TRUNK_FLEXORS.md and BAND_STATIC_OPT.md cite to call the
   -41.749 N*m figure "robust" is doc-prose-only.** `BAND_STATIC_OPT.md` §5 states specific numbers
   (-15.2 to -75.6 N*m, monotonic, no sign flip, across lumbar_extension in [-20,+30] deg) but no array of
   these 8 values exists in any on-disk JSON, and `band_static_opt.py`'s `gravity_only_decomposition()` is
   called exactly once in `main()` — confirmed directly by reading the script. The single-pose value itself
   IS solidly JSON-verified; this specific robustness qualifier is not.
3. **SCAPULA_CLAVICLE.md's "measured" 2:1 scapulohumeral rhythm is largely circular.** Confirmed directly:
   `RHYTHM_RATIO_SCAPULAR_TO_GH = 0.5` is hard-coded in `add_scapula_clavicle.py` (line 113) and imposed as
   the `CoordinateCouplerConstraint`'s slope — the "0.489 vs 0.500" comparison mostly checks whether the
   3D-composition math correctly reproduces a constant the script itself imposed, not real scapular
   kinematics. The doc's own §7.5 admits no real scapular mocap trial exists in this corpus to check against.
4. **CLIMBING_SCENE.md's headline gate table (the "~3% of world record" claim) is never persisted to
   disk anywhere** — confirmed absent twice, independently, by two different processes this session. Both
   times it was reproduced by live-re-executing `scripts/msk/climbing_scene.py anchor-demo`, which is a
   *stronger* check than reading a static file would have been, but means a future reader who doesn't
   re-run the script has nothing to check the claim against.
5. **STATE.md's "~92 cells executed / 20 re-verified" could not be reconciled against the current live
   graph** (998 nodes, only 1 tagged PROVEN) or `bt_memory/LEDGER.jsonl` (an event log, not a per-cell
   tally). Dated 2026-07-18, oldest doc in this ledger — flagged as likely stale rather than false.
6. **YOUTUBE_CASCADE_POC.md's specific "23/23 fetches failed" count has no on-disk artifact anywhere** in
   this repo. The qualitative claim (zero clips downloaded) is independently confirmed by direct filesystem
   check; the exact count is not.
7. **MSK_ELASTIC_BAND.md's neutral-pose ID-tool-bug sub-table (-1252.62 N*m, "1815x too large") has no
   persisted artifact**; only the separately-reported hip-hinge-pose numbers were reproducible from a raw
   `.sto` file.
8. **ORTHOLOAD_INDEX.md's literature-sanity-check band (250-300%BW knee, 210-330%BW hip) cites no
   DOI/PMID** and may trace to the same OrthoLoad investigators' own summary of this identical cohort —
   if so it is same-cohort corroboration, not independent decorrelated literature, despite currently being
   tiered in-vivo-anchored (flagged, not downgraded, because this was not resolved either way this session).
9. **MSK_BUILD_PLAN.md's "checked live" capability inventory was not independently re-run this session**
   (only its single most load-bearing shared claim — the OpenSim version string every other doc depends
   on — was independently re-verified live, see MSK_ENV.md's row in §8).
10. **STATIC_OPT.md's own cited `static_opt_knee_results.json` still contains the STALE, pre-sign-fix
    233.20%BW knee number on disk** — a naive machine cross-check against only that one file (ignoring the
    doc's correction banner) would report the wrong, superseded headline number. The correct 391.10%BW is
    real and independently triangulated (see §1), but this specific file is a live trap for future
    automated verification.
11. **Minor, immaterial:** BAND_NOBAND_CONTROL.md states semimem_l's no-band activation as "1.00000"; the
    raw value is 0.99995988 (rounds to 0.99996) — a real but negligible rounding slip, included here only
    for completeness since the operator asked for every prose/JSON mismatch to be surfaced.

Every row not listed above and not tiered `diagnosed-gap` or `N/A-process` had its headline number
reproduced exactly (to the stated rounding) from a real, currently-existing on-disk file this session.

## §11. Files / methodology

- This doc: `docs/MECHANISM_TRUST_LEDGER.md` (new, this session). Isolation respected: bodytwin only,
  read-only against all model/data files, no git commit, no git push.
- Extends `docs/MECHANISM_JOINT_FORCE_SCORECARD.md` (§1 above independently re-derives and confirms its
  numbers rather than re-quoting them) to the full `docs/MECHANISM_*.md` family (65 files as of this
  session — re-scanned twice mid-session to catch concurrently-added docs; see header).
  `docs/MECHANISM_CELL_INVENTORY.md` and `data/MECHANISM_ANCHOR_GRAPH.json` are a **different**, broader
  whole-body design catalog (not MSK-specific) — cross-referenced only where a doc in this family
  (`STATE.md`) makes a claim about it.
  `docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md` is this project's own prior self-audit at the
  architecture level; this ledger operates at the individual-claim level and is compatible with, not
  contradictory to, that audit's findings (independently spot-checked, §3).
- ~41 of 72 rows were verified directly this session against on-disk JSON/CSV/JSONL/`.sto` files via
  `python3`/`grep`. The remaining ~31 were extracted by two parallel passes, each independently spot-checked
  afterward against their most decisive/surprising findings (4 targeted re-verifications, all 4 confirmed
  exactly, including one live script re-execution reproduced a second time independently — see §10).
- No scripts under `scripts/msk/` were edited; `scripts/msk/climbing_scene.py anchor-demo` was executed
  read-only (reads existing pose JSON, prints to stdout, writes nothing) purely to verify a claim that has
  no other on-disk record.
