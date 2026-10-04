# MECHANISM CELL INVENTORY — the bio build-goal map (wave-1)

Synthesized from `data/mechanism_catalog/cell_designs_wave1.json` → `data/MECHANISM_ANCHOR_GRAPH.json`. **49 cells, adversarially audited (none MEASURED yet — verdicts are on the DESIGN).**

## Certifiability (design-audit verdicts)

| bucket | count | meaning |
|---|---|---|
| REAL | 6 | design survives adversarial verify + external fact-check → ready to EXECUTE |
| WEAKENED | 6 | fixable defect, residual noted per cell |
| RE-SCOPED | 10 | honest weaker claim (archive lacks the ideal instrument) |
| REFUTED | 7 | anchor/leg defect not yet resolved |
| FENCED / HONEST-NEG | 5 | no valid anchor in archive (=PASS, non-promotable) |
| un-audited | 15 | (honest-neg no-anchor cells / node-0) |

**REAL (execute these first):** PARALYSIS-STROKE-STAGE-LOCALIZATION, ORG-WAVEFORM-ANCHOR, ORG-LUNG-GASEXCHANGE, ORG-ENDOCRINE-AXES, ORG-DIGESTIVE-ENZYME, MOL-SAMESAMPLE-MULTIOMIC-CALIBRATION

Each cell = an occluded hidden state cornered by ≥2 decorrelated legs vs a held-out external anchor. `occluded✓` = needs a consequence leg. Verdicts + fixes: `data/mechanism_catalog/verify_verdicts_wave1.json`.


## foundational  (2 cells)

### NODE-0-ENERGY  · **[WEAKENED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** WEAKENED — legs bound; ANCHOR QUESTION: '1st-law closure' is a consistency relation COMPUTED from the legs (heat+work=energy in) -> risks self-consistency not external held-out GT
- **claim:** Cell as free-energy transducer; 1st-law energy closure (heat + O2) is the measurable reduction anchor for the occluded regulatory/efficiency state.
- **hidden state:** cell energy balance / efficiency η(t)
- **legs:** respirometry (O2 flux) [ocr-stats-plos-s4-fibroblast, zdrazilova-2021-mitofit-fibroblast-zenodo] · calorimetry (heat) [weee-indirect-calorimetry-zenodo] · engineering-twin battery/Semenov (coupling, no bio data) [stub]
- **anchor (held-out):** 1st-law energy closure (measured)  ·  **mode:** SYNC
- **couples into CS:** P16 battery vertical / thermal (methodological twin)
- **regime:** single-cell/tissue bioenergetics; 1st-law closure via calorimetry+respirometry.

### KNEE-CELL  · **[WEAKENED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** WEAKENED (n_eff=2) — CORE certifiable claim = OpenSim DOF+ROM+moment chain predicts the GC instrumented-implant MEASURED contact force PER-SUBJECT (legit blind-prediction, non-circular, Fregly Grand Challenge). n_eff=2 NOT 4: leg2≡leg3 = one OpenSim chain; leg1(3dbodytex)+leg4(oai) are DISJOINT cohorts from the GC subjects = population priors, add 0 bits to the per-subject over-determination (leg4 also joint-state mismatch: OAI native cartilage vs GC TKA implant). Coordinator over-credited at n_eff~3 = rescue-reflex, caught by independent re-verify (which coordinator had gated on).
- **claim:** Is a knee pose/motion real and dynamically realizable — 4 decorrelated legs vs the Grand-Challenge in-vivo contact-force anchor.
- **hidden state:** knee pose realizable + real joint load
- **legs:** surface-geometry (optical body scan — INDEP substrate) [3dbodytex] · anatomical DOF+ROM (OpenSim knee model, bundled in GC) [simtk-knee-grand-challenge] · moment-realizability (OpenSim MTU+inverse-dyn — COMMON-MODE w/ leg2) [simtk-knee-grand-challenge, addbiomechanics] · cartilage/contact imaging (knee MRI — INDEP substrate) [oai]
- **anchor (held-out):** simtk-knee-grand-challenge (instrumented-implant in-vivo MEASURED contact force — independent telemetry, held out from the OpenSim predictions)  ·  **mode:** complementary-AND
- **couples into CS:** HUM-HAND / H (SMPL-X per-joint trust)
- **regime:** load/motion regime of the Grand-Challenge gait/stair tasks.


## musculoskeletal  (6 cells)

### MSK-FIBERTYPE  · **[FENCED]** · EMPIRICAL · risk=HIGH · occluded✓

- **verdict:** UN-FENCED (data found) — gse178734 snRNA + gse178733 snATAC decorrelated pair — [cartographer] human+rat muscle fiber-type snRNA (I/II signatures) + matched snATAC = 2 decorrelated modalities. Was FENCED (no data); now executable.
- **residual:** 57+95GB — mouth/ingest-daemon fetch, not a tick-bound subagent
- **claim:** Per-muscle fiber-type (I/IIa/IIx) + satellite reserve via transcript/proteome/EM legs vs histological imaging anchor; satellite sub-axis HYPOTHESIS-only (no functional anchor).
- **hidden state:** fiber-type fractions + satellite reserve
- **legs:** transcriptomic MyHC [gse130977-human-fibertype] · proteomic MyHC [pxd006182-human-fiber-proteome] · EM metabolic-phenotype [empiar-10070-muscle-mito, emdb-sarcomere-cryoet] · satellite state + load-intervention consequence [gse143704-human-muscle-atlas, gse252357-resistance-24h, gse114763-epigenetic-memory] · ncl-sm-myofiber (DEMOTED to segmentation-quality leg — NOT fiber-type GT) [ncl-sm-myofiber]
- **anchor (held-out):** NONE (FENCED — no ATPase/MyHC-IHC fiber-type histology + no force-consequence dataset in archive)  ·  **mode:** complementary-AND
- **regime:** adult human skeletal; satellite/regeneration needs a load/disuse intervention, not a resting atlas.

### MSK-DOF-ROM  · **[un-audited]** · EMPIRICAL · risk=MED · direct

- **verdict:** MEASURED-B — model ROM vs AAOS goniometry executed (mean|Δ|=12.52%, 9/12<20%) — OpenSim .osim `range` verified a-priori solver clamps (radians, 60% diff across authors); AddBiomechanics anchor = marker-fit-quality not ROM. Rajagopal2016 vs AAOS 9.8%, Gait2392 144.6%. Coordinator reproduced ranges + confirmed repo-clean.
- **residual:** simtk files login-gated; certifies non-knee ROM only for the actively-deployed model, not the Gait2392 template
- **claim:** Per-joint DOF+ROM for non-knee joints via geometric + cross-species mocap legs vs AddBiomechanics 273-subject real gait kinematics; extends inherited JOINT-DOF-REFERENCE.
- **hidden state:** per-joint achievable DOF + ROM envelope
- **legs:** anatomical/sim ROM bound (OpenSim family — a-priori solver clamps, n_eff=1) [simtk-per-joint-models, simtk-msk-size-ladder] · cross-species DOF check (independent) [pferd-horse-mocap] · clinical goniometric ROM (non-OpenSim, decorrelated) [aaos-goniometry-rom]
- **anchor (held-out):** addbiomechanics  ·  **mode:** complementary-AND
- **couples into CS:** JOINT-DOF-REFERENCE, H, PER-DOF-OBS
- **regime:** normal joints only; pathology is MSK-ORTHOPEDIC.

### MSK-ORTHOPEDIC  · **[RE-SCOPED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** RE-SCOPED — OrthoLoad-hip instrumented-implant anchor replaces the GRF common-mode; scope shrunk to HIP contact-force magnitude (n_eff=1) — GRF->load common-mode CONFIRMED (the only bridge IS leg-a itself, London Lower Limb static-optimization = same artifact not just same family). OrthoLoad (hip+shoulder instrumented implant) found as the real decorrelated internal-load anchor (mirrors KNEE). dryad-meniscus-cadaver DROPPED (confirmed KNEE-specific, off-regime). All live-HTTP-verified.
- **residual:** OrthoLoad access individually-gated (BT-HOLD, like CARDIAC/IMMUNE); SHOULDER pilot BLOCKED until the bundled simtk-per-joint-models id is split (shoulder=dsem project unbound); hip claim is n_eff=1 magnitude-only (same thinness as MSK-DOF-ROM), no hip injury/ligament dataset in-archive
- **claim:** Extend knee-cell occluded-joint-load paradigm to hip/shoulder+injury via geometry + cadaveric injury-mechanics legs vs GaitRec ~2300 healthy/pathology GRF; anchor STRUCTURALLY WEAKER than knee's instrumented implant.
- **hidden state:** joint contact/ligament load under injury (hip/shoulder)
- **legs:** geometry/moment-arm (London Lower Limb, n_eff=1 thin) [simtk-per-joint-models]
- **anchor (held-out):** orthoload (in-vivo instrumented-implant HIP contact force — real, decorrelated from the OpenSim geometry leg; individually-gated)  ·  **mode:** complementary-AND
- **depends_on:** MSK-DOF-ROM, MSK-CONNECTIVE
- **couples into CS:** sibling of KNEE-CELL; HUM-CONTACT-SIM
- **regime:** pathological/injured joints only (not knee normal-load).

### MSK-XSPECIES  · **[un-audited]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** MEASURED-B — 2-leg decorrelation fold DONE (tautology explicitly avoided) — gait2392 human PCSA 522.85cm2 @ BM=75.2kg sits +0.38 SD from the LOO(n=35, Rajagopal-excluded) amniote allometry, within 95% PI [161,1160]; robust across specific-tension 35-137. Coordinator verified the 40/40-Rajagopal-provenance tautology-catch + reproduced the residual.
- **residual:** agreement is magnitude/intercept match vs a NEAR-ISOMETRIC law (p=0.117, no anatomy beyond dimensional scaling); n_human=1; surrogate PCSA (validated vs published slope)
- **claim:** Cross-species MSK scaling law (moment-arm/DOF/force vs body-mass, mouse→macaque+equine+dog-breed) tested against the PROVEN human-specific-3-DOF JOINT-DOF-REFERENCE finding held out.
- **hidden state:** MSK moment-arm/force-capacity scaling law
- **legs:** geometric scaling (OpenSim gait2392, Delp lineage — Rajagopal REJECTED as anchor-source tautology) [opensim-gait2392-thelen2003] · equine locomotor consequence [pferd-horse-mocap] · genotype→morphometry [dog10k, plassais-canine-morphology]
- **anchor (held-out):** bishop-amniote-msk-allometry (PCSA/muscle-mass-vs-body-mass across 60 amniotes, CC-BY — DECORRELATED from OpenSim)  ·  **mode:** complementary-AND
- **depends_on:** MSK-DOF-ROM
- **couples into CS:** JOINT-DOF-REFERENCE, H
- **regime:** scaling valid over the body-mass range spanned; dog-breed = within-species extreme.

### MSK-CONNECTIVE  · **[RE-SCOPED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** RE-SCOPED — FIX EXECUTED 2026-07-18. (1) species-regime flag added to fairdata-ligament-tendon-mech leg (bovine ex-vivo; registry already said 'bovine knee' in resolver_notes but the cert claim/leg text carried no caveat — now explicit, qualitative structure-function use only). (2) dryad-knee-fe-s193761 CONFIRMED FE-simulation output, not directly-measured cadaveric loading: fetched Dryad record (doi:10.5061/dryad.zcrjdfnpv, 'Subject-specific knee models, data, and results for specimen S193761') — directly-MEASURED elements are kinematics/laxity/imaging (biplanar fluoroscopy, motion capture, 6DOF robotic-tester kinematics, laxity apparatus, CT/MRI/surface scan); ligament forces are Abaqus-Explicit FE-model outputs, optimized (Simplex/Particle-Swarm) to match those measured kinematic targets — calibrated/predicted, not sensor-measured. Validation paper (PMC12391196 / arXiv:2412.19951) states 'differences in the predicted ligament loads and calibrated material properties emerged' — 'predicted' confirms model-output status. Anchor relabeled derived/constructed. (3) searched for a directly-measured substitute: simtk-knee-grand-challenge (instrumented TKA-implant telemetry, real sensor, human in-vivo) exists in-registry and IS a genuine directly-measured contact-force anchor, but it is already KNEE-CELL's anchor (reuse here breaks cross-cell decorrelation) and measures TKA-implant contact force in a REPLACED joint, not native ligament tension or native cartilage (same joint-state mismatch the graph already flags for KNEE-CELL's own oai leg). No open, directly-measured, native-tissue ligament-force or cartilage-load dataset was found in a general web search (classical robotic/UFS in-situ-force studies, e.g. Fujie 1996 PMID:8884481, measure force directly via a universal force-moment sensor but were not deposited as open datasets with DOIs).
- **residual:** the quantitative tendon-tension/load-bearing-capacity number remains uncertified against any independent directly-measured, native-tissue, decorrelated anchor — stays HYPOTHESIS-strength; the omics/imaging/degeneration legs' qualitative structure-function cross-check is the only part of the claim standing. Open task: locate or newly deposit a directly-measured (load-cell/strain-gauge/buckle-transducer) native human tendon-ligament-force dataset.
- **claim:** Tendon-tension/cartilage load-bearing (occluded 'tendon under tension vs slack') via mechanobiology/ECM-proteomics/tensile/ultrasound legs vs cadaveric knee-FE ligament-force model.
- **hidden state:** tendon/ligament tension + cartilage load capacity
- **legs:** mechanobiology under known load [gse150482-tendon-stretch] · ECM proteomics [matrisomedb] · tensile mechanics (consequence) [SPECIES-REGIME: bovine ex-vivo — qualitative structure-function only, NOT quantitative human modulus] [fairdata-ligament-tendon-mech] · in-vivo ultrasound deformation [zenodo-tlf-fascia-ultrasound] · cartilage degeneration [gse104782-oa-chondrocyte, oai]
- **anchor (held-out):** dryad-knee-fe-s193761 [DERIVED/CONSTRUCTED — FE-model output, NOT directly-measured cadaveric loading]  ·  **mode:** complementary-AND
- **couples into CS:** HUM-DEFORMABLE (blocked)
- **regime:** omics legs NOT same-sample-paired with mechanical/imaging legs (thin axis, arch note 6).

### MSK-NEUROMUSCULAR  · **[RE-SCOPED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** RE-SCOPED — grabmyo anchors cross-subject EMG-pattern generalization; 'motor-unit force output' NOT certified (only hyser, un-held-out) = HYPOTHESIS. sEMG common-mode cluster flagged; ninapro/IMU carry decorrelation.
- **claim:** Motor-unit force + NMJ transmission (occluded 'muscle positioned vs force-generating') via same-sample EMG-force + GRF + kinematic + molecular-NMJ legs vs held-out cross-subject GRABMyo; NMJ leg has no decorrelated consequence-anchor (compensated NMJ failure invisible).
- **hidden state:** motor-unit force output + NMJ integrity
- **legs:** same-sample EMG-force (hyser; force sub-claim = HYPOTHESIS, un-held-out single-study) [hyser] · whole-body GRF/IMU (decorrelation weight) [locomotion-emg-kinetic] · kinematic-marker (encoder, genuinely decorrelated) [ninapro-db9-kinematics, dwivedi-semg-finger-kin]
- **anchor (held-out):** grabmyo (anchors EMG-PATTERN generalization only, NOT force)  ·  **mode:** complementary-AND
- **couples into CS:** HUM-HAND, HUM-CONTACT-SIM (occluded-truth-via-consequence instance)
- **regime:** hyser/ninapro legs + grabmyo anchor are cross-cohort, population-level not individual.

**honest gaps (musculoskeletal):**
- no in-vivo instrumented contact-force anchor for hip/shoulder (only knee has Grand-Challenge)
- no same-SAMPLE multimodal pairing in this cluster — all cross-cohort
- MSK-FIBERTYPE satellite sub-axis HYPOTHESIS-only
- MSK-NEUROMUSCULAR NMJ leg has no NMJ-specific consequence anchor
- canine SOD1-DM data exists (wave11) but is a neuro/ALS leg not MSK — boundary call
- several licenses unconfirmed


## energy_metabolism  (5 cells)

### MITO-CELL  · **[REFUTED]** · EMPIRICAL · risk=HIGH · occluded✓

- **verdict:** REFUTED — the same-sample r≥0.7 claim FAILS on the cell's own raw data (agent reproduced from Zdrazilova zenodo 5518059) — raw Pearson r (no shared normalization, n=3 HDF-batches = the ENTIRE paired population): ET-capacity/max-OXPHOS (claim target) = -0.226; robustness median +0.127 / incl-outlier -0.208 / Spearman +0.5 — ALL far below 0.7. n=3 (df=1) statistically near-uninformative anyway. External mechanism: source paper says O2k used multiple FCCP titrations vs XF24's 2 -> XF24 'E' reads low for PROTOCOL reasons not biology (quantity not cross-platform-comparable). Anchor-family-match CONFIRMED (ocr-stats source paper literally 'OCR-Stats...using Seahorse XF Analyzer').
- **residual:** needs (a) a genuinely open NON-Seahorse anchor (SOMMA 31P-MRS candidate is real but dbGaP-gated) AND (b) a same-sample dataset with >3 paired units before the r≥0.7 claim is even testable
- **claim:** Same-human-fibroblast Oroboros(amperometry) vs Seahorse(fluorescence) agree on max OXPHOS O2-flux (r≥0.7), within the 203-line ocr-stats percentile band held out. TRUE same-sample pair.
- **hidden state:** single-cell max OXPHOS flux capacity
- **legs:** electrochemical O2 electrode [zdrazilova-2021-mitofit-fibroblast-zenodo] · optical O2-quenching (SAME sample) [zdrazilova-2021-mitofit-fibroblast-zenodo]
- **anchor (held-out):** ocr-stats-plos-s4-fibroblast  ·  **mode:** SYNC
- **depends_on:** NODE-0-ENERGY
- **couples into CS:** NODE-0-ENERGY; P16 (rate-capability-under-load analogy)
- **regime:** FUNCTION layer only; structure (cristae SA/V→flux) NOT closeable — EMPIAR mouse-muscle mismatched to human-fibroblast respirometry. Structure leg FENCED.

### MEMBRANE-NERNST-CELL  · **[RE-SCOPED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** REFUTED-still (narrow RE-SCOPED at most) — no orthogonal 3rd leg exists AND the archive has no disruption arm — agent downloaded the actual zenodo-bbb-chip archive (456MB) + enumerated: TEER+Papp are same-state SYNC (n_eff=1, confirmed physics+lit) BUT there is NO barrier-disruptor/perturbation arm (only baseline model-char + TMZ/DOX drug-crossing) -> the magnitude/dose-response regression is NOT executable. Orthogonal 3rd leg (ion-selective/patch-clamp paired w/ TEER): searched 3 ways, NONE open. Anchor gse116222 is scRNA-only (no permeability assay); dryad-fd4 is mouse-only. This CORRECTS the prior over-optimistic RE-SCOPED which assumed a disruption arm without checking the raw archive.
- **residual:** only a same-instrument cross-source-type consistency check (iPSC vs HBMEC monolayers) is honestly supportable — explicitly NOT a disruption/dose-response cert. Needs a wet-lab TJ-disruptor arm or a not-yet-existing open ion-selective barrier dataset.
- **claim:** Same hiPSC-BBB-chip TEER(electrical) vs Papp(mass-transport) move oppositely under disruption; directional sign checked vs held-out gut barrier (FD4 + IBD). Targets barrier-INTEGRITY, not literal per-ion Nernst. [RE-SCOPED: hidden_state -> aggregate paracellular barrier-INTEGRITY (not per-ion Nernst — no ion-selective/patch-clamp in archive); test via magnitude/dose-response vs TEER-Papp transfer curve, not directional pass/fail]
- **hidden state:** barrier ionic-transport integrity
- **legs:** TEER (electrical) [zenodo-bbb-chip] · Papp (mass-transport, SAME chip) [zenodo-bbb-chip]
- **anchor (held-out):** dryad-fd4-gut-permeability + gse116222-ibd-barrier-breakdown (direction check)  ·  **mode:** SYNC (TEER+Papp = two transductions of ONE pore-state; 3rd orthogonal leg NONE-FOUND in archive)
- **depends_on:** NODE-0-ENERGY
- **couples into CS:** P16 EIS impedance (same interface-impedance principle)
- **regime:** TEER is aggregate paracellular, NOT ion-species-resolved; true per-ion Nernst (patch-clamp) not in archive.

### ADIPOSE-BAT-THERMOGENESIS-CELL  · **[HONEST-NEG]** · EMPIRICAL · risk=LOW · occluded✓

- **claim:** REGIME-GATE CLOSED (forced honest-negative): no archived dataset pairs cold-challenge indirect-calorimetry with cold-challenge browning transcriptome on same organism/protocol; human BAT omics absent (wave26). Verified via grep, not lazy.
- **hidden state:** UCP1 uncoupled thermogenic capacity
- **legs:** molecular browning trajectory [gse133486-browning-timecourse]
- **anchor (held-out):** weee-indirect-calorimetry-zenodo (regime-mismatched → invalid, flagged)  ·  **mode:** complementary-AND (UNACHIEVED)
- **depends_on:** NODE-0-ENERGY, MITO-CELL
- **couples into CS:** MITO-CELL; P1 conjugate-heat-transfer (regulated inverse of P16 runaway)
- **regime:** 1 real leg only; no valid same-organism/protocol anchor. regime_gate.open=false.

### AGING-ETA-DRIFT-CELL  · **[REFUTED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** REFUTED — η-DNAm coupling NOT detected (measured-null, robust; UNDERPOWERED n=4) — within-donor passage-detrended ρ=-0.21 p=0.31 (Horvath2 pre-reg); robust to η-def/covariate/permutation. 7-clock sweep sign-flips = noise, NOT cherry-picked. FigShare-403 was stale (API works); N-discrepancy resolved (RNAseq n=4 vs DNAm n=6)
- **residual:** n=4 underpowered (hFB14 only 2 pts) — 'not detected' not 'absent'; powered test needs more paired-lineage donors = FENCED-acquire
- **claim:** Same fibroblast lines across replicative lifespan: bioenergetic η(t)=ATP-linked/total OCR drifts monotonically (|ρ|≥0.5) correlating with same-sample DNAm-age; age-tracking checked (existence only) vs Hannum blood. [RE-SCOPED: n=4 HEALTHY-only (drop 3 SURF1), passage-detrended PARTIAL correlation, N reported at donor level (not 345 samples); Hannum = face-validity only (blood, cross-tissue) NOT magnitude anchor]
- **hidden state:** bioenergetic coupling efficiency η(t) vs replicative age
- **legs:** Seahorse ATP-rate (longitudinal same-cell) [sun-picard-fibroblast-lifespan] · DNAm clock (SAME cells) [sun-picard-fibroblast-lifespan]
- **anchor (held-out):** gse40279-hannum  ·  **mode:** SYNC
- **depends_on:** NODE-0-ENERGY, MITO-CELL
- **couples into CS:** MITO-CELL; P16 (DNAm-clock ~ non-destructive SOH proxy for destructive direct measurement)
- **regime:** sun-picard N-discrepancy unresolved (verify GSE179848); DNAm-age magnitude does NOT transfer cross-tissue, only qualitative. Cycling(passage) vs calendar(chronological) analog (P16) — don't conflate.

### CANCER-ROGUE-GROWTH-CELL  · **[RE-SCOPED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** RE-SCOPED-breast design, but EXECUTION -> WEAKENED (anchor correspondence sign-unstable/non-significant; leg2 Warburg clean but insufficient alone) — DepMap∩MetMap Jaccard 0.325 (83.8% shared); Sanger-adversary 78.7% (structural, not source-fixable); flux MCF-7-only. 2 effective legs (fused CRISPR+organotropism) vs msk-met breast subset n=2609. Coordinator independently confirmed the structural overlap (99% all-catalogued).
- **residual:** agent's /mnt/data_root raw cache zero-filled (corruption); exact 409 not disk-reproducible but computed pre-persist + direction independently confirmed
- **claim:** Cell-line CRISPR-essentiality + 13C-Warburg-flux + in-vivo organotropic metastasis (3 decorrelated legs) vs real-patient organotropism msk-met (25775 pt), ρ≥0.3. Flux leg is cell-line-only (no human-tumor ex-vivo flux, unlike node-0's normal-liver anchor).
- **hidden state:** degree growth escaped energy-budget regulation
- **legs:** CRISPR essentiality [depmap-ccle] · 13C Warburg flux [mtbls182, mtbls241] · in-vivo organotropic spread [metmap]
- **anchor (held-out):** msk-met-2021  ·  **mode:** complementary-AND
- **depends_on:** NODE-0-ENERGY, MITO-CELL
- **couples into CS:** MITO-CELL, NODE-0-ENERGY; P16 thermal-runaway (uncontrolled positive-feedback energy capture, structural analogy)
- **regime:** normal-vs-cancer anchor asymmetry: normal tissue has human ex-vivo flux anchor, cancer does not.

**honest gaps (energy_metabolism):**
- MITO-CELL structure-function EM+respirometry pairing FENCED (organism mismatch)
- no per-ion Nernst dataset — MEMBRANE cell is barrier-integrity proxy
- ADIPOSE-BAT not buildable (verified zero calorimetry+browning pairing; human BAT blank)
- sun-picard N-discrepancy unresolved (FigShare 403)
- no human-tumor ex-vivo 13C-flux
- P1/P16 are CS lane labels NOT ANCHOR_GRAPH node-ids — analogies only, never depends_on
- registry double-registration (GSE133486, MSK-MET) — dedup before ingest


## neural_sensory  (9 cells)

### EYE-OPTICS-FORWARD-MODEL  · **[un-audited]** · ENGINEERING · risk=LOW · direct

- **verdict:** corneal-power leg MEASURED on open NHANES (B) — 12685 eyes, mean 43.57D (n=1.3375); my formula matches NHANES native diopter field <0.02D (computation validated). Index-convention finding: 1.3375 clinical vs 1.376 anterior-only = ~5D (posterior surface)
- **residual:** spherocylindrical not OPD-raytrace (need an eye-raytracer, inherited is projector-only); NHANES-native 'anchor' is a formula self-check not the independent OCT-biometry cross-check (wholeeye gated); fellow-eye non-independence
- **claim:** Cornea+lens refractive power via corneal-surface raytrace (OPD/wavefront metric, normal-subset) cross-checked vs an INDEPENDENT paraxial schematic-eye from OCT biometry (different instrument+cohort). [REAL: anchor=wholeeye-oct-biometry resolves the Sim-K same-source tautology; retired keratometric-index clause pruned; CornOrb restricted to 889 normal eyes so keratoconus severity can't masquerade as raytrace error]
- **hidden state:** eye refractive power/wavefront (what the retina receives)
- **legs:** corneal-surface raytrace (OPD enforced; NORMAL 889-eye subset only, exclude keratoconus) [cornorb-topography]
- **anchor (held-out):** wholeeye-oct-biometry (paraxial schematic-eye from OCT biometry — different instrument+cohort than cornorb topography, resolves the tautology; population-level not same-eye)  ·  **mode:** SYNC
- **couples into CS:** LENS-FIELD-RAYTRACE-SIM (real inherited node)
- **regime:** not yet same-eye (cornorb vs oct-biometry different cohorts → population-level SYNC). MUST use OPD/wavefront metric not spot-RMS (position-only-metric-blind-opd-catches); smallest-positive-t root.

### VISION-LOSS-STAGE-LOCALIZATION  · **[RE-SCOPED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** RE-SCOPED — per-subject UKB design (audit's 'sensor-vs-relay' fix was STALE/unconstructible) — ⚠ audit FIX text stale: no relay leg ever existed (HCP-7T=cortex=central, prior-session correction); GBD/IAPB has NEITHER a central/CVI NOR a relay/optic-atrophy category (WHO/IAPB primary-source verified) -> anchors sensor ONLY. Re-scoped to a SAME-SUBJECT UKB retina-OCT + brain-MRI design (N=6446, Sun 2024): retinal thickness vs thalamus(relay-proxy) vs occipital(central). GBD demoted to sensor prevalence sanity-check.
- **residual:** NOT executed: UK Biobank is DAA-gated (not open-HTTP). Coordinator should re-fetch PMID 38716763/40885450 independently before trusting N=6446/108. relay proxy is whole-thalamus not LGN-specific (UKB 1mm res); Chen-2025 open-status unconfirmed. Both relay+central FENCED vs GBD, un-fenceable vs UKB pending execution.
- **claim:** Localize vision-loss stage (sensor/RGC, relay/optic-nerve, central/cortex) by over-determining intact stages, ruling out optical cause via EYE-OPTICS first. [CORRECTED: this localizes sensor(RGC)-vs-CENTRAL(cortex); there is NO relay/optic-nerve source, so 'relay' cannot be adjudicated; GBD/IAPB anchor + disjoint-cohort legs => population-ecological, not per-patient]
- **hidden state:** which visual-pathway stage failed
- **legs:** sensor/RGC structural [idrid, kermany-oct2017, refuge-glaucoma-suite, gse137400-rgc-atlas] · central/cortical retinotopic map topology (HCP 7T = visual CORTEX, not relay/optic-nerve) [hcp-7t-retinotopy]
- **anchor (held-out):** ukb-retina-brain-mri (per-subject sensor→relay→central; GBD/IAPB demoted to sensor-only prevalence sanity-check)  ·  **mode:** complementary-AND
- **depends_on:** EYE-OPTICS-FORWARD-MODEL, BRAIN-CONNECTOME-CENTRAL-COMPUTE
- **couples into CS:** EYE-OPTICS-FORWARD-MODEL
- **regime:** GBD/IAPB anchor is population cause-of-blindness, not per-patient; RP/IRD has no open corpus (sensor leg can't cover it).

### HEARING-LOSS-STAGE-LOCALIZATION  · **[un-audited]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** MEASURED-B — ARHL genetics localize to STRIA VASCULARIS > sensory epithelium (resolves the field question) — GWAS-locus->cochlear-cell-type bridge established via 3 converging peer-reviewed ARHL-GWAS x cochlear-scRNA studies (Wells 2019 44 loci, Trpchevska 2022 48 loci [machine-verified 48==48 vs GWAS-Catalog GCST90105058], Shi 2026 140 loci); field disagreement resolved by Eshel/Elkon 2024 (sensory epithelium). Auditory-nerve stage = Ahmed 2026 (spiral ganglion neurons).
- **residual:** significance only at the 500kb window (direction robust 5/5, magnitude window-dependent); mouse->human naive symbol uppercasing (no ortholog table); cross-platform marker confound (shared ~equally by both compartments); central stage still fenced (no open corpus)
- **claim:** Localize deafness stage via hair-cell structural+genetic vs ABR-latency-pattern + human audiometry; CENTRAL stage UNANCHORED (no corpus). [FIX: narrowed to sensor-vs-relay; 'central' DROPPED (unanchored); GWAS-locus->cochlear/spiral-ganglion cell-type bridge required before ukb-arhl-gwas is a stage anchor]
- **hidden state:** which auditory-pathway stage failed
- **legs:** hair-cell structural+genetic [gse135913-cochlea, dvd-deafness-variation] · ABR + audiometry (cross-species) [zenodo-abr-mouse, nhanes-audiometry]
- **anchor (held-out):** ukb-arhl-gwas  ·  **mode:** complementary-AND
- **regime:** sensor-vs-relay machine-checkable via ABR peak-I; central (cortical deafness) 2/3 anchored only.

### PARALYSIS-STROKE-STAGE-LOCALIZATION  · **[REAL]** · EMPIRICAL · risk=LOW · occluded✓

- **verdict:** REAL design; consequence leg MEASURED (B) — BCI-decode consequence leg strong (z=26.8, central signal survives relay failure = occluded-truth-via-consequence working); soop NIHSS anchor confirmed n=1106
- **residual:** atlas-r2-stroke per-subject BLOCKED (Google-form gate, brief misclassified as open); structural-lesion leg deferred
- **claim:** Localize paralysis relay(SCI/nerve) vs central(stroke). OCCLUDED-TRUTH worked example: DANDI:000147 intracortical BCI in tetraplegic proves central signal survives despite relay failure.
- **hidden state:** central signal survival vs relay-distribution deficit
- **legs:** central signal via movement-decode consequence [dandi-000147-tetraplegic-bci] · lesion structural mapping (stroke MRI+mask) [atlas-r2-stroke]
- **anchor (held-out):** soop-openneuro + nscisc-sci-registry (clinical NIHSS/mRS + ASIA severity, exam-derived — decorrelated from the raw MRI/mask leg)  ·  **mode:** complementary-AND
- **depends_on:** BRAIN-CONNECTOME-CENTRAL-COMPUTE
- **couples into CS:** PER-DOF-OBS
- **regime:** PRO-ACT ALS = negative control (LMN/effector failure, both central+relay intact) — must not misclassify.

### PERIPHERAL-ENTERIC-AUTONOMIC-EDGE-COMPUTE  · **[REFUTED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** REFUTED — anchor & legs are DISJOINT animal cohorts (per-animal LOO structurally impossible) — execution revealed a deeper defect than the anchor∈legs tautology I fixed: anchor sparc-vagus #287 = UCL pigs (Pig15/22/24/25, domestic) vs leg gse154411/gse149212 = Jefferson Yucatan MINIpigs (1534/1643/1705/1729) vs gse130710 = MOUSE — NO subject-ID/breed/institution overlap on 4 axes. A per-animal correspondence needs the SAME physical animal on both sides; it doesn't exist. Forced index-pairing diagnostic: rho=-1.0 but exact-perm p=0.333 (n=3, zero power). Data fully parseable (real per-pig EIT + SAN-fraction extracted) — the blocker is design-level disjointness, not parsing.
- **residual:** needs a SINGLE study with both intrinsic-cardiac-NS cell-type architecture AND VNS-cardiac functional response in the SAME animals (neither #287 nor GSE154411 provides both jointly) -> cartographer hole
- **claim:** Enteric NS + intrinsic cardiac NS show local sensory-inter-motor complement for edge-compute in 2 organ-distinct systems vs DRG relay-only negative control; SPARC vagus VNS→cardiac-effect as functional anchor.
- **hidden state:** local reflex-processing capacity independent of CNS
- **legs:** enteric architecture vs DRG null [enteric-ns-atlas, gse168243-human-drg] · intrinsic cardiac + vagal consequence (sparc-vagus REMOVED — was anchor∈legs) [gse154411-ragp-intrinsic-cardiac-ns, gse130710-san-pacemaker]
- **anchor (held-out):** sparc-vagus  ·  **mode:** complementary-AND
- **regime:** certifies STRUCTURAL precondition (cell-type complement), NOT measured autonomous firing (no isolated-ganglion ephys).

### SOMATOSENSORY-TOUCH-PROPRIOCEPTION-HAND  · **[un-audited]** · EMPIRICAL · risk=HIGH · occluded✓

- **verdict:** MEASURED-B — afferent-type separation measured (186 units) — FA-1/SA-1/SA-2 separable on adaptation-index F=34.5 (z=33 vs perm-null), canonical direction FA>SA-1>SA-2 matches textbook; subject-demean adversary holds (z=18). md5+186-units+real-labels QC-verified.
- **residual:** COORDINATOR caveat: FA-vs-SA axis partly DEFINITIONAL (types named for adaptation). Non-tautological = SA-1/SA-2 gradation + quantitative textbook match. DRG-transcriptome/hand-pressure cross-link not yet folded
- **claim:** Hand-surface pressure (21 ADL) vs DRG mechanoreceptor-type composition predicting frequency-band sensitivity. TOUCH ONLY — proprioception ZERO coverage; NO anchor (no human microneurography) → stays HYPOTHESIS.
- **hidden state:** mechanoreceptor afferent activity ↔ skin pressure
- **legs:** hand pressure (ADL) [zenodo-hand-tactile-adl] · DRG afferent substrate [scp2370-hdrg, gse168243-human-drg]
- **anchor (held-out):** saal-microneurography (single-unit human tactile-afferent recordings — the missing anchor, FOUND)  ·  **mode:** complementary-AND
- **regime:** no human single-unit microneurography anywhere in archive — 2 legs but NO anchor.

### BRAIN-CONNECTOME-CENTRAL-COMPUTE  · **[WEAKENED]** · DATA-ANCHOR · risk=MED · direct

- **verdict:** WEAKENED-still — species/modality conflation ALREADY fixed (h01-only human-EM leg); cross-scale no-shared-subject defect remains — human-brain-cell-atlas drop was real; microns NOW structurally removed from the leg (was prose-only 'not counted' = invisible to the flattener). Residual (pre-fix, unresolved): cross-scale has NO shared subject (macro MRI cohorts / 1 EM specimen / BigBrain histology) -> population/regional only
- **residual:** cross-scale per-subject over-determination is impossible (EM destructive); keep as regional-consistency claim only
- **claim:** Cross-scale connectome (macro MRI-tractography+fMRI vs micro/meso EM-synaptic) vs BigBrain cytoarchitecture — the CENTRAL-COMPUTE baseline consumed by vision/paralysis localization cells.
- **hidden state:** cross-scale structural connectome consistency
- **legs:** macro connectome [hcp, uk-biobank-abcd-adni] · micro/meso EM connectome (human, counted) [h01-human-cortex-em]
- **anchor (held-out):** bigbrain  ·  **mode:** SYNC
- **regime:** no subject has both macro-MRI + EM (EM destructive) → cross-scale check is population/regional, caps per-patient precision downstream.

### FACIAL-MUSCLE-EMOTION-CONSEQUENCE  · **[REFUTED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** REFUTED-still — disjoint-subject test RESOLVED the residual, AGAINST the claim — agent executed the exact open residual (disjoint-subject re-measure of osf-dux3w, md5-verified). In-sample R²=0.34 (NOT the claimed 0.78 — inflated). Raw-coord LOSO catastrophic (-49) = subject-geometry confound; after per-subject centering (fair adversary), expression-SHAPE generalizes (LOSO 0.43, 6/6 folds). DECISIVE: the sEMG-activation-SPECIFIC increment over an expression-identity-only baseline is null-to-NEGATIVE (LOSO -0.097 2/6; L2O -0.119, 3/15 folds positive, t≈-3.6) -> the cell's mechanism (motor-nerve drives INDIVIDUAL muscles, generalizably) FAILS a fair disjoint-subject test; only 'which expression' generalizes (weaker claim).
- **residual:** leg-B rescue blocked: biovid-partb-facial-emg is the modality-correct independent anchor but DUA-gated (Univ Ulm signed agreement, not open); animalfacs-primate = taxonomy manual (0 data); rat-grimace = pain images (no EMG). REFUTED-still: needs an OPEN same-sample sEMG+facial-shape dataset with a real disjoint-subject anchor to even re-attempt.
- **claim:** Facial CN-VII drive (occluded) via sEMG-driven ArtiSynth face-model consequence (R²=0.78) vs independent facial-EMG generalized cross-species (primate/rat FACS). Couples to CS H (per-joint→per-muscle trust). [FACS promoted to CO-anchor (video-AU, EMG-independent) to break the surface-EMG common-mode between leg A and the biovid anchor]
- **hidden state:** facial motor-nerve drive to individual muscles
- **legs:** sEMG->ArtiSynth mechanics consequence [osf-dux3w-semg-artisynth, artisynth-face-models]
- **anchor (held-out):** biovid-partb-facial-emg (human EMG) + animalfacs-primate + rat-grimace-dataverse (video-AU CO-anchor — EMG-independent, breaks the surface-EMG volume-conduction common-mode)  ·  **mode:** complementary-AND
- **couples into CS:** H
- **regime:** R²=0.78 held-out-vs-training UNVERIFIED — confirm before promoting. Both legs surface-EMG → EMG cross-talk is a shared failure mode (common-mode, flagged).

### CHEMOSENSORY-OLFACTION-GUSTATION  · **[HONEST-NEG]** · EMPIRICAL · risk=HIGH · occluded✓

- **verdict:** HONEST-NEGATIVE -> taste anchor FOUND (partial unblock) — nhanes-taste-csx gives a 3708-subject taste-psychophysics anchor (was zero taste data). Still needs a taste-receptor molecular leg + an olfaction 2nd leg
- **residual:** acquire taste-bud/receptor molecular data + odor-perception; then taste half is a 2-sided cert
- **claim:** FORCED HONEST-NEGATIVE: olfaction has 1 uncrossed leg, no 2nd leg, no anchor; taste has ZERO data. Not composable.
- **hidden state:** odor-receptor response ↔ perceived odor (taste: none)
- **legs:** olfactory epithelium (uncrossed) [gse139522-olfactory]
- **anchor (held-out):** NONE  ·  **mode:** n/a
- **regime:** single uncrossed leg; taste is the cluster's clearest zero-coverage gap.

**honest gaps (neural_sensory):**
- taste ZERO coverage
- proprioception ZERO coverage
- central-auditory-loss unanchored (2/3 taxonomy)
- RP/IRD no open corpus
- no same-eye cornea+biometry pairing (EYE SYNC is population-level)
- no isolated-ganglion ephys (edge-compute = structural precondition only)
- no human microneurography (SOMATOSENSORY has no anchor)
- ALS needs own cell
- Bell's palsy no data
- CornOrb Sim-K schema unconfirmed
- ArtiSynth R² held-out status unverified


## organ_systemic  (14 cells)

### ORG-WAVEFORM-ANCHOR  · **[REAL]** · DATA-ANCHOR · risk=LOW · direct

- **verdict:** REAL — MEASURED/PROVEN — live probe: 6/6 sources HTTP 200; VitalDB carries SpO2(6386)/EtCO2(12581)/ECG(9745) case-tracks (machine-verified from api.vitaldb.net/trks, 486k tracks). Confirms the EtCO2 dependency ORG-LUNG needs.
- **residual:** downstream consuming cells (LUNG/HRV) still need their own execution; this is the anchor-availability proof only
- **claim:** Whole-body physiological waveform corpus EXISTS and is largely OPEN (VitalDB CC-BY intraoperative multi-channel; MIMIC-IV Waveform DB ODbL ICU ECG/ABP/PPG/resp; PTB-XL CC-BY 21799 12-lead ECG; PhysioNet HRV/autonomic-aging/Fantasia ODC-BY) -- the shared CONSEQUENCE-anchor infrastructure every occluded-organ cell below draws its external anchor from (operator physiology-waveform framing).
- **hidden state:** whole-body multi-organ operating state (the shared substrate this node provides, not a single organ)
- **legs:** 
- **anchor (held-out):** vitaldb, mimic-iv-waveform, ptbxl-ecg, physionet-hrv-suite, physionet-autonomic-aging, physionet-fantasia  ·  **mode:** SYNC
- **regime:** Prefer the fully-open no-credential subset (VitalDB, PTB-XL, PhysioNet HRV suites) for first cert passes; MIMIC-IV waveform signal itself is open, only its clinical-annotation companion is credentialed. This is aggregate whole-body state, not organ-specific -- usable as anchor ONLY where the organ's function has a distinct waveform channel (SpO2/EtCO2 for lung, RR-interval for cardiac autonomics); does not anchor kidney/liver/endocrine/immune cells.

### ORG-CARDIAC-PLAQUE  · **[un-audited]** · EMPIRICAL · risk=LOW · occluded✓

- **verdict:** BLOCKED at cert (outcome anchor + calcium leg both DUA/RUA-gated, verified live) — plaque-composition LEG measured — inter-patient plaque-macrophage-composition heterogeneity is REAL (Bashore within-study Cochran chi2=2471 df=17 p~0, batch/sorting/tissue held constant, 20k-binomial-null), replicated across 3 vascular beds (fails only in n=3,4 cohorts). Marker-grounded, sorting-confound-audited. Coordinator reproduced Bashore chi2 exactly.
- **residual:** Athero-Express MACE (dataverse DOI:10.34894/4IKE3T, DSA-gated) + COCA calcium (Stanford RUA/Redivis 401, gated) both verified live -> no individual-outcome fold; open EHJ-2026 macrophage-burden->MACE = population-level anchor only
- **claim:** Certify atherosclerotic plaque state (cellular composition + calcification burden) as a predictor of real clinical outcome (rupture/MI/stroke), not just cross-sectional atlas description.
- **hidden state:** plaque composition (foam cell/fibrous cap/macrophage subset ratio) and vulnerability to rupture
- **legs:** plaque_molecular_atlas (hca-heart-atlas dropped: healthy donor hearts, not plaque) [traeuble-plaque-atlas] · calcification_imaging [coca-coronary-calcium]
- **anchor (held-out):** athero-express  ·  **mode:** complementary-AND
- **regime:** Athero-Express outcome anchor is DSA-controlled (measured, not yet fetchable) -- design is sound, acquisition pending. Leg (a) is internally COMMON-MODE (both scRNA/snRNA); leg (b) CT is the real decorrelating force.

### ORG-HRV-AUTONOMIC  · **[RE-SCOPED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** RE-SCOPED — design internally sound (anchor independent of its OWN legs); fix applied = shared-anchor-pool tag — the WEAKENED flag was a GRAPH-LEVEL independence-accounting issue (its physionet anchor pool is verbatim-shared with ORG-WAVEFORM-ANCHOR), not a within-cell defect. Fix = tag both cells with shared_anchor_pool so an 'N independently-anchored cells' count does not over-credit them. Cell design unchanged, still executable.
- **residual:** cross-cell independence must discount the shared physionet-cardiac pool (1 effective anchor across HRV + WAVEFORM, not 2)
- **claim:** Certify cardiac autonomic control (SAN pacemaker firing + vagal/sympathetic modulation) via its measurable population-level consequence: heart-rate variability.
- **hidden state:** SAN pacemaker cell firing rate + vagal/sympathetic tone on intrinsic cardiac ganglia
- **legs:** pacemaker_substrate_atlas [gse130710-san-pacemaker, gse154411-ragp-intrinsic-cardiac-ns] · vagal_stimulation_response [sparc-vagus]
- **anchor (held-out):** physionet-autonomic-aging, physionet-fantasia, physionet-hrv-suite  ·  **mode:** complementary-AND
- **depends_on:** ORG-WAVEFORM-ANCHOR
- **regime:** Fully open end-to-end (no controlled-access dataset in this cell). Leg (a) pools two snRNA sources (pacemaker + ganglion) -- COMMON-MODE on assay; leg (b) SPARC vagus stimulation-response is the true orthogonal force.

### ORG-IMMUNE-TOLERANCE  · **[un-audited]** · EMPIRICAL · risk=LOW · direct

- **verdict:** BLOCKED at cert (anchor individual-GT DUA-gated, verified live) — both open LEGS measured + adversary-forced — central-tolerance stringency (mTEC-II vs mTEC-I depth-deconfounded, z=9.03) + islet beta-depletion (donor-level, Control 40.6%->T1D 17.1%, matches nPOD literature). Exemplary honest KILL of the depth-artifact single-gene AIRE proxy. Coordinator reproduced the islet numbers exactly.
- **residual:** TEDDY(dbGaP phs001442)/DAISY(NIDDK) individual GT DUA-gated -> no individual held-out cert; only a population fold vs krischer-2015 aggregate incidence (open). gse243061 leg-component source-corrupted (verified)
- **claim:** Certify the self/non-self discriminator's calibration point -- central tolerance (thymic AIRE-mediated deletion) vs its failure (autoreactive escape) -- against real longitudinal autoantibody seroconversion preceding clinical autoimmune onset.
- **hidden state:** thymic negative-selection stringency + peripheral autoreactive-clone escape probability
- **legs:** central_tolerance_substrate [e-mtab-8581-thymus-dev, gse243061-aps1-aire-treg] · islet_multimodal_readout [hca-t1d-islet]
- **anchor (held-out):** teddy, daisy  ·  **mode:** complementary-AND
- **regime:** TEDDY/DAISY anchor is dbGaP/registration-controlled (measured, real longitudinal cohort). hca-t1d-islet's internal CyTOF+IMC arms are the BEST same-sample decorrelation found in the whole cluster -- prioritize that pairing over cross-cohort legs elsewhere.

### ORG-TRANSPLANT-GRADIENT  · **[RE-SCOPED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** RE-SCOPED — platform-confound fix specified (rank-based within-platform AUC + pseudobulk); execution BLOCKED (xeno anchor DAC-gated) — fix = rank-based single-sample enrichment AUC on the common-coverage gene-set (2005 2-color array = floor) killing the platform-era confound; +pseudobulk-aggregate any single-cell/spatial before rank-scoring (removes a resolution confound the audit didn't name). Monotonicity self<allo<xeno UNTESTED.
- **residual:** xeno anchor EGAS50000000244 is DAC-gated (verified: EGA page 'DATA ACCESS COMMITTEE') + phs001667 dbGaP controlled -> cannot execute without out-of-band grant. 3 GEO series + IPD-IMGT-HLA are fetchable but can't close the xeno term. Public substitute xeno atlas may exist (gds Count=65) — follow-up.
- **claim:** Certify that the self/non-self rejection signature scales monotonically across the allo->xeno genetic-distance gradient (kidney/heart/gut allograft -> pig->human kidney xenograft), the direct readout of the body's self-cert threshold (architecture note 3).
- **hidden state:** graft-infiltrating immune rejection intensity as a function of donor-recipient genetic/species distance
- **legs:** allograft_rejection_signature [gse145927-kidney-chimerism, gse2596-heart-rejection, gse134662-gvhd-gut] · hla_genetic_distance [ipd-imgt-hla, phs001667-kidney-genomics]
- **anchor (held-out):** egas50000000244-pig-human-xeno  ·  **mode:** complementary-AND
- **depends_on:** ORG-IMMUNE-TOLERANCE
- **regime:** PRE-REGISTERED falsifier: rejection-signature effect size/margin must be LARGER in xeno than allo (monotonic in genetic distance) -- a non-monotonic or reversed result falsifies the 'one discriminator, one gradient' claim, not just this cell. Xeno anchor is EGA-controlled (measured, real). HLA leg genuinely orthogonal (genotype vs expression).

### ORG-INFECTION-REPERTOIRE  · **[WEAKENED]** · EMPIRICAL · risk=MED · direct

- **verdict:** WEAKENED (executed) — infection-response REAL but severity-gradient NULL; OAS↔COMBAT severity-join structurally blocked — OAS Disease=SARS-COV-2: 104.8M seqs / 1977 units / 10 published studies, Vaccine=None (natural infection), CC-BY (machine-verified: sum of per-row unique-seqs == page total). CoV-AbDab (12,918 rows, reconciles banner) confirmed to have NO patient-ID/clonal-freq/repertoire-size field -> relabeled a CDR3-keyed antigen-SPECIFICITY side-input (known-binder lookup only), excluded from repertoire-magnitude. COMBAT anchor: DERIVED per-patient BCR feature tables (78 samples x 48 metrics, WHO-severity) openly on Zenodo 6120249 (no-auth GET verified); raw FASTQ EGA-gated.
- **residual:** OAS raw-count is concentration-SKEWED (Kim_2020 = 61.7% of total) -> must use per-study/per-individual-NORMALIZED metrics, NOT pooled raw count, as the repertoire-magnitude statistic; COMBAT raw VDJ FASTQ remains EGA-gated (only derived features open)
- **claim:** Certify adaptive-immune repertoire response (BCR+TCR clonal expansion + systems-vaccinology transcriptome) against real measured infection/vaccination clinical severity. [FIX: repertoire leg restricted to OAS repertoire-seq; re-verify COMBAT anchor access/license before firm]
- **hidden state:** antigen-specific clonal expansion + systems-level transcriptional response magnitude
- **legs:** antibody_repertoire (OAS repertoire-seq subset only; CoV-AbDab = antigen-specificity side-input, NOT counted) [oas-covabdab] · tcr_and_systems_vaccinology [vdjdb-ireceptor, hipc-immunespace]
- **anchor (held-out):** combat-covid  ·  **mode:** complementary-AND
- **regime:** COMBAT anchor license unconfirmed + portal ECONNREFUSED this pass (flagged in registry) -- treat as PROVISIONAL anchor until re-verified live.

### ORG-KIDNEY-FILTRATION  · **[HONEST-NEG]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** WEAKENED (executed honest-negative) — leg1(transporter expr)↔leg2(eGFR GWAS) correspondence NULL, robust — ckdgen actionability RESOLVED (EBI mirror GCST90103633 md5-verified; portal genuinely NXDOMAIN). KPMP real subsegmental RNAseq range-fetched (4MB/2.9GB, CRC32-verified). leg1↔leg2 Spearman rho=0.098 p=0.766 (200k perm), robust across GTEx cross-check (0.224) + 3 windows (all p>0.48). Void-floor: anchor NOT dead (11/18 metabolites 6/6-consistent p=7e-13) -> null is meaningful.
- **residual:** structural null-space: 10/12 candidate genes (Na-transporters+aquaporins) have NO substrate in ST002785's HILIC-neg 18-metabolite panel (Na/H2O not LC-MS analytes) -> only SGLT/glucose anchor-testable (also null). Anchor is perfusate-only (tissue-metabolism/net-exchange, not pure filtration). Honest-negative: transporter EXPRESSION does not track eGFR GENETICS.
- **claim:** Certify per-nephron-segment filtration/reabsorption capacity (structural atlas + donor-recipient genomic risk) against a directly measured ex vivo perfusion metabolic functional readout.
- **hidden state:** per-nephron-segment transport/filtration capacity
- **legs:** nephron_segment_atlas (aggregated to whole-kidney) [kpmp] · baseline filtration-function genomics (decorrelated) [ckdgen-egfr-gwas-stanzick2021]
- **anchor (held-out):** mw-st002785-kidney-perfusion  ·  **mode:** complementary-AND
- **regime:** Anchor (mw-st002785) is the rare TRUE functional measurement (ex vivo normothermic perfusion of an intact kidney) in this cluster -- CC BY 4.0, open. No in-vivo GFR/clearance dataset located (honest gap).

### ORG-LIVER-ZONATION  · **[WEAKENED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** WEAKENED — re-exec RESOLVED the QC contest: 9/10 robust (NOT the claimed 10/10, NOT my 6/10 collapse), OAT donor-heterogeneous + SLC1A2 genuine null — independent 3rd re-execution (9 donors parsed, n=3040 exact) adjudicated both prior claims. Periportal pole 5/5 robust (never disputed). CYP1A2/RHBG/CYP3A4 correct-sign WITHOUT donor-batch (my earlier pipeline undercounted them). Only OAT strictly needs donor-batch to flip — and it's DONOR-HETEROGENEOUS: P301+P304 (63% cells) drive the positive, P310 SIGNIFICANTLY CONTRADICTS (rho=-0.215 p=2e-4) even post-correction; within-donor sign split 5-vs-3 (LOO stays positive). Permutation-null confirms the flip needs the TRUE donor grouping (not over-fit). SLC1A2 = genuine null (rho=+0.0018 p=0.92) that the prior agent mislabeled a 'win'.
- **residual:** honest count 9/10 robust + 1/10 null; the pericentral pole is REAL but weaker/donor-heterogeneous than a clean 10/10 -> WEAKENED not MEASURED-B. Species mismatch resolved (human GSE124395); bulk-flux consistency partial (2/3).
- **claim:** Certify hepatic lobule zonation (periportal-oxidative vs pericentral-detox metabolic division of labor) against directly measured 13C metabolic flux in intact human liver ex vivo.
- **hidden state:** zone-specific hepatocyte metabolic flux partition
- **legs:** spatial_zonation_map (HUMAN) [gse124395-human-liver-zonation] · mechanical_stiffness_imaging (whole-organ corroboration only) [scidb-mre-brain-liver]
- **anchor (held-out):** mtbls10481-human-liver-13c  ·  **mode:** complementary-AND
- **regime:** Cleanest cert in the cluster: all three sources open/CC-BY, legs genuinely orthogonal assay technology (RNA-seq vs mechanical imaging), anchor is a real flux MEASUREMENT not an inference. Supplementary leg geo-gse93382-gh-pulse (hormone-driven liver transcriptome dynamics) available but shares leg (a)'s RNA-seq modality -- not counted as a 3rd independent leg.

### ORG-LUNG-GASEXCHANGE  · **[REAL]** · EMPIRICAL · risk=LOW · occluded✓

- **verdict:** REAL design; disease-leg MEASURED (B) — IPF alveolar AT1/AT2 loss robust (subject-level p=0.00014, Bonferroni-survives, ext-anchor Aberrant_Basaloid); agent caught pseudoreplication trap (cell-level chi-sq false positives)
- **residual:** disjoint cohorts (mode 10): disease-leg vs VitalDB gas-exchange anchor = population-level fold only, not individual
- **claim:** Certify alveolar-capillary gas-exchange efficiency (cellular atlas + CT structure) against the measured whole-body consequence: arterial SpO2/EtCO2/respiratory waveforms.
- **hidden state:** alveolar-capillary gas exchange efficiency
- **legs:** airway_alveolar_atlas [hlca, gse136831-ipf-copd, gse102580-ionocyte] · structural_ct_imaging [lidc-idri]
- **anchor (held-out):** vitaldb, mimic-iv-waveform  ·  **mode:** complementary-AND
- **depends_on:** ORG-WAVEFORM-ANCHOR
- **regime:** No pulmonary-function-test (spirometry/FEV1) dataset located in archive -- anchor substitutes gas-exchange OUTCOME (SpO2/EtCO2) for airflow-MECHANICS per the occluded-truth rule; related but not identical hidden state, flagged as a partial proxy (honest gap).

### ORG-ENDOCRINE-AXES  · **[REAL]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** REAL (scoped) — analytic-leakage defect RESOLVED: SYNC non-circular by construction — GSE122541=Resuehr 2019 (n=6 nurses). Transcriptomic phase = JTK_CYCLE/MetaCycle (hormone-agnostic, generic algos; input = expression matrix + real sample clock-time only); cortisol/melatonin analyzed separately (cosinor/Rayleigh). NO DLMO/hormone-peak re-referencing of the transcriptome time axis (the leakage pattern hunted + confirmed absent, 3 independent primary-source fetches). SYNC is a genuine decorrelated-instrument confirmation.
- **residual:** scope = n=6 adrenal cortisol/melatonin axis ONLY (small n); structural_crossgland_atlas leg (thyroid/pituitary/parathyroid) has NO functional hormone time-series = honest gap unchanged. Design fixed; n=6 SYNC measurement not yet executed.
- **claim:** Certify the HPA/circadian hormone axis (pilot: adrenal cortisol rhythm) using SAME-SUBJECT paired clock-gene transcriptome vs measured serum cortisol/melatonin -- the cluster's one genuine same-sample calibration case (wave-21 gold standard).
- **hidden state:** circadian phase / disruption state of the HPA axis
- **legs:** peripheral_clock_transcriptome [gse122541-shiftwork] · structural_crossgland_atlas [hormone-cell-atlas, adrenal-cortex-dev, gse142653-human-fetal-pituitary, gse233962-parathyroid]
- **anchor (held-out):** gse122541-shiftwork  ·  **mode:** SYNC
- **regime:** Only the adrenal/cortisol axis closes this way; thyroid/pituitary/parathyroid remain structural-atlas-only with no functional hormone-level time series located (honest gap). SYNC mode chosen deliberately: both sides estimate the SAME target (circadian phase/disruption), redundantly, across genuinely different instruments (RNA-seq vs immunoassay) -- not a tautology since the anchor is a different measurement instrument on the same subjects, not the same number.

### ORG-DIGESTIVE-ENZYME  · **[REAL]** · EMPIRICAL · risk=LOW · occluded✓

- **verdict:** REAL design; MEASURED components (B-grade) — real cell-level detail measured (acinar 0.112, enterocyte 0.39, per-gene expr); anchor cgmacros externally validated (iAUC vs A1c ρ=0.71, BMI-null). Coordinator verified data+counts; ρ agent-reported not yet reproduced
- **residual:** DISJOINT COHORTS (mode 10, like knee): legs & anchor share no subjects -> population-level fold only, not individual over-determination. Reproduce ρ + fetch stomach/gutcellatlas legs
- **claim:** Certify digestive enzymatic processing (chief-cell pepsinogen, exocrine pancreas, brush-border absorption) via the measured whole-body postprandial glycemic consequence.
- **hidden state:** luminal digestion + nutrient absorption rate
- **legs:** digestive_cell_atlas [gutcellatlas, hca-stomach, gse84133-pancreas, gse125970-nutrient-absorption] · enzyme_kinetics [brenda-sabio-digestive]
- **anchor (held-out):** cgmacros  ·  **mode:** complementary-AND
- **regime:** BRENDA/SABIO kinetics are CONSTRUCTED/aggregated in-vitro purified-enzyme measurements, not raw single-experiment data -- used only as supporting mechanism evidence, NEVER as the anchor (measured-vs-constructed discipline). cgmacros is genuinely measured, open, CC-BY-NC-SA.

### ORG-BARRIER-BBBGUT  · **[RE-SCOPED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** RE-SCOPED — cross-species generalization; human functional gut-permeability = FENCED (no lactulose-mannitol/Ussing in archive)
- **claim:** Certify barrier integrity (BBB tight-junction + gut epithelium) via functional permeability assays, then test generalization to a HELD-OUT barrier-organ+disease context (gut IBD breakdown) -- mirrors the reduced-rep generalization criterion in the verified core. [RE-SCOPED: explicit cross-species test (human BBB-chip -> held-out MOUSE gut permeability); same-species human functional validation is a FENCED gap]
- **hidden state:** tight-junction/epithelial barrier permeability
- **legs:** barrier_structural_atlas (human scRNA, common-mode disclosed) [gse163577-vineseq-bbb, gse116222-ibd-barrier-breakdown] · functional_permeability_assay (human hiPSC) [zenodo-bbb-chip]
- **anchor (held-out):** dryad-fd4-gut-permeability (MOUSE — explicit CROSS-SPECIES generalization test of the human BBB-chip mechanism)  ·  **mode:** complementary-AND
- **regime:** This is a cross-ORGAN generalization test (BBB-trained legs -> gut-disease anchor), not a same-organ cert. Falsifier: the barrier mechanism fails to transfer (no signal) when applied to the held-out gut context -- would indicate 'barrier integrity' is organ-specific, not a general discriminator.

### ORG-BLOOD-COAG  · **[HONEST-NEG]** · EMPIRICAL · risk=LOW · direct

- **claim:** Certify hematopoietic lineage integrity (incl. megakaryocyte/platelet output feeding primary hemostasis) against plasma coagulation factor levels -- flagged as the WEAKEST cell in the cluster.
- **hidden state:** hematopoietic differentiation fidelity + coagulation factor output
- **legs:** hematopoietic_atlas [hca-hematopoietic-atlas] · malignant_disruption [van-galen-aml]
- **anchor (held-out):** charge-coagulation-gwas  ·  **mode:** complementary-AND
- **regime:** HIGH-risk/incomplete cert: both legs share scRNA assay technology (no true orthogonal 2nd leg found in archive this pass) -- NOT a genuine complementary-AND yet, effectively single-mechanism-class. No functional clotting-time (PT/aPTT/thromboelastography) dataset located. DEFER promotion until a functional coagulation-assay dataset is acquired; do not treat as PASS on current evidence.

### ORG-SKIN-BARRIER  · **[RE-SCOPED]** · EMPIRICAL · risk=MED · direct

- **verdict:** RE-SCOPED — split into 2 coherent sub-cells (category-error fixed) — (a) STRUCTURAL-AGING: hca-skin-aging (HES1/KLF6 senescence decline, 9 donors) vs gse241132 wound-roadmap — reframed from 'aging predicts healing speed' to a MOLECULAR-PROGRAM correspondence (anchor has NO age covariate, 3 donors — GEO-verified). ⚠ leg+anchor both transcriptomic = single-modality common-mode. (b) LESION-MALIGNANCY: ham10000 dermoscopy vs PAD-UFES-20 (biopsy-proven BCC/SCC/MEL, live-verified) — fixes the neoplasia-vs-wound category error.
- **residual:** (a) single-modality common-mode + no age covariate + body-site match unconfirmed (eyelid leg vs anchor site); (b) biopsy-proof only for malignant classes (restrict scoring to them) + dermoscope-vs-smartphone domain-shift confound must be controlled + ISIC must exclude HAM10000-source (tautology risk, like EYE-OPTICS/FACIAL). Physical 2-cell split deferred; documented here.
- **claim:** Certify skin structural aging + lesion phenotype against a real measured wound-healing process trajectory. Skin is the LEAST occluded organ in this cluster (directly imageable/biopsiable) -- used as an intentional low-occlusion contrast case against the solid-organ cells above.
- **hidden state:** epidermal/dermal structural integrity + lesion malignancy risk
- **legs:** skin_aging_atlas [hca-skin-aging] · dermoscopic_imaging [ham10000]
- **anchor (held-out):** gse241132-wound-roadmap  ·  **mode:** complementary-AND
- **regime:** No transdermal-permeability functional assay (TEER-analog) located for skin specifically (unlike BBB's zenodo-bbb-chip) -- anchor is a PROCESS outcome (healing trajectory), not a barrier-integrity instrument (honest gap). Wound-roadmap access is mixed (GEO open / EGA companion controlled).

**honest gaps (organ_systemic):**
- Same-sample multi-modal pairing (the wave-21 gold-standard cert-calibration axis) is achieved in only 2/14 cells here: ORG-IMMUNE-TOLERANCE (hca-t1d-islet's internal CyTOF+IMC+RNA) and ORG-ENDOCRINE-AXES (gse122541-shiftwork's transcriptome+cortisol+melatonin same subjects). Every other cell's two legs are cross-cohort/cross-study decorrelation -- a weaker cert form than a true same-sample pairing.
- Cross-cutting COMMON-MODE risk: a majority of 'structural/molecular' legs in this cluster (heart, plaque, SAN, ganglion, thymus, hematopoietic, AML, skin, BBB, gut/stomach/pancreas, liver, lung atlases) are droplet-based scRNA/snRNA -- one shared assay-technology mechanism-class (ambient-RNA, dissociation-stress genes, batch/platform effects). Every cell's 2nd leg was deliberately chosen non-transcriptomic (imaging, kinetic assay, electrophysiology, genotype, functional permeability) to force this; within a merged leg-family it remains ONE piece of evidence, not N -- flagged per-cell.
- ORG-BLOOD-COAG has NO genuine 2nd orthogonal leg (both legs are scRNA) and no functional clotting-time (PT/aPTT/TEG) dataset was located in the archive -- this cell should be treated as DEFERRED, not PASS-eligible, until a functional coagulation assay is acquired.
- No pulmonary-function-test (spirometry/FEV1) dataset located -- ORG-LUNG-GASEXCHANGE substitutes the whole-body SpO2/EtCO2 gas-exchange OUTCOME for airflow-MECHANICS; related but not identical hidden state.
- No thyroid- or pituitary-specific FUNCTIONAL hormone-level time series located (only cross-sectional structural atlases) -- ORG-ENDOCRINE-AXES only closes the adrenal/cortisol axis; thyroid/pituitary/parathyroid remain FILLED-structural-only.
- Several designed anchors are MEASURED but controlled-access, not yet fetchable: athero-express (DSA), egas50000000244-pig-human-xeno (EGA DAA), teddy/daisy (dbGaP/registration), phs001667-kidney-genomics (dbGaP), charge-coagulation-gwas individual-level (dbGaP; summary open). Design is sound; acquisition is a separate, tracked step.
- BRENDA/SABIO digestive-enzyme kinetics are CONSTRUCTED/aggregated in-vitro measurements (curated database), not a single raw experiment -- used only as supporting leg evidence in ORG-DIGESTIVE-ENZYME, never as an anchor, per the measured-vs-constructed discipline.
- No transdermal-permeability functional instrument (TEER-analog) exists for skin specifically in the archive (unlike BBB's zenodo-bbb-chip) -- ORG-SKIN-BARRIER's anchor is a process outcome, not a barrier-quantification instrument.
- combat-covid (ORG-INFECTION-REPERTOIRE's anchor) had an unconfirmed license and portal ECONNREFUSED on the most recent registry pass (wave 24) -- treat as provisional until re-verified live.


## molecular_program  (12 cells)

### MOL-GENOME-VARIANT-GROUNDED-ANCHOR  · **[REFUTED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** PROVEN — circularity-refuted AND GENERALIZES across genes (BRCA1/TP53/PTEN/MSH2, median AUROC 0.975) — functional-DMS-vs-ClinVar AUROC 0.94-0.99 across 4 genes, blind-snapshot circularity-refuted on 3/4 (BRCA1/TP53 solid, MSH2 refuted w/ measurable PS3-inflation, PTEN blind inconclusive). Decorrelated (small rho) from the frequency leg everywhere. Coordinator reproduced all per-gene AUROC exactly.
- **residual:** per-gene groundable via MaveDB (~hundreds of DMS score-sets, not all genes); PTEN sparse-benign (n=7) + inconclusive blind; MSH2 PS3-inflation over time; clinical-label anchor tier; ClinVar ancestry-enriched
- **claim:** Variant-effect predictor (FILLED decoder) scored ONLY vs POST-training-cutoff ClinVar (held-out, CAMEO-style); comparator bounded by gnomAD-frequency + functional-assay legs. Frequency-agreement alone ≠ validation (common-mode: frequency is a training feature). [RE-SCOPED: SINGLE-leg (allele-freq vs post-cutoff ClinVar, PS3-subset); ChEMBL/BindingDB DROPPED (mis-slotted, not variant-DMS); 2-leg claim WITHDRAWN pending a real variant-DMS corpus (ACQUIRE-gap)]
- **hidden state:** true variant pathogenicity
- **legs:** population allele-frequency (epidemiological) [gnomad-v4] · functional DMS effect (biochemical, frequency-independent) [mavedb, findlay-brca1-sge]
- **anchor (held-out):** clinvar-variant-summary (held-out P/LP-vs-B/LB; cross-checked 97% vs MaveDB's own ClinVar controls)  ·  **mode:** complementary-AND
- **regime:** restrict ClinVar anchor to post-cutoff submissions or it's a tautology gate; gnomAD license contested.

### MOL-PATHWAY-NETWORK-FILLED-DECODER  · **[un-audited]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** MEASURED — passes essentiality-matched adversary (B+); mode-11 resolved — Reactome co-membership -> DepMap co-essentiality survives an ESSENTIALITY-QUANTIL-MATCHED null 7x (not the pan-essential strawman floor), cross-checked vs BioGRID-raw physical (matched rate-ratio 9x, p~0). Anchor now genuinely external (BioGRID-raw fetched directly, mode-11 bundle avoided)
- **residual:** restrict to more pathways beyond glycolysis/TCA/spliceosome; the confound is real (naive 58x->9x) — always report the matched-null number
- **claim:** Reactome-pathway co-members show DepMap CRISPR co-essentiality above random floor, cross-checked vs BioGRID low-throughput physical interactions. Common-mode caveat: restrict to classical pre-2010 core pathways.
- **hidden state:** true gene-gene functional relationship
- **legs:** curated pathway co-membership (Reactome ONLY) [reactome] · CRISPR co-essentiality [depmap-ccle]
- **anchor (held-out):** string-biogrid (BioGRID raw low-throughput physical MITAB file ONLY, excl STRING combined/coexpr/textmining) — now external to leg A  ·  **mode:** complementary-AND
- **depends_on:** MOL-GENOME-VARIANT-GROUNDED-ANCHOR
- **regime:** Reactome/STRING/BioGRID/DepMap all confirmed-open; KEGG/MSigDB license traps excluded.

### MOL-PROTEIN-STRUCT-DYNAMICS-AF-HONESTY  · **[un-audited]** · EMPIRICAL · risk=LOW · occluded✓

- **verdict:** MEASURED — calibration + apply BOTH pass (A-, generalizes) — CALIBRATION pLDDT-vs-MD-RMSF ρ=0.65 (60 proteins) + APPLY to 8 novel post-cutoff structures ρ=0.608 (matches calibration, beats void-floor 19σ, homology-filtered). The decoupled design works: AF confidence predicts ACTUAL crystallographic flexibility on structures released AFTER training = the honesty layer functioning
- **residual:** n=8 novel small (attrition from SIFTS-lag/duplicate-UniProt, not signal); pLDDT partly self-referential but paired with independent crystal B-factor
- **claim:** AlphaFold static prediction within conformational envelope of MD-ensemble(ATLAS)+NMR/IDP(PED), scored vs CAMEO3D held-out target released AFTER training cutoff — NEVER vs its own training PDB (the tautology gate this cell forbids). = the improved-AlphaFold honesty layer. [FIX: DECOUPLED — calibrate AF-error-vs-flexibility on ATLAS/PED/mdcath, then APPLY to CAMEO via orthogonal pLDDT/B-factor (no ATLAS/PED membership required); machine-count intersection first; add seq-identity filter + IDP void floor]
- **hidden state:** true dynamic conformational ensemble
- **legs:** experimental static structure [rcsb-pdb] · MD ensemble (CALIBRATION stage) [atlas-md, mdcath] · NMR/IDP ensemble (CALIBRATION stage) [ped-idp-ensembles] · flexibility indicators (APPLICATION: pLDDT/B-factor on CAMEO) [alphafold-db, rcsb-pdb]
- **anchor (held-out):** cameo3d-heldout (post-cutoff)  ·  **mode:** complementary-AND
- **regime:** CAMEO never stale by construction; a model NEVER certifies its own FILL; ATLAS CC-BY-NC.

### MOL-SAMESAMPLE-PATCHSEQ-CALIBRATION  · **[un-audited]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** MEASURED — passes null (B+), morphology-free fix applied — ARI 0.177 vs perm-null z=16.6 p=0.002, robust over 16 configs + 5 seeds; external anchor (subclass->morphology textbook) corroborates; confounds fell. FIRST promotion-calibrator to pass execution
- **residual:** HONEST: joint T×E does NOT beat T-type-alone (0.177 vs 0.227) — ephys no incremental signal here; core cert passes but joint-beats-transcriptome NOT shown. n=135 modest; repo-processed not raw DANDI
- **claim:** FLAGSHIP calibration cell: patch-seq ephys+transcriptome joint cell-type call predicts held-out classical (Petilla) morphological type >2x floor. The mechanism that PROMOTES transcriptome-only atlas labels from BT-HOLD → GROUNDED.
- **hidden state:** single-neuron true cell-type identity
- **legs:** electrophysiology [dandi-000020-mouse-patchseq, dandi-000023-human-patchseq] · transcriptome (SAME cell) [dandi-000020-mouse-patchseq, allen-celltypes-patchseq-portal]
- **anchor (held-out):** morphology (Petilla, held-out)  ·  **mode:** SYNC
- **depends_on:** MOL-GENOME-VARIANT-GROUNDED-ANCHOR
- **regime:** DANDI open; striatal patch-seq missing → generality untested outside cortex.

### MOL-SAMESAMPLE-MULTIOMIC-CALIBRATION  · **[REAL]** · EMPIRICAL · risk=LOW · occluded✓

- **verdict:** REAL — MEASURED (B-grade, coordinator-verified ARI=0.358≫null z=266) — ARI≈0.376 (non-ambiguous 6493 cells) ≫ perm-null; TRIPLE-verified (agent + independent replication bit-exact + robust across seed/res z=215-273). My earlier 0.358 was a filter bug (included 'Unassigned' cells) — caught by the replication.
- **residual:** run the CROSS-DATASET design (10x legs vs gse158013 anchor) + independent replication before PROVEN
- **claim:** Same-cell RNA+ATAC cluster predicts held-out CITE-seq surface-protein identity vs canonical CD-marker criteria above floor.
- **hidden state:** single-cell regulatory state (ATAC→RNA→protein chain)
- **legs:** chromatin accessibility [10x-pbmc-multiome] · transcriptome [10x-pbmc-multiome]
- **anchor (held-out):** gse158013-teaseq-trimodal (ADT/protein channel via MANUAL CD-threshold gate, NOT Azimuth-WNN — disjoint from the 10x-multiome legs, resolves train-on-test)  ·  **mode:** SYNC (RNA+ATAC = one coupled 10x-multiome co-assay; decorrelation burden is at the disjoint anchor, met)
- **depends_on:** MOL-GENOME-VARIANT-GROUNDED-ANCHOR
- **regime:** no fully-independent flow-only cohort → CD-marker CRITERION substitutes; true clinical-flow is ACQUIRE gap.

### MOL-SAMESAMPLE-RADIOGENOMICS-CALIBRATION  · **[un-audited]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** MEASURED-B — radiogenomic correspondence executed on n=130 (prolif-sig vs %GG rho=-0.36) — 3-way join (TCIA imaging ∩ GSE103584 RNA-seq ∩ TCIA clinical): 61.6% for full n=211 BUT 100% unique+exact for the n=130 that actually has all 3 channels (0 dup/orphan/whitespace defects). The 81 shortfall is COHORT DESIGN (AMC arm=imaging+clinical only, no RNA-seq), NOT ID-ambiguity -> the 'join-noise fabricates agreement' defect is REFUTED for the used population. CT-precedes-surgery confirmed (all 130 interval>0, median 39d) = radiology independent of tissue channels. ★ Cross-validated vs Bakr 2018 Sci Data descriptor: 211/162/49/130 EXACT match (external over-determination).
- **residual:** cert admissible ONLY for the n=130 subcohort (not the full 211 the original cert implied); histo/RNA-seq are deliberately same-sample (that IS the calibration target — independence not claimed for those, only radiology-vs-tissue)
- **claim:** REGIME-GATE CLOSED: the NSCLC patient-ID imaging↔omics JOIN (falsifier owed from wave21, claim-on-page not executed) must be machine-verified (≥90% unique join) BEFORE any cert. Conditional: radiomic subtype vs RNA-seq subtype vs histopathology.
- **hidden state:** tumor molecular state ↔ imaging phenotype
- **legs:** imaging radiomics [nsclc-radiogenomics] · molecular RNA-seq [nsclc-radiogenomics]
- **anchor (held-out):** nsclc-radiogenomics (histopathology)  ·  **mode:** SYNC
- **regime:** gate closed until join executed — name-collision risk; do NOT train imaging↔omics pair first.

### MOL-DNA-EPIGEN-PROGRAM  · **[un-audited]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** MEASURED-B — joint decorrelated fold DONE: two structurally-disjoint molecular clocks converge on chronological age — blood-trained DNAm clock RANKS fibroblast donor-age cross-tissue: Spearman 0.83 exact-perm p=0.004, jackknife-robust; LOO offset-calibrated MAE beats null 0.25x. Raw-MAE test fails ONLY on a known ~20-27y tissue offset (calibration artifact, offset-immune checks confirm real signal). Non-circular (published clocks applied not re-fit). Distinct from AGING (η-DNAm null)
- **residual:** population-level (cross-cohort, no shared subjects — not an individual same-sample fold); minor inherited QC flag (2 HGPS 'Nyr' regex-miss samples, doesn't affect R2)
- **claim:** Blood-trained DNAm-clock generalizes cross-tissue to predict fibroblast donor-age (MAE<0.5x null), same lineages' Seahorse-OCR decline independently correlates same age ordering. Both legs must clear floors.
- **hidden state:** true biological-age/regulatory-drift state
- **legs:** DNAm clock [gse40279-hannum, methaging-clocks] · independent transcriptomic-age (GSE113957, non-Sun-Picard) [gse113957-fibroblast-age]
- **anchor (held-out):** sun-picard replicative-passage/donor-age  ·  **mode:** complementary-AND
- **depends_on:** MOL-GENOME-VARIANT-GROUNDED-ANCHOR
- **couples into CS:** AGING-ETA-DRIFT-CELL
- **regime:** genuine cross-tissue domain-shift test; Sun-Picard rare open exception (load-bearing).

### MOL-GROWTH-BODYSIZE-NUTRITION-ALLOMETRY  · **[un-audited]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** MEASURED-B — allometry direction executed (exp 0.752, R2=0.92) + sign agrees WASH-B; LOW discriminating power (honest) — causal magnitude anchor = WASH-B RCT (LAZ +0.25 SD CI[0.15,0.36], pre-reg PRIMARY, void-floor passes). Cross-species 519-mammal allometry kept as DIRECTION-ONLY sign cross-check, explicitly NOT pooled into a magnitude with the causal leg (different hidden state: cross-species structural vs within-species developmental — analogy not shared-state). WHO demoted to directional-FLOOR. INCAP dropped. Illusory-openness adversary forced (fetched PubMed, LAZ confirmed pre-reg primary not mined secondary).
- **residual:** WASH-B anchor is published AGGREGATE not raw-microdata-reproduced (Dataverse/GitHub have only code); MAL-ED observational corroboration not re-verified (budget). Cannot certify a single unified causal-magnitude spanning allometry+nutrition legs.
- **claim:** RANDOMIZED INCAP nutrient-dose→growth slope directionally consistent with observational MAL-ED + cross-species allometric exponent (519 mammals) vs WHO growth standards. Sign-disagreement falsifies causal claim.
- **hidden state:** true causal nutrient→growth effect
- **legs:** randomized intervention [incap-guatemala] · cross-species allometry [dryad-mammal-intestinal-allometry, anage-longevity]
- **anchor (held-out):** wash-benefits-bangladesh-rct (open causal RCT, LAZ +0.25 SD — replaces normative-WHO + inaccessible-INCAP)  ·  **mode:** complementary-AND
- **regime:** INCAP is contact-only access — HIGH risk stalls on acquisition.

### MOL-LIFESTAGE-PUBERTY-MENOPAUSE  · **[WEAKENED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** WEAKENED — puberty anchor doesn't contain the GT — rsna-bone-age has ONLY bone-age (Greulich-Pyle) + sex, NO Tanner staging (categorical mismatch); menopause half firmer (SWAN FMP is genuine) but possible ReproGen/SWAN sample-overlap unchecked
- **residual:** rescope puberty anchor to 'bone-age skeletal maturation (radiographic, NOT Tanner)' or source a Tanner-staged+genotyped cohort; confirm SWAN not in ReproGen's 21-study discovery meta (leakage check)
- **claim:** ReproGen genetic timing score (spans menarche+menopause) predicts held-out bone-age pubertal timing AND SWAN menopausal timing vs clinical event (Tanner/FMP). Both life-stages must clear floors independently.
- **hidden state:** HPG-axis timing program (2 life-stages)
- **legs:** physiological marker (bone-age/hormone) [rsna-bone-age, swan-icpsr-dbgap] · reproductive genetic architecture [reprogen-menarche-menopause]
- **anchor (held-out):** clinical timing event (Tanner/FMP)  ·  **mode:** complementary-AND
- **regime:** puberty & menopause are DIFFERENT cohorts — genetics legitimately spans both, not same-sample.

### MOL-REPRODUCTION-DEVELOPMENT-EMBRYOLOGY  · **[REFUTED]** · EMPIRICAL · risk=MED · occluded✓

- **verdict:** REFUTED — cross-atlas claim NOT supported (accuracy vs decorrelated ARI contradict) — accuracy 0.576 z=2.49 beats null BUT decorrelated ARI (z=0.11) does NOT -> symmetric-QC rejects the positive. Root cause: CS7 is a singleton confounded with assay technology (Tyser SMART-seq2 vs Xu 10x) -> apparent signal is within-Xu, not cross-atlas. Label withheld from fitting (verified); mild HVG transductive exposure flagged. Confirms audit's FILLED-not-leg
- **residual:** run the interrupted Xu-only 7-specimen decomposition (removes confounded CS7) to isolate any genuine same-technology signal; a fair cross-atlas test needs multiple specimens/stage/technology
- **claim:** scRNA developmental-stage classifier recovers classical Carnegie stage (morphology-only, co-registered same specimens) above floor on held-out.
- **hidden state:** molecular program driving staged morphology
- **legs:** single-cell transcriptome [tyser-gastrulation, xu-early-organogenesis] · classical morphological staging [tyser-gastrulation, xu-early-organogenesis]
- **anchor (held-out):** Carnegie-stage (morphology, held-out)  ·  **mode:** SYNC
- **regime:** Carnegie predates scRNA → decorrelated though same-specimen; no blind-unstaged validation exists.

### MOL-MICROBIOME-METAGENOMICS-HOST-COUPLING  · **[un-audited]** · EMPIRICAL · risk=HIGH · occluded✓

- **verdict:** MEASURED-B — same-sample LEG strong (z=17.4); host-outcome ANCHOR definitively MARGINAL (3 falsifiers converge, not rescued) — exclusion of iHMP_IBDMDB_2019 is 100% enforceable (per-study directory-skip); its DOI = Lloyd-Price 2019 = the exact anchor study, and it's the ONLY match -> dropping it removes the whole overlap. Anchor now disjoint (curators flag no iHMP2 overlap = audited). Same-sample pairing verified on RAW data (FRANZOSA 220/220 IDs identical)
- **residual:** 3 anchor falsifiers RAN + converged z~1.4-1.6 (real-but-weak/underpowered, none clear 1.645); metabolite-COUPLING step adds NO host-signal beyond raw taxonomy (delta AUROC +0.002) -> same-sample prediction is the solid result, metabolite->host coupling only marginally supported. Investigation CLOSED as honest-marginal.
- **claim:** Same-sample metagenome→fecal-metabolite prediction corroborated by external host outcome (IBD-activity/CGM) on held-out. Same-sample pairing = the microbial-metabolite→host cert-calibration axis.
- **hidden state:** true microbe→host phenotype effect
- **legs:** metagenomic profile (microbiome-metabolome RESTRICTED: EXCLUDE Lloyd-Price/iHMP_IBDMDB_2019=HMP2; retain Yachida/Franzosa-PRISM/Sinha/He/Jacobs/Poyet/Erawijantari/Kim/Mars/Kang/Kostic/Wandro/Wang) [microbiome-metabolome, hmp1] · fecal metabolome SAME sample (same restricted collection, HMP2 excluded) [microbiome-metabolome]
- **anchor (held-out):** ibdmdb-hmp2 (host outcome) — genuinely disjoint once the HMP2/Lloyd-Price-2019 study is excluded from the paired metabolome legs at acquisition (per-study directory-skip, audited; the anchor study is NOT among the retained leg cohorts)  ·  **mode:** SYNC
- **regime:** GutMIND flagship excluded (no verifiable URL), not substituted.

### MOL-PHARMACOLOGY-DRUGRESPONSE-TOX  · **[un-audited]** · EMPIRICAL · risk=LOW · occluded✓

- **verdict:** MEASURED — full MoA cert passes (B+), mega-cluster adversary forced — 629 compounds->138 MoA clusters (L1000-978 + PRISM-viability, CCLE-covariate-regressed) share ToxCast target-activation kappa=0.115 z=4.86 vs perm-null; EXCLUDING the 3 mega-clusters (80% of pairs) z ROSE to 18.85 = signal not an artifact. GCTX-transpose bug caught by own assert; ToxCast read in-memory (infra rule); self-test + from-scratch bit-exact reproduction
- **residual:** clustering silhouette weak (0.032, chaining) — the signal is real but the MoA partition isn't crisp; a cleaner clustering (or supervised MoA labels) would sharpen it
- **claim:** LINCS L1000 transcriptomic-perturbation + CTRP/PRISM viability agree on MoA, cross-checked vs independent ToxCast/Tox21 target-activation for held-out compounds above floor. [FIX: L1000 restricted to Landmark-978; report {LINCS+CTRP}∩{ToxCast} compound overlap; regress out shared CCLE covariates before scoring MoA]
- **hidden state:** true compound MoA/perturbation effect
- **legs:** transcriptomic perturbation (L1000 Landmark-978 MEASURED only, NOT imputed transcriptome) [lincs-l1000] · viability dose-response [gdsc-ctrp-prism]
- **anchor (held-out):** toxcast-tox21  ·  **mode:** complementary-AND
- **regime:** clue.io retiring (use GEO GSE92742 anchor not portal); ToxCast cleanest license (US public domain).

**honest gaps (molecular_program):**
- radiogenomics patient-ID join NOT executed (falsifier owed, wave21)
- no HDX-MS structure-indexed DB (no 3rd dynamics leg for many proteins)
- striatal patch-seq missing
- gnomAD license contested
- KEGG/MSigDB license traps excluded
- human-aging energetics deposit-thin (Sun-Picard load-bearing)
- IGF-1/GH-axis no verified dataset
- INCAP contact-access only (highest single-dependency risk)
- no human pubertal-window hypothalamus scRNA
- no independent flow-cytometry cohort
- GutMIND unverifiable — excluded
- clue.io/GDSC/DrugBank live-access unstable — cite durable accessions
