# DENTAL-CHAIN-PLANNER-20260928 — GAPS_RANKED

Read with `RESULTS.md` (argument and method) and `COVERAGE.json` (machine-readable).
All numbers were read from local artifacts at build time; see `inherited_execution_record` in `COVERAGE.json`. Nothing here is a clinical claim.

## What the measurement says before any judgement

- 1516 jobs, 1020 reports, 1414 producer jobs, 102 planner jobs.
- All 1414 producer jobs target **6 target nodes** drawn from **8 source directories** (8/8 of the catalog, and 0 outside it).
- Draft mechanism nodes: 386 OPEN, 105 ASSUMED, 13 PROVEN, 40 REFUTED.
- Category split of producer work: {'adversarial': 387, 'mechanism': 414, 'calibration': 84, 'integration': 126, 'audit': 371, 'exploratory': 32}. Audit+adversarial is 758/1414 = 53% against a planner instruction of 20 %.
- **Intake gate open: False** — 496 jobs without a report vs ceiling 480; runnable 365; reserve target 320. Unadmitted proposal files: 42.

The confinement is measurable, not an impression: the loop is not exploring "the existing eight source categories" broadly, it is running 1414 jobs against six nodes, and its refill path has been closed.

## The ten highest-value chain gaps

### 1. Occlusal load case: no patient-linked bite force exists  `G01` — CRITICAL, interface_break, stage L04

Every stress, fatigue, wear and design decision in the system consumes a bite force that no public source can supply for a specific patient. 16_occlusion_jawmotion.md answers the load question with "no" and enumerates six adversarial counter-arguments, all of which hold. The consequence is that DENT-LOAD-MAGNITUDE-SPECTRUM carries an 80x epistemic band (250-20000 cycles/yr above 500 N) and 100+ downstream nodes inherit it silently.

- **Smallest decisive next step.** Publish the single acquisition specification (measured bite transform on Teeth3DS+ arch pairs plus a contact-area x force-per-site table) and PAUSE all derived load-spectrum derivations. One specification, then move on.
- **Falsifier.** Bite registration residual exceeds the 0.3 mm contact band, so a patient-linked table cannot be built at design resolution and the system must instead carry an explicitly population-level load with a stated loss of design resolution.

### 2. Process-property calibration has zero measured pairs; AM-Bench is the route  `G02` — CRITICAL, acquisition, stage L07/L08/L10

The material and process layer is the most quantitatively developed part of the graph and the only part with NO measurement behind it. 73 fatigue jobs and 141 grade jobs re-derive the same inherited power law and eigenstrain model. Meanwhile 17_fabricate_test_pairs.md lists public-domain, DOI-tagged, geometry+measurement pairs (NIST AM-Bench AMB2018-01/03/04, AMB2022-04, AMB2025-02/03) that the loop has never used.

- **Smallest decisive next step.** Acquire one AM-Bench pair end to end and fit the process-property relation the model already parameterises. This converts an entire branch from ASSUMED to measurable.
- **Falsifier.** The fitted property is uncorrelated with the declared process variable, so AM-Bench process-property structure does not transfer and the dental ceramic route must be measured directly.

### 3. Biological stage has no admissible mechanism after three refutations  `G03` — CRITICAL, broken_branch, stage L13

Absolute Frost bands: REFUTED. Tooth-referenced mechanostat: REFUTED. Ortho PDL-pressure driver: REFUTED. 55 DENT-PAT nodes remain OPEN. The biological stage cannot change a design decision, so the lifecycle terminates at the tissue boundary.

- **Smallest decisive next step.** Test a stimulus-SPACE form of the mechanostat against a public load+remodeling dataset with a held-out anchor, or declare an empirical risk term. Do NOT produce a fourth Frost variant.
- **Falsifier.** Even the space form fails to recover the sign, proving the mechanostat family is wrong for this geometry.

### 4. Intake gate is closed: 42 written proposals silently unadmitted  `G04` — CRITICAL, operational, stage L15

The controller admits followups only while (jobs without RESULTS.md) < intake_ceiling. That count is 499 against a ceiling of 480, so the gate is false every tick. 63 proposal files exist; 21 were ingested; 42 sit unadmitted. The self-refilling loop has been stalled at the intake boundary, not the queue boundary, and STATUS.json does not show it.

- **Smallest decisive next step.** Gate intake on the runnable backlog (368), consistent with how MASSIVE_QUEUE_POLICY states targets are "valid research backlog, not concurrency", rather than on the raw unfinished count. Then re-run and confirm the 42 files are admitted.
- **Falsifier.** Admitting them pushes the queue past what the hosts can run, in which case the gate should instead be raised explicitly rather than redefined.

### 5. Prosthesis classes: removable is entirely absent, implant is 1 decision in 1414  `G05` — HIGH, scope, stage L05/L13

All 550 draft nodes and all 1414 producer jobs address a single case (implant-supported). 1 of 1414 decisions mentions implant/thread. 0 mention removable. Support physics for a tissue-borne prosthesis (mucosal compression, retention, wear of the retentive element) has no node at all.

- **Smallest decisive next step.** Open the class as a first-class dimension in the coverage matrix and register one removable-support source family with real files, rather than treating fixed/removable/implant as a footnote on one case.
- **Falsifier.** No public dataset supports mucosal support mechanics at the required resolution, in which case the class is recorded MISSING by acquisition, not by omission.

### 6. Graph hygiene: the two highest-scoring gaps are uncreatable MISSING_DRAFTs  `G06` — HIGH, graph_defect, stage L03/L15

`./graph working rank` returns DENT-PROC-STEP-MODEL (score 5.3) and DENT-UNC-PROPAGATION (score 5.2) as the top two items, both with admission_state MISSING_DRAFT and path_authority UNLINKED_GOAL_RELEVANCE. A planner cannot create the drafts it is being told to work, so it will re-derive a toy identity instead. `path_authority` is also UNKNOWN for both.

- **Smallest decisive next step.** Define the two missing drafts as explicit operator lists with their required inputs, so the graph stops ranking uncreatable nodes; and link them to DENT-GOAL so path_authority resolves.
- **Falsifier.** The existing FE-REFERENCE cell already satisfies DENT-PHYS-FIELD-SOLVE, in which case the node should be an alias, not a new draft.

### 7. No metrology loop: nothing measures an as-built part and updates the design  `G07` — HIGH, interface_break, stage L11

The system predicts as-built deviation to a few micrometres (tool radius 0.3-1.0 mm moves the marginal gap 12.3 um) and has never produced or scanned a part. K2 fit prediction carries a 38-77 um systematic anatomical residual whose two proposed causes are REFUTED. VAL-A is fully specified, hashad, and unrun.

- **Smallest decisive next step.** Mill and scan one crown on the license-free STS die and report the predicted-vs-measured as-built marginal gap. This closes design->manufacture->measure in one artefact.
- **Falsifier.** Measured deviation is dominated by tool wear/holding error rather than the modelled causes, so the process model is mis-attributing its own residuals.

### 8. No independent review binding anywhere in the system  `G08` — HIGH, evidence_gap, stage L14/L15

The interface audit found 1 independently reviewed node out of 102. Every report is PRODUCER_OUTPUTS_UNAUDITED. 371 audit + 387 adversarial jobs (54 % of producer work) re-check model arithmetic, yet 0 bind an independent verdict into the graph. The review effort is real and its output is not connected to anything.

- **Smallest decisive next step.** Bind the verdicts that already exist: take the completed REVIEW-* results, bind them as GRAPH_FEEDBACK with the existing content-addressed archive, and make the bound count a reported quantity in STATUS.json.
- **Falsifier.** The completed reviews were never independent in the sense the graph requires, in which case the audits are redundant with the producers and the 54 % allocation is wrong.

### 9. Wear, corrosion and aging are entirely unmodelled  `G09` — HIGH, missing_mechanism, stage L08

The durability layer models fatigue only. 73 fatigue jobs exist; zero wear, zero corrosion, zero aging/temperature-humidity. For a ceramic-restored chain, tribological wear of the antagonist pairing and hydrothermal degradation are first-order over service time.

- **Smallest decisive next step.** Acquire a three-body wear measurement on a declared antagonist pairing and fit a single Archard-type coefficient. The thermocycling reference protocol (5-55 C, 5 min/cycle, ~3000-10000 cycles, PMMA/poppy-seed third body) is the published route.
- **Falsifier.** The coefficient is pairing-specific and unstable across the three slurries, in which case the layer needs a wear mechanism rather than a wear coefficient.

### 10. Manufacturability is a post-check, not a search constraint  `G10` — MEDIUM, interface_break, stage L09/L10

Milling access, sinter shrinkage and cement spacer are quantified as deviations AFTER a design exists. Six MFG nodes carry the most precise numbers in the graph and none has a cell. A design that is geometrically optimal can be unmanufacturable, and nothing in the search knows that.

- **Smallest decisive next step.** Add accessibility and shrinkage as hard constraints in one existing design loop and measure the Pareto volume lost relative to post-hoc rejection.
- **Falsifier.** Post-hoc rejection removes nothing the constraint admits, so the coupling is unnecessary and the current architecture is already correct.

## Where the chain is disconnected even though both endpoints exist

- **FE solver -> crown design.** DENT-SOLVER-HIGHRES (PROVEN, 10.28 M DOF) vs DENT-PHYS-FE-REFERENCE (ASSUMED) -> DENT-DESIGN-J3-CREST-STRESS-SURROGATE (REFUTED). A verified solver and a design consumer both exist, but the field solution that would let the design ask the solver a different question is the MISSING_DRAFT DENT-PHYS-FIELD-SOLVE. The chain is broken at the middle, not at either end.
- **Bite load -> fatigue life.** DENT-LOAD-MAGNITUDE-SPECTRUM (ASSUMED, P5-P95 250-20000 /yr >500 N) vs DENT-K1-LIFETIME / 351 FATIGUE jobs. A 80x epistemic band on the load enters a life model whose own power law is post-hoc and non-independent. The compounding is never evaluated, so no fatigue number in the system is a number.
- **Measured crest band -> insertion torque.** DENT-PROC-CREST-BAND-MEASURED (PROVEN, 370 sites) vs DENT-PROC-INSERT-S4-TORQUE (PROVEN in foam only). The only measured tissue field in the system feeds a torque model that is refuted in the material that matters. A measured input and a measured output are both present and the link between them is absent.
- **Design -> manufactured geometry.** DENT-MFG-* six quantified nodes, all OPEN vs DENT-K2 fit prediction / DENT-VAL-BENCH-CROWN-FIT. The process model predicts an as-built deviation and nothing measures one. The consumer exists and the producer is a specification only.
- **Design variable -> biological outcome.** Every design variable (implant diameter, thread, surface, crown contour) vs DENT-BIO-* (3 REFUTED drivers). No admissible mechanism connects them, so the biological stage cannot currently veto or justify a design change.

## Wider backlog (not in the top ten)

- `B11` **Soft-tissue obstacle in the access volume (DENT-ACCESS-SOFT-TISSUE-OBSTACLE, OPEN, unranked)** (stage L10) — Manufacturing feasibility ignores lip/cheek/tongue. Milling access is modelled by tool radius only.
- `B12` **Imaging beyond CBCT** (stage L01) — No intraoral scan, optical or MRI input. Intraoral scans are the format the design loop actually wants.
- `B13` **Patient-specific manufacturing datasets** (stage L10) — A Ti-6Al-4V lattice mandibular scaffold study (15 fabricated, 384 computational variants) is dental-adjacent and CC BY 4.0; its reproducibility link is unresolved.
- `B14` **Fractography of clinically failed implant-supported crowns** (stage L08) — Failure-mode classification with dominant debonding at the zirconia-abutment interface. Converts the negative fatigue evidence into a design constraint.
- `B15` **Screw preload and superstructure assembly mechanics** (stage L05) — Assembly is absent from the model chain entirely; L11 covers fit, not the screw joint.
- `B16` **External validation protocol** (stage L16) — No held-out clinical site or multi-centre split exists; the system has no external validity claim.
- `B17` **Consent/de-identification gate in the controller** (stage L16) — Any job may read any dataset on disk. A machine-checked gate is cheap and closes a real privacy risk.
- `B18` **Hip and surface post-treatment route for L-PBF Ti-64** (stage L07/L08) — DENT-MATINNOV-H6 is OPEN and is the clearest available design freedom with a public test-bed.
- `B19` **Patent corpus as a design-variable source** (stage L09) — 107797 strict dental families are on disk with SureChEMBL links and were never mined for design variables; they are the cheapest source of genuinely new design freedoms.
- `B20` **Graph hygiene: resolve 7 duplicated draft-ID conflicts** (stage L15) — DENT-DESIGN-MARGIN-FROM-GEOMETRY-SPREAD, DENT-DS-MMDENTAL, DENT-DS-TOOTHFAIRY2 and 4 others are duplicated across drafts; they block packet freshness (STALE_INPUT).

## Prosthesis classes: what each one changes

### fixed

Tooth-borne. The abutment is prepared enamel/dentin with a controllable margin; the load path is crown -> luting cement -> prepared tooth. Needs: margin geometry, cement film, prep taper, and a tooth-side (not implant-side) fix.

Additional missing capability specific to this class:
  - cement film thickness and its effect on fracture load (2 of 3 transfers REFUTED)
  - preparation taper/height-of-contour as a design variable (absent)
  - thermal and chemical effects of luting on the tooth (absent)

Best-supported stage: L10 (MFG has the most precisely quantified, unrun statements)

### removable

Tissue-borne. No implant, no preparation. Support comes from mucosa and the alveolar ridge; the dominant physics is soft-tissue compression, not contact stress. The whole set of cavity/retention mechanics is absent.

Additional missing capability specific to this class:
  - no node, cell, or result of any kind addresses a removable prosthesis
  - mucosal support area and its pressure-relief behaviour (absent)
  - retention/clasp mechanics, wear of the retentive element, and the mucosa-to-hard-tissue transition (absent)
  - the implant-comparison framing is inapplicable: 0 of 1414 producer jobs is removable

Best-supported stage: NONE — the class is unrepresented

### implant_supported

Bone-borne, threaded, and the only class with a measured placement quantity. Needs implant position, bone quality, micromotion, osseointegration, and superstructure screw preload.

Additional missing capability specific to this class:
  - thread-into-bone constitutive law (REFUTED in cadaver bone, 8-13x overestimation)
  - superstructure screw preload and its relaxation (absent)
  - aesthetics/shade for the ceramic layer (absent)
  - remodeling driver after three REFUTED attempts

Best-supported stage: L05 and L12 (procedure transitions are the only genuine assumption-to-measurement closure in the system)

## Branches that should be paused

- **Further derivation of the inherited fatigue power law and the shell eigenstrain model.** 73 fatigue jobs and 141 grade jobs already exist on the same inherited identity. The next useful move is a measurement (T02), not a fifteenth derivation.
- **A fourth generic Frost mechanostat band.** Three are already REFUTED. One acquisition or an explicit empirical risk term, then stop (T03).
- **Further post-hoc identifiability and null-space analysis of the six-mode shell model.** The amplitude is unmeasured; the analysis is finished. T02 supplies the measurement.
- **Additional re-audits of inherited arithmetic.** 54 % of producer work is already audit or adversarial and zero reviews are bound into the graph. T07 binds what exists.

