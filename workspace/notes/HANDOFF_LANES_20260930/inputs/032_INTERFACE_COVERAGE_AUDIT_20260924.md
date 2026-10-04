# Tissue–construction interface coverage: read-only review 2026-09-24

Machine-readable version: `notes/INTERFACE_COVERAGE_AUDIT_20260924.json` (same numbers, all nodes with flags).

## 0. Basis and concepts

- Read: all 37 `notes/expansion/*.jsonl` (444 nodes with claims, plus 122 source records), `orient_B_graph.md` §6–§8, `AGENT_QUEUE.md`, and result directories (see §5). `GRAPH.json` contains 4 047 native nodes: 3 950 bodytwin and 97 manufacturing. None has a DENT prefix.
- `./graph working status|rank|packet` returned **STALE_INPUT** (draft_changed: patents_bq*.jsonl). `build` writes to disk and was therefore not run. The ranking in §4 is based on the skeleton §8 and AGENT_QUEUE, not on a fresh graph.
- The results distinguish four levels:
  1. **defined node**, that is, a row in expansion;
  2. **source-checked value**, that is, `cert_design.data_sources` contains DOI/PMID/PMC, a patent number or figshare;
  3. **tested cell**, that is, `evidence` points to a `results/` file;
  4. **reviewed result**, that is, an independent review bound to the node. Steelman analyses and reruns by the same lane do not count as review.
- Selection: 102 unique DENT nodes have been assessed as belonging to the interfaces (the list is in JSON). The selection is manual. A node is counted under the interface where it does most.

## 1. Key figures (102 unique interface nodes)

| Metric | Count | Share |
|---|---|---|
| Node defined | 102 | 100 % |
| `dent_kind` set (the contract requires the field) | 61 | 60 % (41 lack the field entirely) |
| `design_role` set | 91 | 89 % |
| Unit specified | 94 | 92 % |
| Validity (box or specific note) | 54 | 53 % |
| Uncertainty specified (kind ≠ UNKNOWN) | 60 | 59 % |
| Source in `data_sources` (DOI/PMID/patent) | 35 | 34 % (another 4 have PMID only in the claim text) |
| Test defined | 85 | 83 % |
| Cell run (evidence → results/) | 60 | 59 % |
| Independent review bound in the node | 1 | 1 % (DENT-PAT-CELLCAND-PDL-ROM ← BQ3_jobs/PDL_C3_AUDIT) |
| Status | OPEN 56 · ASSUMED 29 · REFUTED 12 · DEFERRED 3 · PROVEN 2 | |

The two PROVEN nodes are DENT-PROC-INSERT-S4-TORQUE and DENT-PROC-CREST-BAND-MEASURED. The first is proven **only in PU foam**: in cadaver bone (H6), the same chain overestimates insertion torque 8–13×, and the node DENT-PROC-INSERT-CADAVER is REFUTED.

There is no node combining source, uncertainty, a run cell and independent review. Source and run cell overlap poorly. The cell nodes (PROC, K1, K2, K3, CONTACT) usually have empty `data_sources` and carry their sources in the PREREG files. The literature nodes (gingiva, material, patent) have sources but no run cell. The PROC-INSERT nodes also share a common list of four PMID and have no source per number.

## 2. Table per interface

Columns: n = nodes; U = unit; S = source in data_sources; E = uncertainty; R = cell run; status; missing central variables (see §3). "Entirely missing" means that the variable does not exist as a number anywhere in the sources read.

| Interface | n | U | S | E | R | OPEN/ASS/PROV/REF/DEF | Missing central variables (of which entirely missing) |
|---|---|---|---|---|---|---|---|
| I1 bone–implant | 22 | 20 | 8 (+2 claim) | 14 | 19 | 3/11/2/5/1 | 10 (6) |
| I2 soft tissue–implant/abutment | 14 | 14 | 10 | 4 | 1 | 13/1/0/0/0 | 7 (5) |
| I3 tooth–PDL–bone | 13 | 12 | 4 | 7 | 5 | 8/4/0/0/1 | 5 (4) |
| I4 restoration–cement–tooth | 16 | 15 | 3 (+1) | 11 | 10 | 12/2/0/2/0 | 5 (4) |
| I5 surface → biological response | 6 | 6 | 6 | 3 | 0 | 6/0/0/0/0 | 5 (4) |
| I6 thermal (drilling) | 8 | 8 | 0 | 6 | 7 | 3/4/0/1/0 | 4 (2) |
| I7 contact/friction implant–abutment | 7 | 4 | 1 | 3 | 5 | 4/1/0/1/1 | 4 (2) |
| I8 loading | 9 | 9 | 0 (+1) | 6 | 9 | 4/5/0/0/0 | 3 (0) |
| I9 time (healing, remodeling, fatigue) | 7 | 6 | 3 | 6 | 4 | 3/1/0/3/0 | 3 (2) |
| **Total** | **102** | **94** | **35** | **60** | **60** | **56/29/2/12/3** | **46 (29)** |

### 2.1 Node level, selection (complete list in JSON)

**I1 bone–implant**
- DENT-PROC-INSERT-S4-TORQUE: u(t), Ncm, 4 PMID, uncertainty factor 1,42, cell run. **PROVEN (PU foam)**.
- DENT-PROC-INSERT-CADAVER: u(t), Ncm. **REFUTED**: overestimates 8–13×.
- DENT-PROC-INSERT-CORTICAL: **REFUTED**. Measured 2,71 against predicted 1,07.
- DENT-PROC-INSERT-MICROMOTION: J, µm. ASSUMED: the ordering matches 4/4, but the level is 0,29–0,66 of measured.
- DENT-PROC-INSERT-NECROSIS: g, MPa, the boundary 133 MPa according to Bashutski 2009, PMID 19335092. ASSUMED; the threshold lacks preregistration.
- DENT-DESIGN-ANCHOR-PRIMARY-STABILITY: J, N/µm. ASSUMED; 43 % of decisions change.
- DENT-DESIGN-ABSOLUTE-ANCHORAGE-REQUIREMENT: g, µm. OPEN, test UNKNOWN. The threshold 50–150 µm (PMID 9619438) is only in the claim text.
- DENT-DESIGN-J2-BONE-CONTACT-PROXY: J, mm². ASSUMED. The node's honest_gap reads: "no primary stability/osseointegration model".
- DENT-DESIGN-J3-CREST-STRESS-SURROGATE: **REFUTED**. DENT-DESIGN-J3-FE: ASSUMED.
- DENT-MAT-BONE-HU-MODULUS and -TRABECULAR-E: **REFUTED**. -CORTICAL-E: ASSUMED. The bone material is therefore run as a scenario (DENT-DESIGN-THETA-BONE-SCENARIO, ASSUMED).
- DENT-MATINNOV-H2-TI-LATTICE-BONE-MATCHED: OPEN hypothesis. Stress shielding is included as a premise, and the alveolar mechanostat thresholds are UNKNOWN.

**I2 soft tissue–implant/abutment**
- The gingiva (DENT-MAT-GINGIVA-ATTACHED, E 19,75 ± 6,20 MPa according to PMID 21309720) and the phenotype nodes are source-reported but OPEN.
- DENT-GEOM-GINGIVA-THICKNESS-CBCT is the only run cell (SOFT_gingiva): H1 UNDETERMINED, H1d FAIL, H2d FAIL, H3 PASS 4/5, H4 PASS.
- DENT-DESIGN-G-SOFT-THIN-PHENOTYPE and DENT-DESIGN-J-ESTHETIC-RECESSION are hypotheses without tests.
- **There is no mucosa–implant/abutment interface region.** Anatomy has 11 IF regions, and none concerns soft tissue against a construction.

**I3 tooth–PDL–bone**
- DENT-MAT-PDL-CONSTITUTIVE: θ, MPa, interval 0,52–1,65, SIMULATED. ASSUMED and valid for 5–16 N.
- DENT-MAT-PDL-NONLINEAR: ASSUMED. The curve shape fails (2/6 points).
- DENT-PHYS-FE-LINEAR-VALIDITY-LIMIT: at 843 N, PDL strain exceeds 10 % in all linear scenarios.
- DENT-REG-PDL: a primary source for PDL width is missing, and OFJ's 0,2 mm is constructed.
- The IF regions ROOT-PDL, PDL-BONE and CORT-CANC have neither source nor uncertainty.
- DENT-VAL-BENCH-PDL-MOBILITY: DEFERRED.

**I4 restoration–cement–tooth**
- DENT-PHYS-CONTACT-DAMAGE and -CONTACT-DIES: **REFUTED**. The ratio against the best baseline is 1,03, and the gate required ≤ 0,8.
- DENT-MFG-MARGINAL-GAP-LIMIT: g, 120 µm (McLean 1971), OPEN.
- DENT-MFG-CEMENT-SPACER: x, µm, OPEN. Movement in the marginal gap is ≤ 5,9 µm.
- DENT-REG-CONSTRUCT-CEMENT-LAYER: region without source and uncertainty, with the note "konstruktion ej vald".
- DENT-MAT-DEJ-GRADED (tissue–tissue) is source-reported (PMID 16414277, 11077406) but has no cell.

**I5 surface → biological response**
- All six nodes are OPEN and have no run cell. The four patent nodes have patent numbers as sources, but the values are distributions over claimed intervals (roughness 0,7–8 µm, median 2), so not measurements and not biological response.
- The surface → mechanics connection exists in K1 (DENT-K1-SIGMAW-ROUTE REFUTED) and SURF_strength. The surface → biology connection is entirely missing.

**I6 thermal (drilling)**
- DENT-PROC-DRILL-THERMAL: **REFUTED**. The model is approximately 3,5× too low, and only 3/16 conditions fall within a factor of 2.
- The criterion 47 °C/1 min = CEM43 16 min (PMID 6576145, 6547421) is in `results/K3_drilling/PREREG.md` but not in the node's `data_sources`, which is empty for all 8 nodes.

**I7 implant–abutment**
- DENT-K1B-FLIM: **REFUTED**. pred/meas is 0,61, and 11 % are within the band.
- The cone angle exists only as a patent distribution (3,5–15°).
- Friction μ is varied 0,12–0,50 in K1b's Monte Carlo but has no node of its own.
- DENT-PAT-CELLCAND-TAPER-SELFLOCK: DEFERRED.

**I8 loading**
- All 9 nodes have run cells, but none has a source in `data_sources`. The sources are in PREREG and claim.
- Hard-bite frequency and the per-tooth share are OPEN and have the value UNKNOWN.

**I9 time**
- Fatigue exists (K1 and SCG) but is **REFUTED**: DENT-K1-ISO14801-MARGIN, DENT-K1-SIGMAW-ROUTE and DENT-MAT-SCG-LIFE-CROWN.
- Healing and remodeling have no DENT node.

## 3. Central variables that are missing

Labeling: **UNKNOWN** means that no value exists in the project sources read. If a number exists, it is stated together with the source and the place where it appears.

**I1 bone–implant (10 variables, of which 6 entirely missing)**
- BIC fraction (%): UNKNOWN. "BIC" in SEG_cbct means Bayesian information criterion.
- Micromotion threshold: 50–150 µm (Szmukler-Moncler 1998, PMID 9619438). The number is only in claim text and in MICRO_motion/PREREG and lacks a g node of its own.
- Osseointegration rate and secondary stability ISQ(t): UNKNOWN. PROC_insert/REPORT states that the ISQ proxy has not been calculated.
- Clinical window for insertion torque: UNKNOWN.
- Alveolar mechanostat threshold: UNKNOWN. The Frost bands are only in the generic bodytwin cell `results/I1/out/bone_wolff_law_mechanostat`.
- Stress shielding metric: UNKNOWN, only a premise.
- Crestal bone loss: only "upp till 1,45 mm" (Linkevicius 2009, PMID 19885413) in a hypothesis node's claim.
- μ bone–implant: 0,40 (Grant 2007, PMID 16806236), only in MICRO_motion/PREREG.
- Bone–implant contact law after healing: UNKNOWN.
- Healing and loading protocol: UNKNOWN.

**I2 soft tissue–implant (7 variables, of which 5 entirely missing)**
- Biological width at implants: UNKNOWN. What exists is CEJ–bone crest at a tooth, 1,6–3 mm (PMID 21320168).
- Mucosa–abutment interface region: UNKNOWN.
- Transmucosal seal: UNKNOWN.
- Abutment material/surface as x: UNKNOWN.
- Risk of peri-implantitis or mucositis: UNKNOWN.
- Vertical mucosal thickness: exists as the hypothesis's threshold ≤ 2 mm but is not measured.
- Recession model: UNKNOWN.

**I3 tooth–PDL–bone (5 variables, of which 4 entirely missing)**
- PDL thickness from a primary source: UNKNOWN.
- PDL viscoelasticity: UNKNOWN.
- ν_PDL: UNKNOWN.
- PDL and alveolar remodeling: UNKNOWN.
- Confinement in a thin layer: run in BQ3_jobs/PDL_CONFINEMENT. There is no node and no evidence.json.

**I4 restoration–cement–tooth (5 variables, of which 4 entirely missing)**
- Bond strength (MPa): UNKNOWN.
- CTE mismatch: UNKNOWN.
- Cement modulus: **conflict**. The node has 7,7 GPa from one source, while BQ3_jobs/CEMENT_LAYER_MECH found 12,74–13,63 GPa (PMID 37068429). The latter value has not been loaded into any node.
- Cement aging and leakage: UNKNOWN.
- Biological outcome at a gap: UNKNOWN.

**I5 surface → biological response (5 variables, of which 4 entirely missing)**
- Dose–response between roughness and BIC: UNKNOWN.
- Surface chemistry: UNKNOWN.
- Effect of coating: UNKNOWN.
- Biofilm: UNKNOWN.
- Measured roughness: exists only as a distribution over patent claims.

**I6 thermal (4 variables, of which 2 entirely missing)**
- 47 °C/1 min: exists in PREREG but not in the node.
- η and χ: UNKNOWN.
- Validated temperature field: refuted.
- Thermal properties of a tooth: UNKNOWN.

**I7 implant–abutment (4 variables, of which 2 entirely missing)**
- μ: exists as the interval 0,12–0,50 in K1b.
- Preload loss: DEFERRED.
- Microgap in the coupling: UNKNOWN.
- Bone–implant contact law in J3-FE: UNKNOWN.

**I8 loading (3 variables, none entirely missing; all nodes exist but have the value UNKNOWN)**
- Hard-bite frequency.
- Per-tooth share of chewing force.
- Rise time and damping for an implant compared with a tooth.

**I9 time (3 variables, of which 2 entirely missing)**
- Osseointegration time: UNKNOWN.
- Dental remodeling cell: UNKNOWN. Bodytwin has MODEL-BONE-WOLFF-LAW-TRABECULAR, MSK-BONE-REMODELING-MECHANOSTAT and MODEL-BONE-REMODELING-RANKL-OPG, all OPEN, but no DENT node links to them.
- Soft tissue healing: exists only as a hypothesis.

**Total:** 46 central variables are missing or have no source-backed value of their own. Of these, 29 are entirely missing.

## 4. The chain geometry → material → load → response → biological outcome

| Interface | Geometry | Material | Load | Response | Biological outcome | Break |
|---|---|---|---|---|---|---|
| I1 bone–implant | yes (TF2, NV1, BONE_PV) | scenario (HU→E REFUTED) | yes | partial (J3 surr. REFUTED; micromotion level 0,29–0,66; MICRO_motion H1 FAIL ×5,9) | **missing** | material, response→biology |
| I2 soft tissue–implant | partial (at a tooth) | source-reported | **missing** | **missing** | hypothesis | region, load, response |
| I3 tooth–PDL–bone | partial (PDL constructed) | ASSUMED | yes | yes for 5–16 N; breaks at 843 N | **missing** | response→biology |
| I4 crown–cement–tooth | yes | one source + conflict | yes | **refuted** | **missing** | response, biology |
| I6 thermal | yes | η, χ UNKNOWN | yes (torque) | **refuted** (×3,5 too low) | criterion defined | response |
| I7 implant–abutment | yes | yes | yes (ISO 14801) | **refuted** (0,61) | not modeled | response |

**No chain is closed all the way to biological outcome.** Drilling (I6) is closest, because the outcome criterion (47 °C/CEM43) is defined, but the response step is refuted. All mechanical chains reaching response with external testing have failed, except PDL in the low-load regime and insertion torque in PU foam.

**Load-bearing gaps toward DENT-GOAL (skeleton §8.2):**
- `DENT-PHYS-INTERFACES` (F role, depends on FIELD-REP and MAT-RESTORATIVE-PSP, and DESIGN-PROBLEM in turn depends on it) **does not exist as a node**. It is referenced by 5 expansion files but has no row of its own.
- The same applies to DENT-PHYS-FIELD-SOLVE, DENT-PROC-STEP-MODEL, DENT-GEOM-ANCHOR-PDL-BONE, DENT-PROC-SOFT-TISSUE-RESPONSE and DENT-GOAL: they are referenced but have no row of their own.
- Interface physics thus reaches the design problem only through ad hoc assumptions in the cells (Winkler foundation, spring bed for cement, bonded contact). It does not go through a defined node.
- DENT-PROC-SOFT-TISSUE-RESPONSE has reach 0 according to §8.2. The soft tissue interface is therefore a coverage branch, not a load-bearing branch, as long as the skeleton does not change.

## 5. Results without a node or unreviewed (negatives visible)

| Directory | Node | Content | Review |
|---|---|---|---|
| results/MICRO_motion/ | none (PREREG calls it DENT-PROC-MICROMOTION) | H1 **FAIL**: pred/obs 0,13 and 0,22, underestimation ×5,85; H1b PASS. No REPORT.md or evidence.json. Amendments A5–A6 were made after H1 was read | none independent (own steelman) |
| BQ3_jobs/CEMENT_LAYER_MECH | none | P1 **FAIL**, P2 PASS, P3 PASS, P4 **FAIL** | none |
| BQ3_jobs/PDL_FIBROUS | none | Fung-ROM: linearity breaks at u_c/t ≈ 0,019 (C1–C4 PASS according to evidence) | none |
| BQ3_jobs/PDL_C3_AUDIT | DENT-PAT-CELLCAND-PDL-ROM | Independent review of PC1: C1 is a tautology | is the review |
| BQ3_jobs/PDL_CONFINEMENT | none | ongoing (run.jsonl 106 rows, no evidence) | – |
| results/DESIGN_margin/ | none | only PREREG (frozen 10:52 today) | – |
| results/I1/out/bone_wolff_law_mechanostat | bodytwin, not DENT | G1–G4 PASS; generic | – |

Own steelman or rerun is mentioned in 10 interface nodes (K2, SCG, LOAD, BONE_HU). The folder DS_CROWN_PROC_AUDIT_20260923 is an independent review of DS_PROC_MECHANISM and DS_CROWN_ANCHOR but is not bound to any of the 102 nodes.

## 6. Top 5 next steps

The order is based on reach toward DENT-GOAL and cost. It is heuristic, because the graph view is STALE.

1. **Define `DENT-PHYS-INTERFACES` with four children** (bone–implant, implant–abutment, cement–tooth, mucosa–abutment), with dent_kind, units, validity-box and data_sources per number.
   - Reach: the skeleton's F node on the path to DESIGN-PROBLEM and GOAL.
   - Cost: low.
   - *Test that can refute:* the contract check. A child that lacks a unit, validity or source per number after a source search must be UNKNOWN and must not be included in `parameterizes`. If MICRO_motion and K1b with explicit interface compliance still do not move the ratio into [0,5; 2], then the interface is not the explanation for the underestimation.
2. **Bind MICRO_motion as a node (H1 FAIL visible) and obtain an anchor measuring interface slip.** The current anchor measures the abutment's displacement.
   - Reach: DENT-DESIGN-ABSOLUTE-ANCHORAGE-REQUIREMENT, which according to AGENT_QUEUE is load-bearing number 3 with 7 859 paths.
   - Cost: medium or EXTERN.
   - *Test:* |ln(pred/obs)| for slip in bone must be ≤ ln 2. If it fails, g_abs may only be used relatively.
3. **Retrieve and source-check the outcome variables** BIC fraction, ISQ(t), crestal bone loss (mm/year) and alveolar mechanostat threshold, and link the bodytwin nodes for mechanostat and RANKL/OPG.
   - Effect: closes the step from response to biology in I1.
   - Cost: medium.
   - *Test:* a source-backed relation (micromotion or strain → BIC or bone loss) must predict the direction in a held-out study. The opposite direction refutes it.
4. **Resolve the conflict over cement modulus** (7,7 against 12,7–13,6 GPa) and introduce bond strength and CTE as θ. Then rerun the sensitivity analysis in CONTACT_damage.
   - Cost: low.
   - *Test:* if the ratio RMSE(CD)/best baseline remains above 0,8 with both moduli, the cement modulus is not the reason H1 failed.
5. **Test DENT-HYP-K3-MISSING-HEAT-SOURCE** by calculation, and then bench D.
   - Reach: I6 through PROC-STEP-MODEL.
   - Cost: low for the calculation, EXTERN for the bench.
   - *Test:* at least 60 % of the held-out conditions must fall within a factor of 2 (same gate as K3). Otherwise the hypothesis is refuted.

Outside the top 5: interface region and biological width for mucosa–abutment (reach 0 in the skeleton), surface → biology (no cell) and healing/remodeling time.
