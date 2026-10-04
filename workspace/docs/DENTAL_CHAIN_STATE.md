# Dental restoration: fit, tissue response and service life

Can patient anatomy and bite load determine a restoration's fit, tissue response and service life?
The available calculations do not close that path: insertion transfer to cadaver bone fails, added mechanics loses to simpler remodeling controls, and manufactured fit and lifetime lack matched measurements. [1, 2, 4, 7]

Tooth preparation and crown contact describe parts of a fixed restoration, while insertion and bone response describe implant-supported cases. Removable support and retention mechanics are missing from the coverage inventory. Shared geometry does not make those support laws interchangeable. [1]

## Insertion torque: foam agreement does not transfer to bone

Insertion error is the median absolute natural-log prediction/observation error, with lower values better. The foam result meets its recorded magnitude criterion but loses to the geometric-area baseline. Its cadaver transfer fails. [4]

| Quantity or test | Chain result | Comparator or scope | What it establishes |
|---|---|---|---|
| Insertion torque in rigid foam | Error 0.3517 across 44 conditions | Geometric area baseline: 0.3195 | Magnitude criterion met; error exceeds the baseline. [4] |
| Insertion torque transferred to cadaver bone | 1 of 4 points within a factor of 2 | Published density-dependent observations | Transfer fails; foam agreement does not establish bone behavior. [4] |

## Bone loss and tooth movement

The bone-loss levels compare calculations with extracted clinical group means. They do not supply patient-specific bite forces or specimen-matched material laws. Matching the overall level also does not establish the contrasts between clinical groups. [5, 6]

| Quantity or test | Chain result | Comparator or scope | What it establishes |
|---|---|---|---|
| Marginal bone loss after 5 years from loading | Mechanics plus biology: 0.5557 mm; biology alone: 0.2928 mm | Extracted pooled level 0.49 mm; system means span 0.24–0.75 mm; fixed baseline 2.3 mm | Both pass the recorded level criterion and beat that baseline; the span is not a confidence interval. [5] |
| Remodeling contrast groups | Mechanics plus biology: 1/5; biology alone: 2/5 | Recorded primary contrast groups | Both branches remain unresolved despite passing the level test. [5] |
| Thin-minus-thick mucosa contrast | Biology alone: +0.2747 mm | Extracted +0.80 mm, interval +0.42 to +1.18 mm | Magnitude comparison fails. [5, 6] |
| Platform mismatch minus matched platform | Biology alone: −0.2663 mm | Extracted −0.37 mm, interval −0.55 to −0.20 mm | Recorded comparison passes; it does not rescue the other contrasts. [5, 6] |
| Load-proxy remodeling test | 1/3 contrasts pass | Constant-loss baseline: 3/3 | The mechanics modifier is rejected under the recorded rule. [7] |
| Pressure-driven tooth-movement rate | 5/8 contrasts pass | Constant speed: 5/8; uniform-stress law: 5/8 | Finite-element pressure adds no score improvement in this test. [8] |

Remodeling contrast scores use the recorded sign and magnitude rules. The thin-mucosa contrast fails, even though the platform contrast meets its comparison rule. Both remodeling branches remain unresolved. [5, 6]

The load-proxy observations do not establish a statistically significant direction. Constant loss beating the mechanics modifier therefore does not prove that load has no effect. [7]

The tooth-movement score counts contrasts inside extracted acceptance bands. Those include an animal direction observation and reconstructed bands, so the tie with simpler controls is specific to those conventions. [8]

## Biology can change a design without establishing benefit

A biological term changes some feasible design choices under maximum platform mismatch. That coupling exists, contrary to the coverage narrative's claim that biology cannot affect a decision, but its biological justification remains unresolved. [1, 5, 9]

| Quantity or test | Chain result | Comparator or scope | What it establishes |
|---|---|---|---|
| Biology term added to implant design | Choices change at 0% with matched platforms and 2.7% with maximum platform mismatch | 37 jointly feasible sites; assumed mucosal thickness 3 mm | Conditional design effect; no measured benefit. [9] |

## Anatomy and the missing physical interfaces

Differences between canal labels across dataset releases measure release sensitivity, not error against anatomical truth. The comparison has no independently attributed annotation pairs, so it cannot establish an anatomical error floor. [2, 3]

| Quantity or test | Chain result | Comparator or scope | What it establishes |
|---|---|---|---|
| Proposed gap values | 25 proposed values touching 7 of 57 previously unidentified edges; 50 remain unfilled | 8 values explicitly marked unmeasured | Proposals, not completed measurements or closed edges. [2] |
| Located canal distance between releases | Mean absolute change 0.1706 mm; maximum old-minus-new change 5.9362 mm | 397 paired sites from 2581 candidate sites; independently attributed annotation pairs: 0 | Observed label/release changes; anatomical error floor remains unknown. [3] |

The gap values are proposals and are not applied to the dental network. Review of a calculation does not supply independent physical evidence for a complete restoration. Coverage statements describe the supplied inventory, not every possible implementation. [1, 2]

No patient-specific lifetime, healing time or safe procedure setting follows from these results. Anatomy, regional bite force, bone mechanics, evolving contact and longitudinal outcome are not paired for the same subject. [1, 2]

Preparation and drilling fail temperature-magnitude comparisons. Bone-stiffness and micromotion computations also retain failed transfers. Clinical survival alone cannot identify a specimen's fracture mechanism. [1, 2]

The evidence lacks a measured seated cement film, an independently scanned manufactured crown and a measured relation between fabrication process and material properties. Synthetic metrology and calculated gaps do not establish manufactured fit or attribute its error to holding, tool wear or seating. [1, 2]

Wear, corrosion, aging, removable mucosal support and time-dependent screw-preload loss lack the required evidence across the restoration process. The shared tissue constraints supply no dental-specific anatomy–load–remodeling relation that closes these interfaces. [1, 2, 10]

## Where the numbers come from

The mappings retain the displayed precision and the stored source values. Paths beginning with
`dental/` refer to the sibling workspace from the common workspace root. Other paths are relative
to this repository.

1. Coverage and scope: `notes/HANDOFF_LANES_20260930/inputs/029_DENTAL_CHAIN_COVERAGE_20260928.json :: lifecycle_rows[id=L04,L07,L08,L10,L11,L13,L14,L15,L16], prosthesis_class_deltas, not_claimed`; `tasks/free48/sources/BIORESP/dental_physics_material_procedure.md :: sections 1–3`.
   Narrative context: `notes/HANDOFF_LANES_20260930/DENTAL.md :: Existing material, Automatic work and larger gaps`; `notes/HANDOFF_LANES_20260930/inputs/021_DENTAL_CHAIN_HANDOFF.md :: coverage requirements`; `notes/HANDOFF_LANES_20260930/inputs/030_DENTAL_CHAIN_GAPS_RANKED_20260928.md :: gaps G01,G02,G03,G05,G07,G09,G10`. Their proposed next steps are not execution receipts.
2. Proposal accounting: `results/ASSEMBLY_DENTAL_GAP_PROPOSALS/PROPOSALS_V1.json :: _meta.counts.proposals=25, edges_filled=7, unidentified_edges=57, edges_unfilled=50, proposals_measured_false=8`.
   Meaning and limits: same file `:: _meta.status, _meta.method, _meta.reviewed_vs_unreviewed_note, proposals[*].measured, proposals[*].source_self_description, _meta.unfilled_edges[edge_id=D-E-K12,D-E-K26,D-E-K31,D-E-K32,D-E-K35,D-E-K43,D-E-K48]`. “Filled” is the search's proposed coverage label, not physical closure.
3. Release differences: `dental/results/LANE_X58_CANAL_WALL_SPREAD/results.json :: R3.mean_located_abs_change_mm=0.17057548353076715, R3.max_located_old_minus_new_mm=5.936235214522559, R3.paired_sites=397, R3.x8_total_sites=2581, dropout.R1.independent_pairs_admitted=0`.
   Same file `:: external_referent, outcome, dropout.R3` identifies release revisions and missing pairing, not independently measured anatomy.
   The retained `dental/results/LANE_X58_CANAL_WALL_SPREAD/results.json` locator is absent. The source is available at `../dental/results/LANE_X58_CANAL_WALL_SPREAD/results.json` with the same keys and values.
4. Torque: `dental/results/PROC_insert/analysis.json :: H1.median_abs_ln=0.3517008361415458, H1.n=44, H2.B0area.median_abs_ln=0.3195358972281164, H6.n_within2=1, H6.points` (4 entries); `dental/cells/procedure/analyze_insert.py :: LN2, H6.within` defines the factor 2 comparison.
5. Bone-loss output and rules: `dental/results/REMODEL_crestal/verdicts.json :: legs.L1_full.per_anchor.A1.5.pred=0.5556738510648823, legs.L2_bio.per_anchor.A1.5.pred=0.2927571843982158, legs.L1_full.per_anchor.A1.5.obs=0.49, legs.L1_full.per_anchor.A1.5.baseline=2.3`.
   Same file `:: legs.L1_full.n_contrast_pass=1, legs.L2_bio.n_contrast_pass=2, legs.*.n_contrast=5, legs.*.overall, legs.L2_bio.per_anchor.A2.pred=0.2746500257046316, legs.L2_bio.per_anchor.A3.pred=-0.26632729765297647`; `dental/results/REMODEL_crestal/compare_anchors.py :: contrast_rule, groups, overall`.
6. Extracted clinical anchors: `dental/results/REMODEL_crestal/anchors_heldout.json :: A1.obs_by_t.5=0.49, A1.per_system_5y_from_loading` (range 0.24–0.75), `A2.obs=0.8, A2.ci=[0.42,1.18], A3.obs=-0.37, A3.ci=[-0.55,-0.2]`; extraction and DOI fields distinguish the clinical source from the chain's output.
7. Proxy test: `dental/results/MECHANO_sign/verdicts.json :: by_scen.primary.n_L1_pass=1, by_scen.primary.n_baseline_pass=3, by_scen.primary.rows` (3 scored contrasts), `by_scen.primary.n_sig, verdict`; `dental/results/MECHANO_sign/compare.py :: s_o, keep`; `tasks/free48/sources/BIORESP/MECHANO_sign_REPORT.md :: Disclaimer, Next point of stress`.
8. Orthodontic test: `dental/results/ORTHO_rate/ortho_summary.json :: H1.S_M=5, H1.n_contrasts=8, score.B0=5, score.B2=5, H1.pass, H3_FE_vs_uniform`; `dental/results/ORTHO_rate/analyze_ortho.py :: evaluate`; `tasks/free48/sources/BIORESP/ORTHO_rate_REPORT.md :: Rerun and gaps` records animal anchors and band reconstruction.
9. Conditional design: `dental/results/REMODEL_crestal/design_mbl.json :: variants.L2_matched_tmuc3.0.R_changed_frac=0.0, variants.L2_maxPS_tmuc3.0.R_changed_frac=0.02702702702702703, variants.L2_maxPS_tmuc3.0.n_both=37`; `dental/cells/design/design_implant_v2_mbl.py :: main, t_muc, key` defines the assumed thickness 3 mm encoded by the variant key.
   `tasks/free48/sources/BIORESP/REMODEL_crestal_REPORT.md :: Design coupling, Steelman` supplies the conditional interpretation; `tasks/free48/sources/BIORESP/analyze_remodel.py :: fromload` defines the loading origin separately from insertion.
10. Constraint scope: `data/CONSTRAINT_NET_TISSUE.json :: bodytwin.tissue_constraint_net.edges[*].between`; no dental-specific chain closure is supplied by those endpoint definitions.

Status: pending independent review. No claim of clinical validation.
