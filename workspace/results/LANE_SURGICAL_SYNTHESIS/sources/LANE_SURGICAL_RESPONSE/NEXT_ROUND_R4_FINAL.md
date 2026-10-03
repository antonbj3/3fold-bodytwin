# Next round after RESPONSE R4

PENDING_INDEPENDENT_REVIEW. Round4 is complete: joint prediction/innovation gate FAIL, equally informed control TIE, native joint UNKNOWN. Parent target GRAPH_INVERSE/MULTIPHYSICS and BT-CTX-SURG-HEALING/HEMOSTASIS are OPEN. Read new steer first. All previous failures and raw data are preserved.

Capability → the same incision history should predict chemical inventory and directional rupture load to90days, as well as new load/O2/design trajectories from a sufficient local state.
Conflict → HP increases while collagen amount barely decreases; early chemistry is fast but strength recovery slow; different glycosylation branches do not share yield or timescale.
Obstacle → absolute native precursor/glyco pool and local connected bridge state are missing; two external assay endpoints do not identify the first- versus second-order law. Normal skin day0 is no observed initial wound state.
Changed operation → measure the same skin cohort's stoichiometric precursor→HP/other-product flux under LOX blockade and load variation, and feed these rates and connected-fiber history to a shared chemistry/mechanics checkpoint.

## Start from state that was actually run

Read CHECKPOINT_R4.json, SOURCE_REVIEW_R4.md, DERIVATION_R4.md, OBSERVATIONS_R4.json, PREREG_R4_GLYCO_RESERVOIR.json, r4/SCORES_R4_v2.json, r4/glyco_reservoir/summary.json and RESPONSE_PORTS_R4.json. Check peer ports read-only. R3's Q_U/Q_I/Q_M and R1's spatial O2/cell states are preserved but R4 has not changed them.

The new nominal isolated chemistry checkpoint is `r4/glyco_reservoir/control_first_order_D_NS0.03_checkpoint21.json`. It is synthetic, not a biological posterior. `D+2H+A` is exactly invariant in the declared autonomous closure. `evolve(D,H,A,dt,rate,eta,law)` in glyco_reservoir_r4.py gives the future without history replay; the same conventional kinetics operator is TIE. For nonlinear spatial O2/load feedback a new contract is required, no automatic global reuse-PASS.

## Concrete next measurement/construction change

1. Acquire absolute **GG-DHLNL, G-DHLNL, free DHLNL, HP, alternative aldehyde product, collagen moles per constant reference area/volume** at10/14/21/42/61/90 from the same skin-incision model. Record the sample's changed wound geometry and amount of normal skin: a shifting whole-wound sample can imitate inventory change. Hyp/hydroxyproline and HP/hydroxylysylpyridinoline are distinct observables. Separate destroyed samples must have cohort/donor linkage; no invented joint-SD.
2. Freeze a LOX blockade from10, with verified catalytic inhibition and remaining precursor. The BAPN cell-culture study gives **protocol and analogue**, not skin rate. The available native-D10 must suffice for new HP: the lowest readout box requires .01494 HP moles per C10-ref, thus .02989 divalent-equivalent under2:1 stoichiometry if everything becomes HP. Alternative outcomes increase the precursor requirement. Import/loss/sampling are separate ledger flows.
3. Acquisition start: Wietecha ETH bundle DOI10.3929/ethz-b-000409545 for absolute crosslink data if available; publisher source data is missing from the link menu. Latest handle500, previous429, BioC captcha and EuropePMC500 are access errors. Sricholpech2012 DOI10.1074/jbc.M112.343954 has Figure4 at0/14/28days: its **independent day14** can distinguish the rate laws before new skin prediction. Acquire the original image/full text and freeze dates/rates before fit; the endpoint means used are explicitly approximate. No author contact is included.
4. Implement D_GG,D_G,D_free,H,A with shared birth/removal fluxes and observed catalyticLOX. Do not present GG loss as HP flux. Loss and HP yield must have the same assay lineage. First versus second order were both exact on external endpoints but gave different skin curves; more net does not solve this. Freeze genotype/load/late-time holdout and strongest conventional multitype-renewal model with the same information.
5. Measure connected collagen crossings across the actual cut front, terminal anchor/pullout survival and traction as a function of age/chemistry. SKIN_TOUGHNESS_GAP PORTS_R2_FINAL.native_healing_traction and the other biological damage ports are still null. Strength is not Γ, SHG eccentricity or low-strain modulus. Wietecha **has** Fig8k ultimate load day21; reuse it, correct R3's incorrect deformability-only description. It is not the same biochemical sample and normal skin comes from an earlier cohort.
6. Do load0/1 perturbation after21 with chemistry×angle×bridge state and ledgers. Autonomous post-initiation chemistry should persist at LOX stop and be load-independent within its contract; mechanofeedback should give an observable load response. Couple measured rates to the remaining local O2 state only when reference volume, cell/vascular damage port and feedback are measured. Charge acquisition/setup/history/storage/cert/rebuild/fallback and compare with the same conventional cohort/cohesive operator. No such global gain exists now.

## Preserved falsifiers and limits

HP inventory42/10≈2,220, box1,747–2,847. Denominator-only selective removal requires massratio.356–.488 but observedbox.852–1.012: the specified model fails under a fixed assay reference, not all turnover.

Early mechanistic fits pull rates to bounds and miss HP by at least≈.00475/.00602 mol/mol. Reservoir nominal HP42=.03687 passes20% but is fast k=3/dag, misses early chemistry and worsens strength-RMSE17,77 againstR3native8,31. No support for identified late maturation there.

Cell-free analogue: control k=.06950/dag, η=.2833; LH3-knockdown k=.04474/dag,η=.92. Native D_NS is unknown. For the same syntheticD_NS=.03, firstorder LH3 gives HP42=.04111 but strength-RMSE11,36 and earlyday21 conflict. This is not a selected winner or native genotype transfer. All16second-operation joint gates FAIL. The inverseD_NS requirements.03246–.13525 are diagnostic fitted requirements, never forward inputs. No parameter was selected from late chemistry or strength.

R4's explicitly synthetic defaults: stop initiation/loss at10 in the first reservoir test, feedbackcapacity.2mol/mol, b=.001, f=.5 andd=.03/dag for selective continuation, glycoD_NS=.03 and independentcontrol/firstorder as nominal.75%strength reference and(.043−.003)=.040scale are reused fromR3. HP42 is already known: no blind prediction claims.126/126checks concern mathematics/causality/provenance, biologyUNKNOWN.

## Reproduktion

COMMANDS_R4.md states the order and new exclusive destinations. Do not run original scoring/plot/verify/finalize again on fixed files: copy lane inputs/code to a new lane-local replica or give new suffixes. mechanisms_r4.py requires new --out; score_r4.py uses controlling predictions_v2 and v2-summary. Initial empty dense-output slice, the string-threshold error and initial score metadata are preserved.2CPUthreads, no subagent/cloud/queue/graph/product/peer write. GRAPH_FEEDBACK_R4.json is local coordinator material, no admission.
