# BT-HX-Q036 — preregistration

## Scope and question

This is a first runnable mechanistic model for the core question: which states accumulated during ischemia are needed to predict the early response after reperfusion? The modeled tissue is a regional skeletal-muscle compartment, with normalized flow, oxygen availability, energy state, acidosis/lactate, calcium, ROS, mitochondrial permeability transition (mPTP), edema, irreversible injury, and contractile function. It is not a clinical prediction model and it is not fitted to internal data.

The primary output is the state of contractile function and its mechanistic intermediates 15 min after reperfusion, following a 90 min low-flow ischemic challenge. The model also reports flow, oxygen, ATP, lactate, acidosis, calcium, ROS, mPTP opening, edema, injury, and function throughout the pre-ischemic, ischemic, and early reperfusion phases.

## Builds on

- `node-id: Q036` — the bounded question in `inputs/QUESTION.md`: regional flow, metabolism, and post-reperfusion function. The named input capability identifiers `node-id: K01`, `node-id: K03`, `node-id: K09`, and `node-id: K12` are not present as files or internal records in this directory, so no node-specific semantics are invented.
- `inputs/NIGHT_PREAMBLE.md`: requirement to work only in this directory, distinguish source/derivation/hypothesis, use public literature, and avoid invented measurements.
- `agent.log`: the interrupted session only inspected the directory and inputs; no earlier model, preregistration, or result file was present.
- `ALLOW_WEB`: public literature lookup is permitted.

## Not redone


## Verified primary anchors

These are literature anchors, not calibration data. The values were looked up in the primary articles; none is marked `UNVERIFIED`.

1. El Baradie et al., *Scientific Reports* 11, 6152 (2021), DOI `10.1038/s41598-021-85753-x`, PMID 33731782, PMCID PMC7969970. In the mouse single-hindlimb ischemia experiment, ischemia was 90 min. The Results text and Fig. 2–3 report approximately 1% of baseline flow during ischemia and approximately 10% of baseline flow immediately after release in the vehicle group; these are normalized perfusion units (% of pre-ischemic baseline) and are reported approximations, not digitized observations. Table 2 reports control gait score `1.70 +/- 0.67` and NIM-811 gait score `2.70 +/- 0.82` (unitless Tarlov score). The functional measurement was obtained after recovery, so it is an independent early/clinical anchor rather than a direct 15 min model label.
2. McAllister et al., *American Journal of Physiology Regulatory, Integrative and Comparative Physiology* 295, R681–R689 (2008), DOI `10.1152/ajpregu.90303.2008`, PMID 18509099. In pig muscle, postconditioning or mPTP inhibition reduced infarction and was associated with lower mitochondrial Ca2+ and higher muscle ATP after ischemia/reperfusion. This supports the proposed Ca2+ -> mPTP -> ATP/function chain but supplies no coefficient for this normalized model.
3. Naparus et al., *European Journal of Pharmacology* 686, 90–96 (2012), DOI `10.1016/j.ejphar.2012.04.045`, PMID 22575519. Human rectus abdominis strips underwent 3 h hypoxia/2 h reoxygenation; mPTP inhibition increased ATP content and reduced injury. This supports retaining a human-relevant ATP/function state, not a transferred numerical rate.
4. Tran et al., *PLoS ONE* 7, e43410 (2012), DOI `10.1371/journal.pone.0043410`, PMID 22912870. In a mouse tourniquet model, mitochondrial superoxide, mPTP opening, and apoptosis increased after ischemia/reperfusion; antioxidant or mPTP interventions reduced the linked injury signals. This supports the ROS -> mPTP coupling.
5. Pottecher et al., *Journal of Vascular Surgery* 57, 1100–1108.e2 (2013), DOI `10.1016/j.jvs.2012.09.020`, PMID 23332985. The reported skeletal-muscle maximal oxidative capacity was `4.08 +/- 0.38` versus sham `5.98 +/- 0.56 umol O2/min/g` after 3 h ischemia/2 h reperfusion, and ROS was `3992 +/- 706` versus `1812 +/- 322` arbitrary units. These values motivate tracking oxidative capacity and ROS, but are not inserted as rates because tissue, species, and assay differ.

## Hypothesis

A flow history alone is insufficient. The minimum mechanistically informative state at reperfusion is the joint state of residual regional flow/oxygen delivery, ATP/energy deficit, calcium overload plus ROS-triggered mPTP opening, and the resulting edema/injury/functional reserve. ATP and calcium/ROS are the latent injury states; flow is the delivery state; function is the downstream observable. Inflammatory signaling is deliberately not needed for a first 15 min mechanistic prediction and is left for the next resolution step.

More specifically, with the same external flow history, two tissues with different ATP and Ca2+/ROS states should have different early function and mPTP trajectories. Blocking mPTP should preserve function and improve microvascular flow recovery without changing the external flow schedule.

## Frozen prediction and acceptance criterion

The primary predicted quantity is the normalized change in contractile function from the pre-ischemic baseline to 15 min after reperfusion, `DeltaF15 = F15 - Fpre`. The preregistered sign is negative for the full ischemic state and exactly zero for the flow-only ablation; secondary predictions are `M15 > 0.10`, `A15 < 0.90`, and a persistently reduced actual flow `q15 < q_ext`. These are predictions of the model structure, not clinical effect estimates.

The following criterion is frozen before the first model run:

- With no ischemia, the state must remain at the healthy equilibrium (maximum absolute state deviation below `1e-8` over 30 min) and function must remain `1.0` within numerical tolerance.
- With 90 min ischemia at the verified low-flow challenge followed by external reperfusion, at 15 min the full model must have `mPTP > 0.10`, `ATP < 0.90`, and `function < 0.90` relative to the pre-ischemic baseline. These are directional, normalized internal checks, not empirical accuracy claims.
- In the same challenge, an 85% mPTP-opening blockade must lower `mPTP` by at least `0.05` and improve function by at least `0.05` at 15 min. This is a falsifiable mechanistic placebo/intervention check.
- A flow-only ablation (metabolic states held at healthy values) must not create the mPTP/ATP/function injury pulse. Its function must remain within `0.02` of baseline at 15 min while the full model must be below `0.90`.
- A numerical solution must finish within 1 GB and produce finite, nonnegative fraction states (calcium may exceed 1 but must remain below 5).

Failure of any criterion is reported as a failed preregistration, not repaired by changing the criterion. The current run is a model-structure check; it does not validate the model against human or internal measurements.

## Frozen model structure

The implementation uses one regional compartment and the following state vector:

`[q, o, A, Gly, L, H, C, R, M, D, E, F]`

where `q` is normalized regional flow, `o` normalized tissue oxygen availability, `A` ATP fraction, `Gly` available glycolytic substrate fraction, `L` lactate proxy, `H` acidosis proxy, `C` normalized intracellular calcium, `R` ROS proxy, `M` mPTP-open fraction, `D` irreversible injury fraction, `E` edema proxy, and `F` contractile-function fraction. All are dimensionless fractions/proxies; time is minutes and all rates are `min^-1` or minutes as listed in `model.py`.

The external flow schedule is `q_ext=1` before ischemia, `q_ext=0.01` during the 90 min challenge, and `q_ext=1` after release. Actual flow relaxes toward `q_ext` multiplied by an injury-dependent no-reflow conductance. Oxygen relaxes toward delivery and is consumed by respiration. ATP is produced by oxygen-dependent respiration and glycolysis and consumed by a demand that falls during oxygen deprivation. Glycolytic substrate, lactate, and acidosis track anaerobic metabolism. Calcium rises when flow is absent and is cleared by flow. ROS is generated by reoxygenation and damaged mitochondria and cleared by flow-dependent antioxidant capacity. mPTP opening is a bounded opening/closing process gated by calcium excess, ROS, and ATP deficit. Injury and edema feed back on flow. Function relaxes toward an ATP-, calcium-, mPTP-, edema-, and injury-dependent target.

The implementation is allowed to make the following explicit assumptions: normalized Monod oxygen dependence, first-order relaxation, mass-action-like gates, and a bounded target function. No parameter is claimed to be a measured human rate. Source IDs in `model.py` distinguish `definition`, `literature_anchor`, and `assumption`.

## Null models and errors

- Null/placebo 1: external flow is restored but all metabolic states are held at their healthy values (`flow_only`).
- Null/placebo 2: mPTP opening is reduced by 85% without changing the ischemia duration or external flow (`mptp_blocked`).
- Wrong-model controls: a no-oxygen-delivery response and a no-mPTP response are not accepted as explanations if function remains healthy while the full state is injured.
- A run is an error if integration fails, states become nonfinite, the no-ischemia equilibrium drifts beyond the frozen tolerance, or a sensitivity multiplier produces an invalid state.

## Planned sensitivity and next resolution step

After the first run, vary `vmax_resp`, `k_mptp_open`, and `k_no_reflow` independently to `0.5x` and `1.5x`; report function, mPTP, and actual flow at 15 min. The next resolution step is calibration/validation against a time series containing regional flow, oxygen or oxygen-consumption proxy, ATP/metabolites, and an independent function anchor. Inflammation, endothelial barrier dynamics, and subject-specific perfusion geometry are deliberately deferred.
