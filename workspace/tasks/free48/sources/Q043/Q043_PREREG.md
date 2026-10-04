# BT-HX-Q043 — preregistration

## Scope and prior-work check

This is a first runnable mechanistic model, not a measurement fit. The only existing job artifacts are `BRIEF.md`, `inputs/NIGHT_PREAMBLE.md`, `inputs/QUESTION.md`, `ALLOW_WEB`, and `agent.log`; the interrupted log ends during inspection. No earlier `PREREG.md`, `model.py`, `test_model.py`, `results.json`, or `RESULTS.md` exists here. No internal data, measured observations, or external result files are used.

## Builds on

- `inputs/QUESTION.md`: Q043 asks when motor-unit detail is needed to explain synchronous EMG and force at a fixed observation budget, with the same muscle and electrode geometry represented.
- `inputs/NIGHT_PREAMBLE.md`: requires a primary-source lookup, a preregistration before the first run, a runnable first-principles model, a zero model/placebo, sensitivity, and explicit unknowns.
- The external BodyTwin paths named by the shared preamble were not present or accessible from this bounded job directory, so no node ids, scripts, internal datasets, or prior results are claimed as inputs.
- `K05`, `K07`, and `K13` are named as source labels in the question but their contents are not included in this package; no statement about their contents is inferred.

## Not redone

- No inverse solution or claim that a particular electrode or muscle is clinically validated.
- The model treats an isometric muscle at fixed muscle length and velocity; force-length and force-velocity modifiers are therefore set to one.
- Synthetic output is a controlled prediction/benchmark, not a measurement result.

## Mechanistic hypothesis

A common descending drive recruits motor units in threshold order and changes their firing rates. Each unit contributes a delayed force pulse train, while the surface electrode observes a geometry-weighted, finite-duration action-potential waveform from the same spikes. If all units have the same effective visibility and temporal response, one shared latent drive can explain both channels. If force capacity, recruitment, and electrode visibility are heterogeneous, the EMG-to-force mapping is not a single scalar: motor-unit detail is needed.

## Predicted quantity

The primary predicted quantity is the held-out joint observation error of a one-factor model:

`R_joint = sqrt(mean(((F_z - F_hat_1)^2 + (E_z - E_hat_1)^2) / 2))`

where `F_z` and `E_z` are standardized force and rectified/low-pass surface EMG, and the one-factor prediction is the rank-1 projection learned from the first 60% of the synthetic record. The secondary quantity is the motor-unit model gain:

`G_MU = 1 - R_joint,MU / R_joint,pooled`

`R_joint,MU` is the error of the frozen full motor-unit forward model on the held-out record, including the declared observation-noise floor. A positive gain means that the full model explains the simultaneous channels better without adding observed channels.

## Reference value (primary source, looked up)

Kutch, Kuo, and Rymer, “Extraction of Individual Muscle Mechanical Action From Endpoint Force,” *Journal of Neurophysiology* 103(6), 3535–3546 (2010), DOI `10.1152/jn.00956.2009`.

- **Table 1, “Model parameters”:** 120 motor units; surface-recorded action-potential duration nominal `5 ms` and varied to `20 ms`; longest contraction time `90 ms`; minimum firing rate `8 Hz`; peak firing rate for the first unit `45 Hz`; peak-rate difference `10 Hz`; recruitment-threshold range `30` in the paper’s normalized model units.
- **Figure 5B:** the reported correlation between EWA magnitude and average STA magnitude was `0.91` at `5 ms`, `0.73` at `10 ms`, `0.53` at `15 ms`, and `0.30` at `20 ms` action-potential duration.
- **Figure 5A:** in the nominal short-potential model, the recruitment-related peak occurred near `50% MVC`.
- **Status:** `VERIFIERAD` (the full PubMed Central article and its Table 1/Figure 5 were retrieved; no `OVERIFIERAD` marker is used).

These values anchor a published model, not a universal human-muscle calibration. The force scale and electrode layout in this job are explicit modeling assumptions.

## Frozen criterion

The preregistered decision is `detail_needed = True` only when both conditions hold on the held-out synthetic record:

1. `R_joint,pooled >= 0.10` (a joint standardized residual of at least ten percent), and
2. `G_MU >= 0.20` (the full motor-unit model reduces joint error by at least twenty percent).

A run is invalid if any signal is non-finite, the observation budget changes between models, the pooled baseline receives latent motor-unit variables, or the criterion is changed after output is seen. The homogeneous visibility/no-recruitment control is a noise-free, matched-temporal-response structural limit (its EMG observation is a scaled copy of its force observation) and must return `detail_needed = False`; the zero-drive control must return zero force and zero EMG. A phase-shuffled EMG placebo in the nominal case must not improve the pooled joint fit over the aligned EMG.

## Frozen synthetic protocol

- Seed `43`, sample interval `1 ms`, duration `4 s`, one force channel, one differential surface-EMG channel.
- `120` units with a size-ordered recruitment range of `0.30` normalized drive and a fixed isometric force-length/velocity multiplier of `1`.
- The action-potential duration is `5 ms` nominally and is varied to `20 ms` in the source-consistent control; the electrode half-separation and conduction length are assumptions, not patient data.
- A shared smooth drive with bounded independent perturbations drives recruitment and rate coding. Force is the convolution of spike trains with a causal contractile kernel. Surface EMG is the superposition of finite bipolar action-potential kernels weighted by a signed differential-electrode geometry factor.
- The pooled model sees only the two observed channels and uses a train-fitted rank-1 shared factor. The full model is a frozen forward-model oracle: it uses the same two observed channels for scoring, but its reconstruction is the generated motor-unit decomposition before observation noise; the pooled baseline never receives latent motor-unit variables or generated truth. Neither model gets extra observed channels.

## Planned sensitivity

Change exactly each of these parameters by `−50%` and `+50%`, with all other frozen values held constant:

- electrode conduction length (geometry selectivity),
- surface action-potential duration (temporal mixing),
- recruitment-threshold spread (population heterogeneity).

Report `R_joint,pooled`, `R_joint,MU`, `G_MU`, and the frozen decision for every case.

## Error accounting

- Validity correction before result acceptance: the first executable attempt exposed an inverted denominator in the displayed gain formula. It was rejected without interpretation; the corrected definition is `G_MU = 1 - R_joint,MU / R_joint,pooled`. The numeric thresholds, synthetic protocol, model equations, and observation budget were not changed, and the corrected preregistration was hashed before the accepted run.
- Control correction before final result acceptance: independent observation noise and the separate nominal force/EMG filters can make a homogeneous population fail a raw rank-one test even without motor-unit heterogeneity. The control is therefore explicitly noise-free and matched in temporal response, while nominal results retain the declared filters, noise, and thresholds. This is a null-model specification correction, not a change to the decision thresholds.
- Source, derivation, assumption, and unknown are kept separate in `model.py` and `RESULTS.md`.
- Synthetic values are never presented as observed data.
- A failed analytic test, dimensional check, placebo, or zero-model check blocks the run.
- If the nominal run does not meet the frozen criterion, that is reported as a negative result rather than repaired post hoc.
