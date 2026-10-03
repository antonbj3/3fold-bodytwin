# BT-HX-Q127 — preregistration

Status: frozen 2026-09-25 before `model.py` was written or run.
Scope: one representative adult type-II muscle fibre, with a longitudinal fibre identity and a satellite-cell precursor pool. This is a first mechanistic model, not a clinical or motor-memory model.

## Builds on

- `inputs/QUESTION.md`: the requested state variables, `dM/dt = J_fuse - k_loss M`, the warning that fusion yield is dimensionless per event, and the required outputs.
- `inputs/NIGHT_PREAMBLE.md`: the bounded, one-thread, public-literature-only execution and the separation of source, derivation, and hypothesis.
- The only pre-existing session artifact was `agent.log`; it records an interrupted directory check and contains no model or results. No earlier deliverable was present.
- Public reference: Cumming, K. T., Reitzner, S. M., Hanslien, M., Skilnand, K., Seynnes, O. R., Horwath, O., Psilander, N., Sundberg, C. J., and Raastad, T. (2024), “Muscle memory in humans: evidence for myonuclear permanence and long-term transcriptional regulation after strength training”, *The Journal of Physiology* 602(17):4171–4193, DOI `10.1113/JP285675`. The Ovid full-text page was fetched and the values below occur in its “Number of myonuclei” Results and Figs. 3 and 6: https://www.ovid.com/journals/jphy/fulltext/10.1113/jp285675~muscle-memory-in-humans-evidence-for-myonuclear-permanence. The page was found; this reference is not marked OVERIFIERAD.

## Hypotes

A load episode first increases a fibre’s precursor-to-fibre fusion flux when its cytoplasmic domain is stressed. The number of nuclei can subsequently remain above the no-history control during detraining even when fibre area returns toward baseline. A slowly decaying intracellular priming state can alter per-nucleus synthesis on retraining, but that state is not labelled motor memory and is not inferred from the nuclear count alone.

The model will make no claim that retained nuclei cause faster functional retraining. It will report nuclear retention, nuclear density, fusion-event flux, and a size-matched/history-ablated counterfactual separately.

## Previous session

The directory initially contained only `BRIEF.md`, `ALLOW_WEB`, `agent.log`, and `inputs/`. The interrupted log had no substantive run. This preregistration is therefore the first frozen model specification; no prior parameter or result was reused.

## Quantity and reference

Primary predicted quantity:

`R_excess = (M_D - M_C,D) / (M_1 - M_0)`, dimensionless in [0,1], where `M` is myonuclei per fibre, `M_0` is the pre-training value, `M_1` is the end of first training, `M_D` is the start of retraining after detraining, and `M_C,D` is the concurrent no-first-training control at the same detraining time point. This measures retention of the training-associated excess rather than the total baseline count.

Published means used for the reference:

- `M_0 = 2.4 ± 0.5` myonuclei/type-II fibre at baseline.
- `M_1 = 3.3 ± 0.7` myonuclei/type-II fibre after 10 weeks of first training.
- `M_D = 3.2 ± 0.6` myonuclei/type-II fibre in the previously trained arm at the start of retraining.
- `M_C,D = 2.4 ± 0.5` myonuclei/type-II fibre in the contralateral control arm at the same point.

The frozen derived reference is `R_ref = (3.2 - 2.4)/(3.3 - 2.4) = 0.8888889`. It is a derived ratio of published means, not a percentage directly reported by the authors. The source reports no significant type-II change during detraining (`-2 ± 5%`, `P = 0.276`) and a 33% higher trained-arm count than control at retraining start; the ratio is the prespecified operationalization of those facts.

## Frozen acceptance criterion

The deterministic protocol mean must satisfy `abs(R_pred - 0.8888889)/0.8888889 <= 0.30` (within ±30%). This criterion is applied to the primary retention ratio only. A result outside the interval is a failed prediction, not a reason to retune the parameters. A non-finite state, negative count, invalid dimension check, failed analytical test, or inability to complete the fixed protocol is also a failure. No claim about motor memory or strength is accepted by this criterion.

## Protocol and equations

The fixed schedule is 70 days first training, 112 days detraining, and 70 days retraining, with one-day Euler steps and a fixed 0–1 normalized load input. The stochastic run uses 2,000 fibres, seed `127`, and integer lineage events; the deterministic expectation is the primary prediction.

For fibre `i`, `M` is nuclei/fibre, `S` is activated precursors/fibre, `P` is a contractile-area-equivalent pool, `A` is fibre area, `G` is a dimensionless load-signal state, and `E` is a dimensionless intracellular history/priming state. The area-equivalent pool makes the protein balance dimensionally explicit.

- Load transport: `dG/dt = k_on u - k_off G`.
- Precursor balance: `dS/dt = a_base + a_load G - k_cycle S - J_fuse`.
- Domain: `D = A/M` (area/nucleus), and `p_fuse = sigmoid((D/D_star - 1)/w_D)`.
- Fusion: `J_fuse = k_fuse S p_fuse` (events/fibre/day). Each event consumes one precursor and contributes `y_fuse` nuclei; therefore `dM/dt = y_fuse J_fuse - k_loss M`, with no extra day factor multiplying `y_fuse`.
- Nuclear loss: each lineage has hazard `k_i = k_loss (1 + q_atrophy max(0, 1-A/A_0)) (1 - p_protect I_acquired)`. The aggregate approximation uses the corresponding mean hazard; the stochastic run samples binomial losses.
- Protein balance: `dP/dt = r_nuc M (1 + r_load G + r_history E G) - k_deg P`.
- Area mechanics: `dA/dt = k_inc P - k_turn A`.
- History: `dE/dt = k_learn G (1-E) - k_forget E`. It is an intracellular state variable, not a claim about central nervous-system memory.

New nuclei are assigned a unique lineage identifier, an acquisition day, and a position on the fibre. Candidate positions are sampled around a load/damage site and relaxed by local repulsion; positions are retained across phases so fibre identity is not reset. The stochastic transition draws `F ~ Poisson(J_fuse dt)` and `L ~ Binomial(M, 1-exp(-k_i dt))`, respecting the available precursor count. Integer events are used only for the stochastic lineage illustration; no fractional nucleus is reported as an observed count.

## Frozen parameter table

| Parameter | Value | Unit | Basis |
|---|---:|---|---|
| `M0` | 2.4 | nuclei/fibre | Cumming et al. 2024 type-II baseline mean |
| `A0` | 4029 | µm²/fibre | Cumming et al. 2024 type-II baseline fCSA mean |
| `fiber_length` | 10000 | µm | assumption; only axial spacing/density conversion |
| `P0` | 4029 | µm²-equivalent/fibre | assumption: initial pool equals initial area |
| `S0` | 0.10 | precursor/fibre | assumption; uncalibrated starting pool |
| `E0`, `G0` | 0 | dimensionless | assumption |
| `r_nuc` | 16.8 | µm²/(nucleus·day) | derived assumption `k_deg A0/M0` |
| `r_load` | 0.20 | dimensionless | assumption for load-enhanced synthesis |
| `r_history` | 0.35 | dimensionless | assumption; ablation is run |
| `k_deg` | 0.010 | day⁻¹ | assumption |
| `k_inc` | 0.010 | day⁻¹ | assumption |
| `k_turn` | 0.010 | day⁻¹ | assumption |
| `k_on` | 0.120 | day⁻¹ | assumption |
| `k_off` | 0.120 | day⁻¹ | assumption; steady load signal is 1 |
| `a_base` | 0.001 | precursor/(fibre·day) | assumption |
| `a_load` | 0.030 | precursor/(fibre·day) | assumption |
| `k_cycle` | 0.050 | day⁻¹ | assumption |
| `k_fuse` | 0.035 | event/(precursor·day) | assumption; yield is separate and dimensionless |
| `D_star` | 1200 | µm²/nucleus | assumption for domain-sensitive fusion |
| `w_D` | 0.20 | dimensionless | assumption for sigmoid width |
| `y_fuse` | 1.0 | nuclei/event | assumption for one satellite-cell fusion event |
| `k_loss` | 0.00035 | day⁻¹ | assumption; not fitted to the reference |
| `q_atrophy` | 1.0 | dimensionless | assumption |
| `p_protect` | 0.50 | dimensionless | assumption for acquired-lineage protection |
| `k_learn` | 0.030 | day⁻¹ | assumption |
| `k_forget` | 0.002 | day⁻¹ | assumption |
| `damage_sigma` | 0.12 | normalized fibre length | assumption |
| `repulsion` | 0.35 | dimensionless relaxation weight | assumption |
| `dt` | 1 | day | fixed numerical step |
| `stochastic_fibres` | 2000 | fibres | fixed Monte Carlo size |
| `seed` | 127 | dimensionless | fixed random seed |

The parameter table is frozen even where a value is deliberately an assumption. A failed criterion will be reported as a failed criterion; post hoc parameter changes are outside this preregistration.

## Counterfactual and falsification

The same retraining schedule is run from the no-history control. The primary history contrast is trained-history versus no-history at matched retraining load; a second run sets `r_history=0` while retaining the acquired nuclei. If the excess retention ratio fails ±30%, if the no-fusion/no-loss limit fails, or if the history contrast is indistinguishable from the size-only counterfactual, the mechanistic memory interpretation remains UNKNOWN.

## Not redone

This bounded run deliberately reuses only the question's state variables, the public reference values, and the equations listed above. It does not reuse internal BodyTwin data, external solver runtime, individual trajectories, measured fusion rates, or unpublished geometry. It does not reinterpret nuclear retention as motor memory. The post-run edits to this file only add this explicit scope section; no parameter, schedule, reference, or acceptance criterion changed.

## Not treated as resolved

No internal BodyTwin data, individual trajectories, measured fusion-event rates, fibre geometry, fibre-type mixture, or epigenetic measurements are used. The model does not infer causation from cross-sectional biopsy means, and it does not convert nuclear retention into functional motor memory. The next resolution step is longitudinal fibre-identified PCM1/DAPI imaging with a registered biopsy coordinate, damage-site labels, and a matched no-training control.

## Planned outputs

`model.py` performs the deterministic and fixed stochastic runs, finite ±50% sensitivity analysis for the three prespecified parameters `k_fuse`, `k_loss`, and `r_history`, and writes `results.json`. `test_model.py` checks the zero-flux invariant, the analytic exponential-loss limit, and the dimensional audit. `RESULTS.md` reports the frozen comparison and limitations. No literature value is treated as an internal measurement.
