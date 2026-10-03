# BT-HX-Q083 — preregistration

## Hypotes

A fusion–fission event leaves a memory that is greater than total quantities in the merged organelle. To predict a later perturbation, the first model needs only per-daughter volume, per-daughter concentrations of soluble cargo and protein/regulatory state, membrane potential, and an explicit partitioning rule. Global totals and the number of organelles are a placebo/null model.

## Predicted quantity and frozen ground truth

The ground truth is the difference in expected response between two daughters after the same fusion–fission sequence and the same perturbation. The run must report:

1. time from fusion to fission (`tau_post`);
2. absolute potential gap between the daughters (`gap_psi`);
3. the ratio between the hyperpolarised and depolarised daughter's next fusion risk;
4. how much of this separation disappears if only total quantities are retained.

The canonical mitochondrial event is frozen at `tau_post = 77 s`, `gap_psi = 12 mV` and `riskkvot = 6`. Model parameters and initial states are stated below before the first run.

## Reference values from primary sources

The values are looked up in published literature and have not been changed after the run.

| Source | Primary observation | DOI and location | Use |
|---|---|---|---|
| Twig et al. 2008, *EMBO Journal* | INS1 cells: post-fusion/connected state 77 ± 71 s; COS7: 87 ± 78 s | `10.1038/sj.emboj.7601963`, Fig. 2B | `tau_post` |
| Twig et al. 2008 | In 85 % of 48 fissions, one daughter depolarised by +5,9 ± 3,7 mV while the other hyperpolarised by −6,5 ± 3,8 mV | `10.1038/sj.emboj.7601963`, Fig. 3B | `gap_psi = 5,9 + 6,5 = 12 mV` |
| Twig et al. 2008 | The hyperpolarised daughter fused within 3 min six times more often; mean interval 105 s | `10.1038/sj.emboj.7601963`, Fig. 4 | `riskkvot` |
| Eisner et al. 2017, *PNAS* | 1,4 ± 0,1 fusions/min in a freshly isolated adult ventricular cardiomyocyte; fast mixing <12 s and one representative slow mixing 70 s | `10.1073/pnas.1617288114`, Fig. 2C and Fig. 3A | time scales for soluble mixing |
| Luo et al. 2015, *Journal of Neurophysiology* | A single open Ca²⁺ channel triggers vesicle fusion with approximately 5–6 % probability | `10.1152/jn.00879.2014`, Fig. 5C | vesicle gate |
| Kim et al. 2012, *EMBO Journal* | Syt1 vesicles: `k_fusion = 3,7 ± 0,9 × 10^-3 s^-1`; without Syt1 `1,16 ± 0,25 × 10^-4 s^-1` | `10.1038/emboj.2012.57`, Fig. 3E and Fig. 3F | vesicle fusion rate |

The sources above are verified through PubMed/PMC text. `OVERIFIERAD` is not used.

## Frozen equations and states

Each compartment `i` has state

`S_i = (V_i, n_i, psi_i)`,

where `V_i` is volume in `um^3`, `n_i[s]` is quantity in `nmol` for species `s`, and `psi_i` is membrane potential in `mV`. Species used are `substrate`, `enzyme`, `damage`, `opa1` and `cargo`; the quantities are normalised quantities, not fabricated measured data.

During an open fusion pore, mass-conserving flux is used

`J_s = K_s (c_i,s - c_j,s)`,

`dn_i,s/dt = -J_s`, `dn_j,s/dt = +J_s`, `c_i,s = n_i,s / V_i`.

With `tau_s` defined as the time for the concentration difference to become `e^-1`,

`K_s = 1 / (tau_s (1/V_i + 1/V_j))`.

The same volume-weighted relaxation law with `tau_psi` applies to the potential. After `tau_post`, the merged volume is divided into `q` and `1-q`. For soluble species, deterministically `n_daughter1 = q n_total`. For protected regulatory/protein species (`enzyme`, `damage`, `opa1`, `activation`), a lobe-first partition is used: `a V_left + b V_right = q (V_left + V_right)`, where `a` is filled from the left lobe as far as possible and `b` afterwards. Each species has a small random partition deviation in the dissolved part with a known effective copy-number parameter. The potential is divided as

`psi_1 = psi_mean - 2(1-q) delta_psi`,

`psi_2 = psi_mean + 2q delta_psi`,

which preserves the volume-weighted mean potential for all `q`; at `q = 0,5`, the gap becomes `2 delta_psi`. It is a model rule, not a claimed measured value for every cell.

For mitochondria's future energy response, a simple enzyme/potential model is used

`J_ATP,i = V_i k_ATP c_enzyme,i c_substrate,i /(K_M + c_substrate,i) exp((psi_i - psi_ref)/psi_scale) exp(-c_damage,i) - V_i k_leak`.

For vesicles, the corresponding event model is used

`p_release = 1 - exp(-k_fusion dt activation calcium_factor)`

and expected cargo `cargo * p_release`. `activation` is a slowly mixing regulatory state. Synaptic vesicles that do not divide after exocytosis are handled as the limiting case where all remaining cargo is transferred to the target membranes; the reused topology is the same.

## Frozen parameter and initial values

| Parameter | Value | Unit | Source or assumption |
|---|---:|---|---|
| `V0` | 1,0 | `um^3` | normalised initial case |
| `tau_substrate` | 12 | s | Eisner Fig. 3A, fast mixing |
| `tau_enzyme` | 600 | s | assumption: protein/regulatory memory slower than soluble cargo |
| `tau_damage` | 600 | s | same assumption as `tau_enzyme` |
| `tau_opa1` | 600 | s | assumption; OPA1 is a relevant regulatory protein according to Twig |
| `tau_psi` | 5 | s | assumption for electrical coupling during fusion |
| `tau_post` | 77 | s | Twig Fig. 2B |
| `delta_psi` | 6 | mV | mean magnitude from Twig Fig. 3B |
| `q` | 0,5 | fraction | symmetric base case; sensitivity ±50 % |
| `k_ATP` | 0,08 | `nmol um^-3 s^-1` | assumption for normalised response |
| `K_M` | 0,8 | `nmol um^-3` | assumption |
| `psi_ref` | −180 | mV | assumption |
| `psi_scale` | 25 | mV | assumption |
| `k_leak` | 0,01 | `nmol um^-3 s^-1` | assumption |
| `K_opa1` | 0,60 | `nmol um^-3` | assumption for OPA1-dependent fusion risk |
| `n_opa1` | 4 | dimensionless | assumption: cooperative selectivity |
| `k_fusion_base` | 0,010 | `s^-1` | assumption for basal next-event hazard |
| `risk_window` | 180 | s | 3-minute window according to Twig Fig. 4 |
| `p_single_channel` | 0,06 | fraction | Luo Fig. 5C |
| `k_fusion_vesicle` | 0,0037 | `s^-1` | Kim Fig. 3E |
| `vesicle_dt` | 1 | s | defined prediction window |
| `substrate_perturbation` | 0,5 | fraction of initial substrate | frozen control case |
| `calcium_perturbation` | 0,5 | fraction | frozen control case |
| `vesicle_tau_post` | 30 | s | assumption for a division and reformation window |
| `vesicle_tau_cargo` | 30 | s | assumption: soluble vesicle cargo mixes quickly |
| `vesicle_tau_activation` | 600 | s | assumption: priming/SNARE state may persist |
| `vesicle_V0` | 0,010 | `um^3` | normalised vesicle volume |
| `partition_copy_number` | 10000 | copies | assumption for random partition |
| `lobe_species` | `enzyme, damage, opa1, activation` | category | assumption: protected/regulatory species retain lobe memory; soluble species are pooled |

Mitochondrial initial case before fusion:

| Dotter | `V` (`um^3`) | `substrate` (`nmol`) | `enzyme` (`nmol`) | `damage` (`nmol`) | `opa1` (`nmol`) | `psi` (`mV`) |
|---|---:|---:|---:|---:|---:|---:|
| A | 1,0 | 1,20 | 0,80 | 0,10 | 0,90 | −175 |
| B | 1,0 | 0,50 | 0,40 | 0,60 | 0,30 | −150 |

The vesicle case uses two states with the same volume and different `cargo` and `activation`; values are normalised and stated in the results JSON.

## Acceptanskriterier

The model is accepted only if all the following conditions hold:

1. mass conservation for every species and the potential's volume-weighted mean: relative error ≤ `1e-10`;
2. canonical mitochondrial `tau_post` lies within `50–100 s`;
3. `gap_psi` lies within `9–15 mV` and hyperpolarised/daughter risk ratio within `4–8`;
4. equal initial states give daughter-response difference ≤ `1e-10` (placebo);
5. the total-quantity model reduces the separation needed to explain a response gap by at least `50 %` compared with full state;
6. the analytical limit (complete mixing, equal volumes, no randomness) gives both daughters pooled concentration within `1e-10`;
7. no negative concentrations, no NaN and all reported unit strings match.

## Counter-tests and errors

Total-quantity null: retain only `sum(n)` and the number of daughters. Potential null: retain full composition but set all daughter potentials to the fusion mean. Placebo: two identical initial daughters. Without fusion: no topology event. An implementation that does not conserve mass, uses a calibrated number without a source, changes frozen criteria after running, or reports fabricated measured data counts as an error.

## Sensitivity

For `tau_enzyme` (`tau_damage` and `tau_opa1` change together), `delta_psi` and `q`, 0,5×, 1,0× and 1,5× are run. The daughter-response gap, full-state-minus-null error, potential gap and mass error are reported. No new calibration is performed after sensitivity has been run.

## Builds on

- `inputs/QUESTION.md`, Q083; stated entry points K03, K04 and K09.
- `inputs/NIGHT_PREAMBLE.md`, especially the requirement to distinguish source, derivation and hypothesis.
- `BRIEF.md`, bounded-task requirements and delivery order.
- No local BodyTwin anchor files or internal data were in the job packet, and read-only search outside the job folder was unavailable in the sandbox. No reuse of unknown internal node implementations is therefore claimed.
- Public primary sources: Twig 2008 (`10.1038/sj.emboj.7601963`), Eisner 2017 (`10.1073/pnas.1617288114`), Luo 2015 (`10.1152/jn.00879.2014`) and Kim 2012 (`10.1038/emboj.2012.57`).

## Not redone

- No external musculoskeletal solver run, no body-motion model and no estimate of joint force.
- No training fit, cloud run, LOSO or batch sweep.
- No claimed quantitative mtDNA copy, protein structure or individual measured data.
- No subtle cell-specific extrapolation from the published values; they are used as calibration/accountability references.
