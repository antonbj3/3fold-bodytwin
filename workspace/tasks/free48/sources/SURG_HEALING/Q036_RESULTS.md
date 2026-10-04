BT-HX-Q036

# Results

## Core answer

The first executable model says that the early reperfusion response cannot be determined by flow alone. The smallest mechanically informative state vector is:

1. actual regional flow `q` and available oxygen `o` (delivery state);
2. ATP/energy deficit `A` (reserve state);
3. calcium overload `C` together with ROS `R` (mPTP trigger state);
4. mPTP `M`, damage `D` and oedema `E` feeding back to flow and function `F` (reserve and functional state).

It is this connected state that the model tests. The primary predicted quantity is `DeltaF15 = F15 - Fpre`: the full model is predicted to give a negative value, while flow-only is predicted to give zero. Inflammation is not included in the first resolution, because the question concerns early response and the mechanical chain can be described without an inflammation compartment. No measured data has been used to calibrate any parameter.

## Reproduction

- `PREREG.md` and `PREREG.sha256` were created before the first run; the final hash is `1042448a06adcbf76e41cea35a367b537d884c7e25e134a9a81aed0e69eb0885`.
- `python3 model.py --output results.json` runs the preregistration hash check and writes `results.json`.
- `python3 -m unittest -v test_model.py` runs six tests, including an analytical mPTP limiting case.
- `python3 -m py_compile model.py test_model.py` is the available syntax check; no separate lint or typecheck configuration exists in the directory.

## Model from first principles

Time is in minutes. The state vector is

`y = [q, o, A, Gly, L, H, C, R, M, D, E, F]`.

All states are normalised quantities or proxies. `q`, `o`, `A`, `Gly`, `M`, `D`, `E` and `F` are kept within `[0,1]`; the calcium proxy `C` is kept within `[0,5]`; `L`, `H` and `R` are unbounded concentration proxies. The bounds are numerical/model definitions, not measured data.

### Flow and oxygen

The external flow separator is `q_ext=1` before ischaemia, `q_ext=0.01` during ischaemia and `q_ext=1` after release. The mPTP/oedema-related conductance loss is

`b = 1 - (1 - k_nr M)(1 - k_ed E)`,

`qdot = (q_ext (1-b) - q)/tau_q`.

This makes actual flow different from having restored external perfusion. Oxygen delivery and consumption are

`x = o/(K_M + o)`,

`R = V_max x A (1-M)/(1+2H)`,

`odot = k_on(q q_art-o) - k_use R`.

### ATP, glycolysis and metabolic proxies

The hypoxic/glycolytic port is

`h = max(0,min(1,(K_M-o)/K_M))`,

`G = g_max Gly h A (1-M)/(1+2H)`.

`R0 = V_max o0/(K_M+o0)` is the healthy respiratory equilibrium rate. Energy demand and ATP balance are

`U = R0 (0.35+0.65F)(0.85+0.15x)`,

`Adot = R + G - U + k_rec q(1-A)`,

`Glydot = -G + k_gly q(1-Gly)`,

`Ldot = y_L G - k_L q L`,

`Hdot = y_H G + k_H(1-x)(1-A) - k_acid q H`.

ATP thus cannot simply be restored by opening flow; substrate and previous energy loss are part of the state.

### Calcium, ROS and mPTP

Calcium and ROS dynamics are

`Cdot = k_Ca,in(1-q)(1+0.8(1-A)) + k_Ca,rel(1-q)(1-A) - k_Ca,clear q(C-1) - k_Ca,mito(C-1)_+`,

`Rdot = k_ros,reox q x(1-A)(1-M) + k_ros,mito M + k_ros,anox(1-q)(1-A) - k_ros,clear(0.2+0.8q)R`.

Calcium and ATP excess are defined as `Ca_ex=(C-1)_+` and `E_atp=(1-A)_+`. The mPTP port is

`g_Ca = Ca_ex/(K_Ca+Ca_ex)`,

`g_ros = R/(K_ros+R)`,

`g_E = E_atp/(K_E+E_atp)`,

`trigger = g_Ca [k_priming+(1-k_priming)g_ros][0.30+0.70g_E]`,

`Mdot = (1-b_mPTP) k_open trigger(1-M) - k_close M`.

Here `b_mPTP=0` in the untreated model and `b_mPTP=0.85` in the placebo/intervention case. ROS and mPTP feed back to irreversible damage and oedema, and these feed back to flow.

### Function as a downstream observable

`F` relaxes towards a functional capacity multiplicatively limited by ATP, calcium, mPTP, oedema and damage:

`F_target = (0.15+0.85A) exp(-0.65 Ca_ex)(1-0.65M)(1-0.45E)(1-D)`,

`Fdot = (F_target-F)/tau_F - k_stroke(Ca_ex+M)[0.5+0.5(1-q)]`.

This makes `F` a functional proxy, not an observed Tarlov or LVDP value.

## Parameter table and unit check

The complete table with value, unit, source and assumption is in `model.py` (`PARAMETER_META`) and serialised in `results.json` (`parameter_table`). The most important classes are:

| Parameter | Value | Unit | Source/assumption |
|---|---:|---|---|
| `q_art` | 1.0 | normalised flow [1] | definition |
| `q_ischemia` | 0.01 | normalised flow [1] | E1, approximately 1% of baseline during ischaemia |
| `tau_flow` | 0.35 | min | assumption |
| `k_oxygen_on` | 1.20 | min^-1 | assumption |
| `k_oxygen_use` | 0.70 | min^-1 | assumption |
| `vmax_resp` | 0.12 | min^-1 | assumption |
| `km_oxygen` | 0.20 | normalised oxygen [1] | assumption |
| `gmax` | 0.055 | min^-1 | assumption |
| `k_ca_influx` | 0.035 | min^-1 | assumption |
| `k_ros_reox` | 0.180 | min^-1 | assumption |
| `k_mptp_open` | 0.450 | min^-1 | assumption |
| `k_mptp_close` | 0.160 | min^-1 | assumption |
| `k_no_reflow` | 1.35 | dimensionless | assumption, direction supported by E1 |
| `k_injury` | 0.008 | min^-1 | assumption |
| `k_edema` | 0.004 | min^-1 | assumption |
| `tau_function` | 0.650 | min | assumption |

The unit test in `check_units()` verifies that every parameter has a registered unit and that every ODE term has the time unit `min^-1` or an expression explicitly containing `/min`. `results.json` contains the complete check and `unit_consistent=true`.

## Verified literature anchors

No anchor is marked `OVERIFIERAD`.

- **E1:** El Baradie et al., *Scientific Reports* 11, 6152 (2021), DOI `10.1038/s41598-021-85753-x`, PMID 33731782, PMCID PMC7969970. 90 minutes of hindlimb ischaemia; the results text/Fig. 2–3 states approximately 1% of baseline flow during ischaemia and approximately 10% immediately after release in the control group (flow: % of pre-ischaemic baseline). Table 2 states gait score `1.70 +/- 0.67` in control and `2.70 +/- 0.82` with NIM-811 (dimensionless Tarlov score). Function was measured after recovery and is therefore not a 15-minute value for this model.
- **E2:** McAllister et al., *Am J Physiol Regul Integr Comp Physiol* 295, R681–R689 (2008), DOI `10.1152/ajpregu.90303.2008`, PMID 18509099. The abstract reports that postconditioning/mPTP inhibition reduced muscle infarction and gave lower mitochondrial Ca2+ and higher ATP; mechanistic motivation, not a transferred coefficient.
- **E3:** Naparus et al., *Eur J Pharmacol* 686, 90–96 (2012), DOI `10.1016/j.ejphar.2012.04.045`, PMID 22575519. Human rectus-abdominis strips: 3 h hypoxia and 2 h reoxygenation; mPTP inhibition preserved ATP and reduced damage.
- **E4:** Tran et al., *PLoS ONE* 7, e43410 (2012), DOI `10.1371/journal.pone.0043410`, PMID 22912870. Mouse tourniquet model: mitochondrial superoxide, mPTP opening and apoptosis increased after I/R; the ROS coupling is source-supported.
- **E5:** Pottecher et al., *J Vasc Surg* 57, 1100–1108.e2 (2013), DOI `10.1016/j.jvs.2012.09.020`, PMID 23332985. After 3 h ischaemia/2 h reperfusion, Vmax `4.08 +/- 0.38` against sham `5.98 +/- 0.56 umol O2/min/g` and ROS `3992 +/- 706` against `1812 +/- 322` AU were reported. These values motivate including metabolism and ROS in the model but have not been inserted as normalised rates.

## Run and results

The default case is 5 minutes pre-ischaemia, 90 minutes `q_ext=0.01`, and 30 minutes of external reperfusion schedule. All values below come from `results.json` and are model output, not measured data.

| State | q | o | A | C | M | F |
|---|---:|---:|---:|---:|---:|---:|
| pre-ischaemia | 1.0000 | 0.9423 | 1.0000 | 1.0000 | 0.0000 | 1.0000 |
| end of ischaemia | 0.00063 | 0.00063 | 0.0000 | 2.7895 | 0.6717 | 0.0000 |
| start of reperfusion | 0.00063 | 0.00063 | 0.0000 | 2.7895 | 0.6717 | 0.0000 |
| 15 min after release | 0.06594 | 0.06561 | 0.0000 | 2.4218 | 0.6654 | 0.0000 |
| 30 min after release | 0.07077 | 0.07055 | 0.0000 | 2.2425 | 0.6600 | 0.0000 |

At 15 minutes, the full model's flow is only 6.59% of baseline, which is a mechanical result and not a claimed measurement. mPTP remains high despite the external reperfusion schedule, and function reaches zero. This is a severe 90-minute scenario, not a general prognosis.

### Counter-cases

- `flow_only` keeps metabolic states healthy but lets the same external flow schedule pass. At 15 minutes, `q=1.0000`, `A=1.0000`, `M=0.0000`, `D=0.0000`, `E=0.0000` and `F=1.0000`. Flow alone thus cannot create the model's damage.
- `mptp_blocked` uses 85% mPTP opening blockade. At 15 minutes, `q=0.6768`, `o=0.6646`, `C=1.1774`, `M=0.1590`, `D=0.2970`, `E=0.3760` and `F=0.0646`. The difference from the full model is a mechanistic control outcome, not a clinical NIM-811 effect.

### Frozen criterion

`results.json` states `all_passed=true` for all frozen checks:

- healthy equilibrium: maximum absolute deviation `0.0`, function `1.0` before and after 30 minutes;
- full model at 15 minutes: `M>0.10`, `A<0.90`, `F<0.90`;
- 85% mPTP block: `M` decreases by at least `0.05` and `F` improves by at least `0.05`;
- flow-only: `F` is within `0.02` of `1.0` while the full model is below `0.90`;
- all numerical solutions are finite and the states are projected to valid model bounds.

## Sensitivity, ±50%

The following was run independently for `vmax_resp`, `k_mptp_open` and `k_no_reflow`, with other parameters and the entire flow history fixed. Values are at 15 minutes after release.

| Parameter | 0.5x: `M` / `q` | 1.0x: `M` / `q` | 1.5x: `M` / `q` | 1.5x minus 0.5x, `M` |
|---|---:|---:|---:|---:|
| `vmax_resp` | 0.6637 / 0.06988 | 0.6654 / 0.06594 | 0.6657 / 0.06575 | 0.00197 |
| `k_mptp_open` | 0.4675 / 0.2745 | 0.6654 / 0.06594 | 0.7558 / ~0 | 0.2882 |
| `k_no_reflow` | 0.5921 / 0.4789 | 0.6654 / 0.06594 | 0.6727 / ~0 | 0.0806 |

In this long, severe scenario, ATP and function reach the lower model bound for all three parameters. Therefore, `M` and `q` are the informative sensitivity outcomes here; `vmax_resp` has little effect at the particular 15-minute point because flow/mPTP damage dominates. It is a limitation of the design, not a claim that respiration is unimportant.

## What is not validated

There is no internal training or test data, no subject-specific perfusion, and no independent measurement of all latent states in the same experiment. Empirical validation and K01/K03/K09/K12 semantics are therefore explicitly `UNKNOWN`. The model must not be used for patient prognosis or effect size. The verified E1 function is also a Tarlov score after recovery, not the same endpoint as `F(15)`. No inflammation effects, endothelial barrier parameters, long-term glycogen dynamics or age/sex-dependent coefficients are calibrated.

## Next resolution step

The next step is a time series with regional flow and oxygen turnover, ATP or metabolite proxy, an mPTP/ROS or damage marker and an independent functional anchor at the same time points. The parameters should first be identified against that time series and then tested on a held-out individual/cohort, with age, sex and tissue dependencies as explicit hierarchical variables. Inflammatory signal, endothelial conductance and glycogen reserve should then be connected as separate modules.
