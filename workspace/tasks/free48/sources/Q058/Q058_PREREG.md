# BT-HX-Q058 — fryst preregistration

## Status and limit

This is written before the first run of `model.py`. The model is a small, deterministic and normalized isometric muscle model. It does not represent a particular person, a particular joint or a measurement data set. `F_max = 1` is used as the power unit; all relative errors are therefore fractions of maximum power. No training or post-adjustment is made after the criteria here are frozen.

## Hypotes

A switch between levels of detail is biologically memory-preserving if the coarser model (a) accommodates all slow states that may still affect the future and (b) initializes the fast states with the correct value, or causes them to approach their quasi-stationary value. If only the slow fatigue/capacity pools is transmitted, a residual central adaptation or metabolic fatigue history will produce an observation error after a later perturbation.

## Predicted quantity

Primary prediction is the post-perturbation error

`E_switch = max_{t_p <= t <= t_end} |F_switch(t) - F_continuous(t)| / F_max`

where `F_continuous` is a high-detail reference that never switches model, and `F_switch` switches from high to low at `t_s` before a short perturbation at `t_p`. Secondary prediction is that the error maximum decreases when the fast memory poke has caused the `n_f >= 5` time constants to relax, and that a mode that transfers even the snabba/slow-state-minnena ones is closer to the reference than a naive aggregated mode.

## Reference values from published primary source

Look up 2026-09-25 in the open full text for Potvin & Fuglevand (2017), *A motor unit-based model of muscle fatigue*, PLoS Computational Biology 13(6):e1005581, DOI **10.1371/journal.pcbi.1005581**:

- Fig. 1B and Methods Eq. (1): representative twitch-contraction time is in the range **30–90 ms**; the model therefore uses `c0 = 50 ms` as a middle value, not as new measurement data.
- Methods Eq. (11) and Fig. 1A–B: `ΔCT/CT = 0.379 × ΔF/F`; In case of `20 %` loss of power, this becomes **7,6 %** longer contraction time.
- Methods Eq. (12): the published time constant of the central firing-rate adaptation is **22 s** (`τ_h = 22 s`).
- Fig. 2: published simulated persistence at 20 % MVC is **511,5 s**; Fig. 3: corresponding to 50 % MVC is **95,5 s**. These two values ​​are the literature simulated reference results, not the data used for calibration in this model.
- Fig. 7: the published model loses about **1,4 % MVC/s** during the first 20 s and reaches about 50 % of initial power after **70 s** in its 100 %-MVC-fall.

The reference values ​​are thus verified; `OVERIFIERAD` is not used. The mechanistic parameters in the table in `model.py` explicitly distinguish between published values ​​and assumptions. This model must not claim to reproduce 511,5 s or 95,5 s.

## Fryst kriterium

A switch case counts as memory conservation in the first iteration if all of the following apply:

1. `max(0, F_switch(t_s^+) - F_continuous(t_s)) <= 0,02 F_max` and the corresponding negative leap is `>= -0,02 F_max`.
2. `E_switch <= 0,05 F_max` during the entire 60 s post-perturbation window.
3. all normalized states lie in `[0, 1]` (capacity must not become negative) and both models use the same `F_max`.
4. the naive shift (`m_s` alone) is compared against a reference with both fatigue pools; the difference is reported even if it passes the criterion.

A case that fails 1 or 2 is `FAIL`, not a case that gets a new threshold value. A passed naïve case is not proof that the memory is complete; it is only a numerical condition for the defined perturbation.

## Driving protocol and placebo/counter test

- initial rest `0–5 s`, excitation pulse `u = 0,65` during `5–35 s`, rest `35–55 s`, reactivity `55–80 s`, and perturbation `u: 0,65 -> 0,90` during `80–80,1 s`; then the same drive until `140 s`.
- highly detailed model: activation, excitation–contraction filter, central adaptation, two fatigue pools and fatigue-dependent contraction time.
- low-detailed model: activation and a slow fatigue pool; central adaptation, fast pool and filter memory are missing.
- standard switching occurs at `t_s = 79,5 s`, immediately before the perturbation at `t_p = 80,0 s`.
- placebo: same initial condition but switching at `t=0`; it should be almost identical to the height model because the history-dependent state is then zero.
- No data adaptation, no randomization and no use of internal BodyTwin-/the reference model-data.

## Sensitivity

Frozen ±50 % for `τ_f`, `k_f` and `τ_s`. Each parameter is changed alone, other parameters are kept frozen. Report `E_switch`, capacity jump and relative sensitivity `(Y(+50%)-Y(-50%))/(2Y0)`. This is a sensitivity analysis, not a measurement range claim.

## Builds on

- `inputs/QUESTION.md`, Q058 and its inputs K03, K09, K10; the substrate only contains labels, not their file contents.
- Potvin & Fuglevand (2017), DOI 10.1371/journal.pcbi.1005581, especially Fig. 1–3, Fig. 7 and Methods Eq. (1), (4), (5), (6), (10)–(14).
- Liu, Brown & Yue (2002), DOI 10.1016/S0006-3495(02)75580-X, in support that fatigue and recovery must be treated as separate dynamic processes; its 97 % statement is not used as a calibration.

## Not redone

No internal data, external solver runtime, patient parameters, cloud runs, LOSO-svep or assumptions that the 0,05 criterion is a clinical validation limit. No claimed measurement data is generated; `UNKNOWN` means that the quantification requires measurements that are not in the package.
