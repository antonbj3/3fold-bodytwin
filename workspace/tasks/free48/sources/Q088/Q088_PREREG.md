# BT-HX-Q088 — preregistration

## Scope and build status

This is a first runnable, deliberately reduced mechanistic model of a skeletal-muscle fiber after an energy-demanding bout. It couples a radial oxygen field, ATP/phosphocreatine (PCr) balance, Na/K-ATPase-dependent ion restoration, and a preregistered functional recovery index. The modeled entity is an equivalent cylindrical fiber territory, not a measured human subject. No internal data are used.

Builds on:
- `inputs/QUESTION.md`, Q088: ATP recovery, ion balance, geometry, matched time courses, state history, geometry, and independent response data are the required representation.
- The interrupted-session directory contained no implementation, tests, results, or parameter files; the earlier `agent.log` only recorded the initial directory read. Nothing is silently reused.

Not redone:
- No individual-level fitting, no invented observations, no external musculoskeletal solver execution, no internal datasets, and no claim that the model has been validated in the target subject.
- Geometry is swept as a mechanistic hypothesis. Literature muscle cross-sectional areas are not silently converted into a measured single-cell radius.

## Hypothesis and primary prediction

The joint corner of low oxygen boundary tension, high contractile energy demand, and large diffusion radius should create the longest recovery because oxygen limitation reduces oxidative ATP resynthesis, PCr restoration, and ATP-dependent ion pumping at the same time. The primary frozen quantity is:

`T90_F = time after exercise until F(t) >= 0.90`

where `F = 0.50*PCr/PCr_rest + 0.30*ATP/ATP_rest + 0.20*(1-I)` and `I` is the normalized ion-gradient disturbance (`0` at rest, larger is worse). The secondary quantity is the initial recovery PCr resynthesis rate `V_PCr0`.

Before the reported run, the joint stress corner is defined as `pO2_boundary=2 kPa`, `work_demand=0.90 mM/s`, `R=40 micrometre`; the control is `pO2_boundary=13.3 kPa`, `work_demand=0.60 mM/s`, `R=20 micrometre`. The prediction is a joint-corner `T90_F` at least **90 s longer** than control and a `V_PCr0` no greater than **50%** of control. This is a prediction, not a measured effect size. If the functional threshold is not reached by 600 s, the observation is right-censored and the lower bound is recorded as 600 s.

## External reference values (looked up before implementation)

The values below are anchors for scale and for the qualitative sign checks. They are not silently treated as parameters for this fiber.

1. Layec et al. (2013), *Journal of Applied Physiology* 115:803–811, DOI `10.1152/japplphysiol.00257.2013`:
   - **Table 3**, free flow vs reactive hyperemia: PCr recovery tau `33 +/- 21 s` vs `27 +/- 10 s`; inferred peak ATP synthesis `28.7 +/- 13.3` vs `41.2 +/- 13.6 mM/min`.
   - **Table 4**, tissue reoxygenation mean response time: `70 +/- 15 s` vs `24 +/- 15 s`.
   - **Methods**, rest calibration: ATP concentration assumed to be `8.2 mM`; the paper's equation (3) uses `Km=30 micromolar` for the ADP-control calculation.
2. Heskamp et al. (2021), *The Journal of Physiology* 599:1533–1550, DOI `10.1113/JP280771`:
   - **Table 2**, continuous exercise distal vs proximal tibialis anterior: `k_PCr=0.44 +/- 0.26` vs `1.50 +/- 0.57 min^-1`; `V_PCr=5.2 +/- 3.1` vs `23.3 +/- 8.9 mM/min`; `k_O2Hb=5.4 +/- 3.8` vs `7.8 +/- 4.4 min^-1`.
3. Kushmerick et al. (1992), *PNAS* 89:7521–7525, DOI `10.1073/pnas.89.16.7521`:
   - **Abstract/full-text first page**, fast-twitch fiber metabolite values: ATP `8 mM`, total creatine `39 mM`, PCr `32 mM`, Pi `0.8 mM`, ADP `8 micromolar`. These values motivate the initial PCr/ATP pool scale, not a fiber-specific calibration.
4. Piiper and Scheid (1986), *Respiration Physiology* 64:241–251, DOI `10.1016/0034-5687(86)90118-0`:
   - **Abstract**, Krogh-cylinder and solid-cylinder comparison: skeletal muscle is better represented by the Krogh geometry, with capillary number/fiber number about `2` and capillary-to-fiber radius ratio about `0.1`. This motivates the radial diffusion topology and the explicit geometry hypothesis. The implementation uses an equivalent surface-supplied cylinder as a control volume; it does not claim a literal central-capillary geometry.
5. Segal and Faulkner (1985), *American Journal of Physiology* 248:C265–C270, DOI `10.1152/ajpcell.1985.248.3.C265`:
   - **Abstract**, rat muscle context: calculated critical oxygen-diffusion radius `1.19 mm` at `20 degrees C` and `0.51 mm` at `40 degrees C`. This is a whole-muscle context check, not a substituted single-fiber measurement.

All five sources were located in public literature records/full text. No value is marked `UNVERIFIED`.

## Frozen protocol

- Exercise phase: `250 s`; recovery observation: `600 s`; fixed numerical step: `0.25 s`; one radial finite-volume model with 20 nodes; one process/thread.
- Rest pools: ATP `8.2 mM`, PCr `32 mM`; initial ion disturbance `I=0`.
- Normoxic control boundary: `13.3 kPa`; low-oxygen stress: `2 kPa`; radius sweep: `20, 25, 30, 40 micrometre`; work sweep: `0.30, 0.60, 0.90 mM/s`.
- Oxygen concentration uses dissolved-O2 conversion `0.030 mol m^-3 kPa^-1`, oxygen diffusivity `1.5e-9 m^2/s`, Robin capillary-transfer coefficient `1.0e-5 m/s`, Michaelis constant `0.05 mol m^-3`, and an explicit O2-cost coefficient. The ion-pump coefficient is `0.10 s^-1`. These are declared assumptions in the model parameter table, not measurements.
- The same initial state, time grid, and functional score are used for every condition. The radial profile is volume weighted; no subject-specific tuning is permitted.

## Frozen criteria and failure rules

A run is mechanically valid only if all states are finite, `0 <= ATP <= ATP_rest`, `0 <= PCr <= PCr_rest`, `0 <= I <= 1`, dissolved oxygen is nonnegative, and the bounded oxygen update is nonnegative. A mechanistic pass requires all of:
- increasing boundary oxygen at fixed work and radius does not increase `T90_F` or reduce `V_PCr0`;
- increasing work from `0.60` to `0.90 mM/s` at fixed oxygen and radius does not decrease `T90_F` or reduce `V_PCr0` in the pre-registered stress range;
- increasing radius at fixed oxygen and work does not improve either recovery quantity;
- the joint-corner contrast meets the numerical prediction above, using the 600-s lower bound for a censored functional threshold;
- the analytic zero-consumption limit in `test_model.py` matches the closed-form flat oxygen solution.
A failure is reported as a failure, not repaired by changing a threshold. The source-scale check is descriptive: baseline `V_PCr0` is compared with the published `5.2–41.2 mM/min` span, but no claim of calibration is made if it falls outside.

## Controls and falsification

- Null model: same ATP/PCr and ion equations, but `O2_factor=1`, no radius dependence, and no diffusion delay. A stress result that is no worse than this null model is treated as evidence that the added oxygen/geometry mechanism is unnecessary for that contrast.
- Sham control: zero work demand from a rested state must not create a PCr deficit or ion disturbance beyond numerical tolerance.
- Anoxia control: `pO2_boundary=0` must suppress oxidative ATP recovery in the model; if it does not, the oxygen coupling is rejected.
- A larger radius is not declared beneficial merely because it changes the score: the direction must follow the solved oxygen field and the frozen criterion.

## Required next resolution step

Use matched longitudinal measurements of boundary/perfusion oxygen, PCr or ATP, an ion-balance proxy, and a pre-registered functional outcome under at least low/normal oxygen crossed with low/high work and measured radius or cross-sectional area. Fit only hierarchical parameters, retain subject and fiber heterogeneity, and test the joint interaction against the null model. Until then, all numerical combinations are model predictions, not measured recovery differences.
