# PREREG — BT-HX-Q090

## Freezing and scope

This is a first executable, quasi-static network model of a molar-loaded jaw. It is frozen before the first `model.py` run. The model predicts no unique individual's measurement. No internet check or data collection is performed.

## Builds on

- `inputs/QUESTION.md`: core mechanism, control case, outputs and falsification requirements.
- `inputs/NIGHT_PREAMBLE.md`: requirements for preregistration, model run, sensitivity and prohibited fabricated measured data.
- `BRIEF.md`: format, boundaries and deliverables.

## Not redone

- No direct code anchor or verified individual geometry exists in the task.
- No weighted model, network tolerance, damage law or validated tissue parameter is copied in.
- The literature comparison is used only as a rough external ground truth; it is not a substitute for measured data.

## Hypothesis and mechanism

Greater compliance in the tooth–PDL–bone path should transfer more of the same general chewing load to the TMJ paths and thus reduce PDL load per support area. For an idealised symmetric jaw with two molar gaps and two bilateral TMJ paths, in the linear limit

`phi_tooth = C_total / (2*C_tooth + 2*C_TMJ)`, where `C = 1/k` and `C_total` is the lateral compliance of the entire parallel connection.

The PDL load distribution between horizontal and inclined fibres follows `E_axial = Σ(f_i E_i cos²(alpha_i))`; buckled fibres cannot carry compression. A viscous Darcy term provides a kinetic complement and no irreversible damage is modelled.

## Primary ground truth and frozen references

- Primary mechanical output: `phi_tooth = (F_tooth_L + F_tooth_R)/F_close` at `F_close=100 N` and zero external moment.
- Secondary outputs: `delta_tooth_um`, `p_PDL_kPa`, `F_TMJ_L_N`, `F_TMJ_R_N` and `theta_TMJ_deg`.
- Control prediction: displacement for a single 100 N-static tooth path, compared with an indirect mobility anchor.
- Published comparison: **Bien, A. J. & Topp, H. A. (1970), “Changes in the periodicity of the tooth mobility pattern during mastication”, Journal of Periodology — UNVERIFIED, FROM MEMORY.** The remembered value is `delta_ref = 0.08 mm` under nominal `100 N` load. Load protocol and measurement precision have not been checked; external agreement is therefore `UNKNOWN` even if an arithmetic interval test is performed.

## Frozen acceptance criterion

1. **Arithmetic reference ground truth:** `0.5 <= delta_model_100N/delta_ref <= 2.0`, corresponding to `[0.04, 0.16] mm` for the memory anchor. This can only be reported as `PROXY_PASS` or `PROXY_FAIL`; external validity becomes `UNKNOWN` because of the unverified source and unknown load protocol.
2. **Load equilibrium:** the sum of vertical forces and net moment sum below `1e-6` times the magnitude of the respective input.
3. **Mechanistic ground truth:** increased PDL compliance reduces `phi_tooth` and increases the TMJ share; a shared scaling of all paths' compliance must not change the load shares.
4. **Analytical limiting case:** the linear global solution must match `K q = F` within `1e-10` relative.
5. **Sensitivity:** for `E_PDL`, `t_PDL` and `E_TM disc`, `+/-50 %` is reported one parameter at a time.

## Null model, placebo and falsification

- Null model: a serial tooth path alone cannot change the load distribution; the code path verifies that such a result is not incorrectly reported as a network distribution.
- Placebo: if all support compliances are scaled equally, `phi_tooth`, `phi_TMJ` and the ratio between the forces must remain unchanged within `1e-10`.
- The hypothesis is falsified for this output if reduced PDL compliance does not increase the load share in the TMJ paths, if load equilibrium is not solved, or if the analytical limit fails.
- A result cannot "succeed" by selecting new material data after the run; all assumptions and uncertainties are reported.

## Frysta modellantaganden

- Symmetric, idealised immediate load onset; no inertia in primary load distribution.
- Molar support: fibrous compressible PDL in series with compressible alveolar bone.
- TMJ: compressible layers (disc, condylar and temporal-bone cartilage) in series; no tensile opening or open-joint status in the primary run.
- `F_close=100 N` is a normalising scenario, not a claimed measurement.
- Quasi-statics is accepted only after the transient's acceleration flow has been checked against the load in a separate 1D step test.
- Irreversible damage and fracture are omitted; the work is reported, but no damage levels or data are invented.
