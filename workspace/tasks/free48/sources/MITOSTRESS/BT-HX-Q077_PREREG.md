# BT-HX-Q077 preregistration

## Scope and build-on

This run addresses the core mechanism in `inputs/QUESTION.md`: topologically distinct crista geometries are compared while inner-membrane area, intracristal-space volume, ATP-synthase and ANT site densities, and kinetic constants are held fixed. The model resolves ADP in the intermembrane space (IMS), ADP/ATP exchange by ANT, ATP synthesis in the matrix, ADP gradients, and diffusion and junction transport times.

`BRIEF.md`, `inputs/NIGHT_PREAMBLE.md`, and `inputs/QUESTION.md` were read. No prior implementation or result file for Q077 existed in this task directory at inspection time. The referenced repository anchor graph and source repositories were not written to or used as data. The only external anchor is the public primary research article below.

### Builds on

- `BRIEF.md`: bounded deliverable, no invented measurements, preregistration before execution.
- `inputs/NIGHT_PREAMBLE.md` §§0–7: first-principles scope, one-thread/local execution, source-vs-derivation separation, and next-step reporting.
- `inputs/QUESTION.md`: Q077 observables: ATP flux, ADP gradients, and transport times under matched area/volume and matched transporter numbers.
- Literature node: Adams, Afzal, Jafri & Mannella, *Cells* 14, 257 (2025), DOI `10.3390/cells14040257`, PMCID `PMC11853683`.
- Node-id: `UNKNOWN` for an internal BodyTwin node; no Q077 node id was present in this directory and none is inferred.

### Not redone


## Hypotheses and primary prediction

H1: With the same total membrane area and IMS volume, a topology that provides access at both ends of a long crista will have a smaller standing ADP depletion and a higher steady-state ATP-synthase output than a one-ended topology with the same total junction conductance.

H2: A branched topology whose terminal paths are shorter than the long trunk will reduce the ADP gradient and raise ATP output relative to the one-ended straight topology, even when the total membrane area and total number of transporter sites are unchanged.

H3: the maximum end-to-end diffusion time will track the longest graph distance to a cytosolic junction, so adding an opposite junction or short branches will lower the modeled transport time.

The preregistered primary quantity is:

`R_flux = J_ATP(two_trans_CJ) / J_ATP(one_CJ)`

at `L = 0.9 um`, with identical `A_IM`, `V_IMS`, ATP-synthase density, ANT density, and kinetic constants. Secondary quantities are maximum and mean fractional ADP depletion, matrix ADP, ATP output in mol/s and molecules/s, and diffusion/junction transport times.

## Reference value and source

Status: **VERIFIED**, not `UNVERIFIED`.

Primary source: Adams R, Afzal N, Jafri MS, Mannella CA. “How the Topology of the Mitochondrial Inner Membrane Modulates ATP Production.” *Cells*. 2025;14(4):257. DOI: `10.3390/cells14040257`. PMCID: `PMC11853683`.

- Figure 4 and Section 3.1.2 report long-crista ATP-synthase output of approximately `0.60 × JMAX` for detached or nearly bottlenecked long cristae and approximately `0.90 × JMAX` at maximum connectivity. The quantity is dimensionless relative flux, not a direct measurement of ADP concentration.
- The Discussion reports an approximately `17%` increase in relative ATP-synthase flux for the observed cardiomyocyte topology compared with one narrow CJ per crista.
- Figure 1 and Section 2.1 provide the modeling context: cytosolic ADP `0.037 mM` for the displayed map, matrix ATP/ADP starting values `0.44/0.72 mM`, and membrane potential magnitude `172 mV`.
- Section 2.1 reports a maximum reduced-model flux of approximately `100 molecules ATP/ms/um2`; this is used only as a scale check and is not imposed as a fitted target.

The model uses SI concentration units (`1 mM = 1 mol/m3`) and reports both SI and molecule-based flux. The literature value is an external plausibility anchor; matching it is not a success criterion.

## Frozen criterion

The primary prediction passes only if all of the following hold in the default run:

1. `R_flux >= 1.10`.
2. `max_fractional_ADP_depletion(two_trans_CJ) < max_fractional_ADP_depletion(one_CJ)`.
3. `max_diffusion_time(two_trans_CJ) <= 0.30 × max_diffusion_time(one_CJ)`.
4. The total membrane area, IMS volume, and ANT and ATP-synthase site counts agree between the compared topologies to numerical tolerance.
5. The analytical unit checks and the no-sink/equal-boundary limit test pass.

A topology effect is called unresolved if condition 1 fails, if the no-reaction limit develops a gradient, or if the matched-geometry bookkeeping fails. A negative result is retained; no parameter is retuned to force a pass.

## Counter-test and placebo

- No-reaction limit: with zero ANT/synthase sink and equal cytosolic boundary concentrations, the IMS ADP profile must be spatially uniform.
- Fast-diffusion placebo: increasing `D` by 50% must not increase the topology contrast beyond the fixed numerical tolerance used in the sensitivity report; convergence toward the no-crista limit is expected.
- A no-crista/flat-membrane reference is computed from the same total area and volume and is used only as a flux scale, not as a measured baseline.

## Fixed numerical scope

- One-dimensional finite-volume diffusion on graph edges with a local ADP sink.
- One fixed cytosolic ADP boundary and a Robin exchange conductance at each CJ.
- Matrix ADP is solved from steady-state ANT/ATP-synthase stoichiometric balance.
- Total adenylate pool is fixed; ATP and ADP are exchanged with one-for-one stoichiometry.
- Three geometries are run: one-ended straight, two-ended straight, and a branched tree.
- Sensitivity is run at `-50%`, baseline, and `+50%` for `D`, ANT turnover, and total CJ area. No Monte Carlo or parameter fitting is part of this bounded run.

## Model equations to implement

1. `D d2C/dx2 - q_ADP(C) = 0` in each IMS branch.
2. In the frozen low-substrate mass-action branch, `J_ANT = k_ANT K_ref (e u - e^-1 v)`, where `u=[ATP]_m[ADP]_IMS/(K_ATP K_ADP)` and `v=[ATP]_cyt[ADP]_m/(K_ATP K_ADP)`; `e=exp(2 fP F |psi|/(RT))` for ADP-in/ATP-out exchange. The saturation denominator is deliberately omitted rather than fitted.
3. `J_AS = V_AS,pmf ( [ADP]_m/K_ADP [Pi]_m/K_Pi )/(1+[ADP]_m/K_ADP+[ATP]_m/K_ATP+[Pi]_m/K_Pi)`.
4. `sum_i J_ANT,i A_i = J_AS A_IM` at steady state.
5. `tau_diff = max_x d_graph(x,CJ)^2/D`; `tau_CJ = V_IMS/G_CJ`.

All concentrations are in `mol/m3`, lengths in `m`, areas in `m2`, volumes in `m3`, diffusivity in `m2/s`, kinetic coefficients in `s^-1` or `m/s`, and flux in `mol/(m2 s)` unless explicitly labeled total flux.
