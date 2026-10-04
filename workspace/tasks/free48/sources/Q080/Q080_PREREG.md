# BT-HX-Q080 — preregistration

## Scope and build history

Q080 asks whether the spatial patch distribution of receptors changes the spatial activation distribution and local peaks during a controlled ligand pulse when total receptor number is held fixed. The first model is a deterministic, effective two-dimensional membrane model. It is a mechanism test, not a measurement or a receptor-specific biological claim.

**Builds on**

- `K01`, `K02`, and `K03` as named in `inputs/QUESTION.md`; the supplied inputs do not include their underlying scripts or data.
- The local job directory was checked first. The interrupted session contained only `BRIEF.md`, `inputs/`, `ALLOW_WEB`, and `agent.log`; no runnable predecessor was present.
- Public primary sources were checked before this preregistration: Briddon et al. (2004), Suzuki et al. (2005), and Calebiro et al. (2013), cited below.

**Not redone**

- No internal data, measured response data, or hidden BodyTwin data are used.
- No receptor subtype is claimed to have the modeled response. The model is a transferable GPCR/membrane first-principles scaffold.
- Any source value is an anchor for a parameter, not evidence that a particular cell or receptor has that value.
- The only comparison changes receptor patch distribution. Domain, patch locations, total receptor number, ligand pulse, binding/dissociation rates, diffusion coefficients, and grid are identical between the two primary cases.

## Hypothesis and predicted quantity

At fixed total receptor number, concentrating receptors into patches should increase the maximum local activation density during a finite ligand bolus relative to a spatially uniform distribution, while increasing spatial nonuniformity. Total accumulated response may be similar or lower because high-density patches consume the finite ligand pulse more rapidly; the primary claim does not require a particular sign for the total response.

Primary predicted quantity:

`peak_ratio = max_x[A_patched(x, t_pulse)] / max_x[A_uniform(x, t_pulse)]`

where `A = Rbar*c` is bound-receptor activation density, `c` is the bound fraction field, and `t_pulse` is the end of the controlled pulse. Secondary quantities are the coefficient of variation of `A`, peak-to-mean activation, total bound activation, and time-integrated activation. The primary directional prediction is `peak_ratio > 1` and `CV_patched > CV_uniform`; exact numbers are not preregistered.

## Reference values looked up

1. **Primary source:** Briddon SJ, Middleton RJ, Cordeaux Y, Flavin FM, Weinstein JA, George MW, Kellam B, Hill SJ. “Quantitative analysis of the formation and diffusion of A1-adenosine receptor-antagonist complexes in single living cells.” *PNAS* 101:4673–4678 (2004). DOI: `10.1073/pnas.0400420101`.
   - **Figure 6b:** fitted binding isotherm for the membrane-localized ligand/receptor component: `Kd = 33 nM`; fitted `Bmax = 75 nM`, reported in the text as approximately `55 receptors/µm²`.
   - **Figure 5 and text:** free XAC-BY630 diffusion `D = 2.51 × 10^-6 cm²/s = 251 µm²/s`; the membrane complex component has `τD2 = 17.2 ms` and `D = 9.1 × 10^-9 cm²/s = 0.91 µm²/s`.
   - **Use:** `Rbar = 55 receptors/µm²` and `Kd = 33 nM` are numerical anchors. The pulse field, kinetic rates, and domain are assumptions, not copied measurements.

2. **Primary source:** Suzuki K, Ritchie K, Kajikawa E, Fujiwara T, Kusumi A. “Rapid hop diffusion of a G-protein-coupled receptor in the plasma membrane as revealed by single-molecule techniques.” *Biophysical Journal* 88:3659–3680 (2005). DOI: `10.1529/biophysj.104.048538`.
   - **Figure 4d and Conclusions:** smaller membrane compartment `L = 210 nm`; median residency time `45 ms`; median `D(25 µs, 100 ms) = 0.23 µm²/s`; the paper estimates the within-compartment free diffusion coefficient as `4.5–6 µm²/s`.
   - **Use:** This supports treating membrane patch scale and receptor mobility as separate timescales. It does not set the simulated patch width as a measured value.

3. **Primary source:** Calebiro D, Rieken F, Wagner J, et al. “Single-molecule analysis of fluorescently labeled G-protein-coupled receptors reveals complexes with distinct dynamics and organization.” *PNAS* 110:743–748 (2013). DOI: `10.1073/pnas.1205798110`.
   - **Figure 2D/I:** median lateral diffusion coefficients of `0.052 µm²/s` for β1AR and `0.039 µm²/s` for β2AR under the reported single-molecule conditions.
   - **Use:** `D_R = 0.04 µm²/s` is a conservative literature-anchored default for the relative receptor-state field, not a universal GPCR constant.

## Frozen model protocol

- Square membrane domain: `4.0 µm × 4.0 µm`, 21 × 21 finite-volume cells, `dx = 0.2 µm`.
- Mean total receptor density: `Rbar = 55 receptors/µm²`, giving `N_total = 880 receptors` in the domain.
- Primary patterns:
  - `uniform`: relative receptor density is exactly one everywhere.
  - `patched`: four fixed Gaussian patches centered at `(1.2,1.2)`, `(2.8,1.2)`, `(1.2,2.8)`, and `(2.8,2.8) µm`; default `sigma_patch = 0.30 µm`; the field is renormalized discretely so its integral is exactly `N_total`.
- Controlled pulse: a spatially localized Gaussian source centered at `(2.0,2.0) µm`, `sigma_ligand = 0.35 µm`, active for `0.10 s`, then switched off. Default source rate is `15000 nM/s`; this is an explicit model assumption and is not measured data.
- Ligand field: `L` in nM; default `D_L = 25 µm²/s` effective 2D value and `k_clear = 0.5 s^-1`.
- Binding: `k_on = 0.02 s^-1 nM^-1`; `k_off = 0.66 s^-1`; their ratio is the frozen `Kd = 33 nM` anchor. The rate pair is an assumption, not a fitted measurement.
- Receptor mobility: `D_R = 0.04 µm²/s` for free and bound relative receptor-state fields. Bound-state mobility is not given a separate value in the core model.
- Simulation endpoint: `1.2 s`; explicit no-flux finite differences with a conservative stability step.

## Frozen first-principles equations

Let `r(x,y,t)` and `c(x,y,t)` be dimensionless free and bound relative receptor-state fields, `L(x,y,t)` the nM ligand field, and `Rbar` the receptor-density anchor. Define `R(x,y,t)=Rbar*(r+c)` and `A(x,y,t)=Rbar*c`.

\[
\partial_t c = k_{on} L r-k_{off}c+D_R\nabla^2c
\]

\[
\partial_t r = -k_{on} L r+k_{off}c+D_R\nabla^2r
\]

\[
\partial_t L = D_L\nabla^2L-k_{clear}L+J(t)\exp[-(x-x_0)^2-(y-y_0)^2)/(2\sigma_{ligand}^2)]
\]

The source is present only during the pulse. The conserved receptor count is `∫R dA`; the conserved quantity in the implementation is `∫(r+c)dA`. At constant ligand and zero diffusion, the analytic binding limit is

\[
c^* = \frac{k_{on}L}{k_{off}+k_{on}L},\qquad K_d=\frac{k_{off}}{k_{on}}.
\]

The implementation must preserve these equations and must report receptor-mass error, ligand non-negativity, and unit dimensions.

## Frozen criteria and controls

**Primary pass criterion:** at pulse end, `peak_ratio >= 1.20` and `CV_patched - CV_uniform >= 0.10`. This is a model criterion, not a claim about biology.

**Interpretation bands:** `peak_ratio >= 1.20` with the CV condition is a mechanistic prediction supported by the run; a value between `1.05` and `1.20`, or a CV increment below `0.10`, is an inconclusive small effect; `peak_ratio < 1.05` is a failed directional prediction. In every case, report total response and uncertainty-free deterministic limitations rather than retuning the protocol.

**Mass and validity checks:**

- Initial and final total receptor count must agree within relative error `1e-10`.
- Free ligand and bound receptor fields must remain non-negative within `1e-10` nM or normalized units.
- The analytic `D_L=0`, `D_R=0`, constant-ligand limit must equal `c*` cell by cell within `1e-10`.
- A zero-pulse control must produce zero activation and a no-pulse control must remain at zero.
- The unit checker must report consistent dimensions for `L`, `Rbar`, `D_L`, `D_R`, `k_on`, `k_off`, source, `dt`, and `dx`.

**Countermeasures / placebo:** the uniform pattern is the same-total spatial placebo; the zero-pulse run is the temporal placebo; the zero-diffusion equilibrium run is the analytic negative control; a mass-reset check catches accidental pattern renormalization. No threshold, parameter, or criterion may be changed after seeing the primary result.

## Frozen sensitivity analysis

Run the patched-versus-uniform comparison at `−50%` and `+50%` for exactly three controlling parameters, all else held fixed:

1. `D_L`: effective ligand diffusion, default `25 µm²/s`.
2. `k_off`: receptor dissociation rate, default `0.66 s^-1`.
3. `sigma_patch`: patch width, default `0.30 µm`.

Report `peak_ratio`, CV increment, total-response ratio, and peak timing for every sensitivity point. These are deterministic structural sensitivities, not confidence intervals.

## What counts as an error

- A missing or non-primary source for a claimed reference value; any number presented as measured without a source/derivation/assumption label.
- Violation of any frozen criterion, receptor conservation, non-negativity, equation, or unit check.
- Changing the pattern normalization, pulse dose, total receptor count, or threshold after inspecting results.
- Treating the model output as experimental data, a validated receptor-specific parameter set, or evidence about a clinical/biological population.

## Next resolution step

Use experimentally measured, same-cell receptor coordinates and a calibrated ligand concentration/dose field; fit the patch width, effective ligand diffusion, and kinetic rates with preregistered holdout cells. Resolve ligand layer thickness and receptor–receptor cooperativity before interpreting a quantitative biological effect.
