# BT-HX-Q084 — preregistration

## Scope and question

This is a bounded, first-principles model of passive transport of a finite-size tracer in a crowded compartment. The primary estimand is the passage-time contrast between a spatially resolved crowding model and a matched scalar adjusted-diffusivity model, together with the bound/free spatial distribution.

The model has three scales: micron-scale excluded obstacles, nano-scale excluded obstacles with a hydrodynamic drag closure, and reversible binding. It is a mechanistic model, not a claim that any supplied biological or clinical measurements exist.

## Builds on

- `inputs/QUESTION.md`, Q084: the required representation of passage time, binding kinetics, spatial distribution, packing fraction, particle interactions, and matched geometric obstacles.
- `inputs/NIGHT_PREAMBLE.md`: read-only use, explicit separation of source versus derivation versus hypothesis, and the requirement to write only in this task directory.
- `inputs/NIGHT_PREAMBLE.md` §1: no prior task artifact was present in this directory. The referenced BodyTwin anchor graph, scripts, docs, and prior result maps were not available through the bounded workspace, so no node id or external implementation is reused.
- Public primary source checked before coding: Destrian, Moisan, Mège, Ladoux, Goyeau & Chabanon, “Cytoplasmic crowding acts as a porous medium reducing macromolecule diffusion,” *PNAS* 123(4), e2519599123 (2026), DOI `10.1073/pnas.2519599123`, PMCID `PMC12846847`. Fig. 3C reports iso-osmotic free-GFP diffusivity of approximately `16–20 µm² s⁻¹`; Fig. 3D reports about `20%` lower diffusivity in low-porosity regions than in high-porosity regions. The text accompanying Fig. 4D–G reports less than `10%` reduction from micro-obstacles and up to `70%` reduction from nano-obstacles. These are verified source values, not measurements generated here.

## Not redone


## Hypothesis and estimands

At a fixed mean excluded nano-volume fraction and identical binding rates, local crowding plus reversible binding changes first passage relative to a uniform adjusted diffusivity. The primary estimand is

`Rpass = tau_mechanistic / tau_adjusted - 1`.

Positive `Rpass` means the mechanistic model is slower. Secondary estimands are the effective diffusivity ratio `D_eff / D_alpha0`, the bound occupancy, and the spatial range/CV of the late-time free and bound occupancy. A uniform model can reproduce the mean local diffusivity but is expected to miss bottlenecks and spatially varying binding residence times.

## Frozen mechanism and parameters

The implementation will use a one-dimensional finite-volume Markov diffusion model with cell width `Δx`, reflecting inlet, absorbing outlet, and two states per nonterminal cell: free and reversibly bound.

- `D_alpha0 = 20.0 µm² s⁻¹`: assumed iso-osmotic free-cytosol diffusivity for the nominal tracer; it is an explicit model input, not a fit.
- `D_water = 87.0 µm² s⁻¹`: source reference for free GFP in water, used only as context.
- `L = 10.0 µm`, `N = 161`, tracer radius `R = 2.3 nm`, characteristic nano-pore radius `a = 4.0 nm`, micro excluded fraction `phi_m = 0.08`.
- `phi_n(x) = 0.18 + 0.08 cos(2πx/L)`: a deterministic high/low spatial obstacle field; all `phi_n` values remain below 0.26.
- Relative volume `V_r = 1.0`; solvent correction `D_alpha(x) = D_alpha0 / (1 + k_visc(1/V_r - 1))`, with `k_visc = 1.0`.
- For scale `i`, accessible fraction is `chi_i = 1 - phi_i`; the local path-length closure is `tau_i = 1 + 0.5 phi_i / chi_i`.
- Kozeny–Carman-like permeability is `K_nano = a² chi_n³ / [2 phi_n² + chi_n phi_n (1 + interaction_strength phi_n)]`, with `K_nano = ∞` at `phi_n = 0`; `interaction_strength = 4.0`. The particle hindrance is `H = 1 + R/sqrt(K_nano) + (R/sqrt(K_nano))²/9`, so the empty-space limit is exactly `H = 1`.
- The local free diffusivity is `D(x) = D_alpha(x) (chi_m/tau_m) (chi_n/tau_n) / H`.
- Binding is `free → bound` with `k_on(x) = k_on0 (1 + interaction_strength phi_n(x)) chi_n(x)` and `bound → free` with `k_off(x) = k_off0`; `k_on0 = 0.12 s⁻¹` and `k_off0 = 0.03 s⁻¹`. This is the minimal coarse-grained encounter/interaction closure.
- The adjusted-diffusivity comparator is `D_adj = spatial mean of D(x)` over the same transient (non-absorbing) grid cells, with the same `k_on(x)` and `k_off(x)`. Thus the comparison holds mean free mobility and chemistry fixed and isolates heterogeneity and binding-state transport.
- Observation time for the spatial distribution is `t_obs = 10.0 s`; passage time is the expected first exit time from the absorbing outlet, not a time to a fitted percentile.

## Frozen criteria and falsification rules

1. **Mechanistic separation:** pass if `abs(Rpass) >= 0.10` at the nominal field. Otherwise report that a scalar adjusted diffusivity is not distinguishable at the preregistered threshold.
2. **Spatial signature:** pass if the late-time bound occupancy range (max minus min over spatial bins) is at least `0.10`, or its spatial coefficient of variation is at least `0.10`. If neither occurs, the binding/spatial component is rejected as unresolved.
3. **Analytical sanity:** with no obstacles, no interaction, no binding, and a constant `D`, the expected first passage of the discrete walker must agree with `L²/(2D)` within `2%` for the limiting one-sided diffusion problem, and the no-obstacle mean diffusivity must equal its input within `1e-12` in the finite-difference construction.
4. A failed criterion is reported as a failure; no post hoc parameter tuning, measured data, or fabricated observation is permitted. Large numerical error, negative rates, inaccessible-volume fraction below zero, or a nonconvergent linear solve is an execution error, not scientific evidence.

## Reproducibility

The run must use one process, the frozen defaults above, deterministic grid construction, and a fixed random seed only for any optional diagnostic. Numerical outputs are written to `results.json`; the implementation and analytical test are in `model.py` and `test_model.py`.
