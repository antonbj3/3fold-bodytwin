# BT-HX-Q049 — frozen preregistration

**Status:** frozen before the first model run, 2026-09-25. This is a bounded first mechanistic model, not a clinical or tribometer calibration.

## Hypothesis and predicted quantities

The core hypothesis is that a mechanically correct, transport-capable bearing interface requires at least two *functionally distinct* thin layers to remain explicit: (1) an interfacial lubricating/transport-conditioning film and (2) a depth-specific superficial zone. A third thin layer at the cartilage–bone boundary is predicted to control fluid exchange with the substrate. The model tests this rather than assuming that every visible layer is necessary.

The primary predicted quantity is the short-time surface fluid-load support

`psi_s(t) = p_surface(t) / sigma_total(t)`.

The secondary quantities are the top exudation flux, bone-side exchange flux, total compressive stress, solid stress at the surface, and the time for the surface pressure to fall halfway from its post-step value. The decision is based on changes in both mechanics (`psi_s`, solid stress, total stress) and transport (fluxes and drainage time), not on a friction number alone.

Predictions registered before implementation:

- **P1:** retaining a finite interfacial film will reduce the rate of pressure loss and preserve higher `psi_s` than replacing it by a zero-thickness, freely draining interface.
- **P2:** replacing the superficial zone by bulk properties will change `psi_s` or the pressure gradient even if total tissue thickness is held fixed, because the surface zone has a distinct modulus and permeability.
- **P3:** removing the calcified/barrier layer will increase bone-side fluid exchange and shorten the pressure-drainage time.
- **Null/placebo:** if a homogeneous, single-layer control reproduces the full-layer outputs within the frozen tolerances, the corresponding layer is not resolved as necessary by this first model.

## Scope and facies

The model is a one-dimensional, small-strain, layerwise consolidation column normal to the bearing surface. The column contains, from the free surface toward the substrate, an interfacial film, a superficial layer, bulk cartilage, and a calcified/barrier layer. The film is represented by a series mechanical compliance and a hydraulic boundary conductance; it is not treated as a resolved porous control volume. The substrate is a rigid, pressure-sink boundary. No lateral contact geometry, cyclic waveform, chemistry, friction law, or tissue damage law is claimed.

The fixed numerical case is a stress step normalized to 5% of the bulk constrained modulus, with observation at 1 s, 10 s, and 60 s. The normalization is a modeling choice within a small-strain regime, not a patient-specific load.

## First-principles equations

For each resolved porous layer, with compressive solid strain `e_i` positive, fluid pressure `p_i`, and hydraulic conductivity `K_i = k_i / eta_i`, the reduced mass balance after an instantaneous total-stress step is

`S_i dp_i/dt = -dJ_i/dz`,
`J_i = -K_i dp_i/dz`.

Here `M_i` is the constrained solid modulus, `S_i = alpha^2 / M_i` is the reduced storage coefficient, `k_i` is intrinsic permeability, and `eta` is viscosity. The coefficient `alpha=1` is an assumption for this first run. The step initializes the near-incompressible confined state at `p_i(0+) = sigma_total`; the subsequent pressure loss is the modeled transport phase. The top and bottom boundary fluxes are resistance laws, with ambient pressure set to zero:

`J_top = -G_top p_top`,
`J_bone = G_bone p_bone`.

The total compressive stress is prescribed and common to the series layers. The solid stress and strain in each layer are

`sigma_solid,i = sigma_total - p_i`,
`e_i = sigma_solid,i / M_i`.

Thus `psi_s = p_surface / sigma_total` starts at one in the immediate analytical limit and decreases only through drainage. Flux signs are positive downward, from the free surface toward bone.

The interfacial film contributes a series compliance `h_f / M_f` and a hydraulic boundary conductance `G_top`; removing the film sets its thickness and modulus to zero and substitutes the preregistered direct-drain conductance. The barrier contribution is resolved through its finite thickness, stiffness, and permeability, followed by the bone-side conductance.

## Literature anchors (verified before implementation)

1. **Park, Krishnan, Nicoll & Ateshian (2003), “Cartilage interstitial fluid load support in unconfined compression,” Journal of Biomechanics 36:1785–1796, DOI `10.1016/S0021-9290(03)00231-8`.** Primary source. Table 1 reports peak `dWp/dW` of **94 ± 4% at the articular surface** and **71 ± 8% near the deep zone** for bovine cartilage (dimensionless fluid-load support). The same table reports 79 ± 11% and 69 ± 15% for human cartilage. These values are an anchor, not a direct prediction for this 1-D column.
2. **Krishnan, Park, Eckstein & Ateshian (2003), “Inhomogeneous Cartilage Properties Enhance Superficial Interstitial Fluid Support and Frictional Properties, But Do Not Provide a Homogeneous State of Stress,” Journal of Biomechanical Engineering 125:569–577, DOI `10.1115/1.1610018`.** Primary source. Table 1 gives representative femur depth-dependent `H−A` values of 0.32–0.73 MPa and `kz` values of 0.39–0.58 × 10^-15 m^4/(N·s); Table 2 gives contact `Wp/W` values of 0.93–0.98 at 1 s. These are used as parameter anchors/proxies, not as a claim that the tested slices equal the model layers.
3. **Flannery et al. (1999), DOI `10.1006/bbrc.1998.0104`.** Primary source establishing the lubricating surface-associated protein context. It is not used to assign a film thickness or modulus.

The numerical reference quantities above were located in the primary articles; no `OVERIFIERAD` tag is used. Film thickness, film modulus, direct-drain conductance, barrier properties, viscosity, and the 5% load are explicitly assumptions rather than disguised measurements.

## Frozen parameter policy

The run uses SI units. The superficial modulus and permeability are the femur values from Krishnan Table 1 as a representative proxy; bulk values use the deep-zone femur values. The surface and barrier thicknesses are deliberately simple frozen values so the layer effect is identifiable. `G_top`, `G_bone`, `eta`, and `alpha` are assumptions. No internal BodyTwin/external musculoskeletal solver data are used. The parameter table in `model.py` repeats units and provenance for every numeric parameter.

## Frozen comparisons and criteria

The following are computed before any qualitative interpretation:

- **Immediate analytical limit:** for a homogeneous, sealed column, an instantaneous stress step has `p=sigma_total`, zero flux, and `psi_s=1`; its solid stress is zero at the instant of the step. Numerical error must be below 1%.
- **Conservation:** cumulative fluid-content change must agree with integrated boundary flux within 1% for a sealed test.
- **Resolution:** outputs at 16 and 32 cells per resolved layer must agree within 1% for `psi_s(1 s)` and 5% for the two boundary fluxes.
- **Full-layer plausibility anchor:** `psi_s(1 s)` for the frozen full case must lie in **[0.90, 1.00]**. This is intentionally a broad interval around the verified 0.93–0.98 contact-model range and is not a clinical validation threshold.
- **Superficial-layer necessity:** call the superficial layer resolved as necessary only if its explicit-versus-merged comparison changes `psi_s(1 s)` by at least **0.02 absolute** or changes `J_bone(10 s)` by at least **25% relative**.
- **Film necessity:** call the film transport-relevant only if explicit-versus-no-film changes `psi_s(1 s)` by at least **0.02 absolute** or `J_top(1 s)` by at least **25% relative**. This is a model-structural conclusion because film properties are assumed.
- **Barrier necessity:** call the barrier transport-relevant only if explicit-versus-no-barrier changes `J_bone(10 s)` by at least **25% relative** or `psi_s(10 s)` by at least **0.02 absolute**.
- **Sensitivity:** `k_bulk`, `M_superficial`, and `G_top` are varied individually by ±50%. A result is called robust only if the sign of the reported layer effect is unchanged across all six perturbations; effect sizes, not significance, are reported.

A layer that misses its threshold is reported as “not resolved as necessary by this model,” not as biologically absent. A failure of the analytical, conservation, resolution, or plausibility criterion invalidates the run as a mechanism test.

## Counterfactual and stopping rule

The zero-model/placebo is a homogeneous porous column with the same total tissue thickness, mean constrained modulus, mean permeability, and no explicit film or barrier. It is run through the identical solver and compared before layer attribution. The run stops after the frozen scenarios, convergence checks, and three one-factor sensitivity sweeps; it does not tune parameters to force a layer to pass.

## Builds on

- `BRIEF.md`, `inputs/QUESTION.md`, and `inputs/NIGHT_PREAMBLE.md` in this directory.
- The primary literature anchors listed above; no node IDs or prior BT-HX result files are present in the bounded directory.
- The interrupted-session log contains only directory reads and no model output to reuse.

## Not redone

No external solver runtime, internal data, external repository writes, literature values not listed above, or unverified measurements are used. The result is a first 1-D mechanistic model; its conclusion is conditional on the stated boundary conditions and assumed film/barrier parameters.
