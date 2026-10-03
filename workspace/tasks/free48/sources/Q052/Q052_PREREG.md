# BT-HX-Q052 — fryst preregistration

Status: `FROZEN` before the first run. No parameter or model code changes after the hash is created.

## Scope and goal

Q052 asks for a first executable mechanistic model of how process-dependent porosity couples to stiffness, permeability and tissue response for the same sample. The model is a single 1D representative unit for a trabecular scaffold sample. It is deliberately a mechanistic baseline model, not a clinical validation.

## Builds on

- Node `Q052` in `inputs/QUESTION.md`; the requirement for the same microstructure, mechanics and transport.
- Local instructions and requirements in `BRIEF.md` and `inputs/NIGHT_PREAMBLE.md`.
- Public primary study: Chao et al. (2021), *Frontiers in Bioengineering and Biotechnology* 9:779854, DOI `10.3389/fbioe.2021.779854`. Figure 13A verifies `K = 1.87e-8 m²` at `φ = 0.80` and `D = 800 µm`; figure 14A-B gives natural bone permeability `K = 1.50e-10 m²`; figure 15 shows the direction in the cell attachment experiment.
- Public mechanistic background: Karande, Ong & Agrawal (2004), DOI `10.1007/s10439-004-7825-2`, for the coupling between pore geometry and nutrient diffusion; Alias & Buenzli (2018), DOI `10.1007/s10237-018-1031-x`, for pore geometry control of osteoblast behavior.
- No internal data and no assumed measurements are used.

## Not redone

- No body movement model, training data, animal data or patient data is used.
- No claimed validation of the tissue response is made; the response law is a hypothesis.
- The process variable is an explicit dimensionless AM process mode, not a calibrated machine model.
- External BodyTwin anchors according to the preamble (`MECHANISM_ANCHOR_GRAPH.json`, `scripts/msk/`, `docs/MECHANISM_*`, `bt_memory/`) were not readable in this session and are therefore not reused.

## Hypothesis and reference

**H1:** For a given process sequence, a common microstructure `s → {φ, d}` arises. Increased pore volume and pore size increase hydraulic and diffusive transport, but decrease load-bearing solid volume and can therefore give a non-monotonic response to mechanical loading and tissue response.

Facit: `Q052-porositet–mekanik–transport–respons`.

## Frozen prediction target

The primary prediction is hydraulic permeability `K` for the nominal case `s = 0.5`, giving `φ = 0.80` and `d = 800 µm`.

The reference value is **VERIFIED** in Chao et al. (2021), figure 13A: `K_ref = 1.87e-8 m²`, at `φ = 0.80`, `d = 800 µm`, Ti-6Al-4V scaffold, water at 37 °C. This is a measured/reproduced CFD result in a primary study, not an internal BodyTwin value. The natural bone reference `1.50e-10 m²` (figure 14) is context, not a requirement on the same geometry.

## Frysta ekvationer

1. Processgeometri:
   `φ(s) = φ_min + (φ_max - φ_min)s`
   `d(s) = d_min + (d_max - d_min)s^q`.
2. Pore geometry and tortuosity:
   `S_v = 6φ/d` and `τ = 1 + b_τ(1-φ)`.
3. Permeability, Darcy/Kozeny-Carman form:
   `K = K_g φ^3 d^2 / τ²`.
   `K_g = 0.0690521240234` dimensionless is a simple geometric calibration to the verified nominal reference; it is not an independent validation.
4. Mechanics:
   `E_eff = E_s χ_E (1-φ)^m_E (d_ref/d)^m_d`.
   `E_s = 110 GPa`, `χ_E = 0.75`, `m_E = 2.0`, `m_d = 0.35`; this is a strut/network approximation, not a universal empirical curve.
5. Darcy velocity and wall shear:
   `u = K g_p/µ`; `τ_w = µu/(d/2)`, where `g_p` is the pressure gradient.
6. Diffusion and advection:
   `D_eff = D_0 φ/τ²`.
   Steady-state physical field: `D_eff c'' - u c' - λc = 0`, `c(0)=c_in`, `c(L)=c_out`.
7. Hypothetical, dimensionless response potential:
   `R = response_max · f_O(c_mid) · f_E(ε_loc) · f_S(τ_w)`,
   where `f_O=c_mid/(c_50+c_mid)`, `f_E=exp[-0.5((ε_loc-ε_peak)/ε_width)^2]`, and `f_S` is a log-normal shear window. `R` is not measured data and must not be interpreted as cell count.

## Fryst parameterregister

| Parameter | Value | Unit | Source or assumption |
|---|---:|---|---|
| `phi_min` | 0.70 | 1 | assumed process endpoint |
| `phi_max` | 0.90 | 1 | assumed process endpoint |
| `d_min` | 0.60e-3 | m | Chao et al. design level |
| `d_max` | 1.00e-3 | m | Chao et al. design level |
| `q` | 1.0 | 1 | assumed linear process coupling |
| `E_s` | 110e9 | Pa | Chao et al. material data |
| `chi_E` | 0.75 | 1 | assumed network correction |
| `m_E` | 2.0 | 1 | assumed bending scaling |
| `m_d` | 0.35 | 1 | assumed architecture effect |
| `d_ref` | 0.80e-3 | m | nominal reference |
| `b_tau` | 0.50 | 1 | assumed tortuosity model |
| `K_g` | 0.0690521240234 | 1 | calibrated to Chao Fig. 13A |
| `D_0` | 2.0e-9 | m²/s | assumed diffusive scale at 37 °C |
| `mu` | 1.45e-3 | Pa s | Chao et al.; converted from MPa·s |
| `L` | 12e-3 | m | Chao et al. sample height |
| `g_p` | 2.0e4 | Pa/m | assumed perfusion |
| `lambda` | 0.05 | 1/s | assumed dimensionless consumption rate |
| `c_in`, `c_out` | 1.0, 1.0 | 1 | normalized boundary conditions, assumed |
| `grid_points` | 101 | 1 | numerical choice |
| `sigma_app` | 2.0e6 | Pa | assumed compression load |
| `k_t` | 0.12 | 1 | assumed stress concentration |
| `epsilon_peak` | 0.003 | 1 | assumed mechanotransduction window |
| `epsilon_width` | 0.004 | 1 | assumed width |
| `c_50` | 0.5 | 1 | assumed normalized half-saturation |
| `tau_opt` | 0.5 | Pa | assumed shear window parameter |
| `tau_ratio_width` | 1.5 | 1 | assumed shear window width |
| `response_max` | 1.0 | 1 | hypothetical maximum potential |

## Frysta acceptanskriterier

The model is accepted for this first run only if:

1. All unit checks pass and all outputs are finite.
2. `|K_nominal/K_ref - 1| <= 0.25`.
3. Analytical limiting cases pass: `φ=0 ⇒ K=0`; without flow and consumption the field is uniform; at fixed `φ,τ`, `K∝d²` applies.
4. Nominal `E_eff` lies within a factor of 1.5 from Chao et al. figure 8/9's summarized interval `2.6–4.0 GPa`.
5. In a process sweep, the correlation between `s` and mean concentration is at least 0.80; a fixed random placebo coupling has median absolute correlation below 0.30.
6. No negative concentrations, permeabilities or response potentials.

The criteria are frozen here. A failure is reported as a result and is not changed retroactively.

## Countertests and error limits

- Zero porosity is a blocking transport limit.
- The uniform field without advection and consumption tests the diffusion limit.
- The `K(d)` test tests quadratic pore size scaling.
- Placebo: decouple K from the process sequence through random permutation of K across the same geometry samples and compare transport correlation.
- Failure to keep process geometry constant while s changes is an error.
- Negative K, D, concentration or R is an error.
- If the reference value or figure reference cannot be reproduced, it must be marked `UNVERIFIED`; here it is verified.

## Execution and provenance

`model.py` must expose the parameter register, `simulate`, `transport_field`, `sensitivity` and `run_analysis`. `test_model.py` must contain at least one analytical limiting case. `results.json` must separate `literature_reference`, `derived`, `assumption` and `hypothesis`, and `RESULTS.md` must start with `BT-HX-Q052`.
