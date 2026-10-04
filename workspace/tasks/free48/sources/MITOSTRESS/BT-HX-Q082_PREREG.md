# BT-HX-Q082 — preregistration

## Hypothesis and reference

An explicit moving-control volume model can move both neutral species and positive ion through a membrane without introducing a numerical dilution or a drift term. Reference is the relative deviation between the model concentrations and an analytical limiting case for a uniformly initially distributed solution when only the geometry of the membrane is changed. Positivity, balance residual for substance and charge, and surface-bound number are recorded at each step.

## Predicted quantity and frozen criterion

The primary prediction is that proper volume displacement keeps a uniform concentration constant to machine absolute precision. The following criteria are frozen before the first run:

1. `max_relative_concentration_error <= 1e-10` in the analytical limiting case.
2. `max_normalized_mass_residual <= 1e-10` when all internal flows are closed.
3. `max_normalized_charge_residual <= 1e-10` when all internal flows are closed.
4. `min_inventory >= 0`, `min_concentration >= 0`, `min_surface_bound_count >= 0` i huvudscenariot.
5. With a leftward concentration gradient, both the net transfer of neutral and positive ion must be greater than zero at the end time, and the surface-bound number must be greater than zero.
6. The placebo that changes membrane volume without bringing with it the displaced solution volume substance must fail criterion 1; otherwise the test is trivial.
7. A thirty-factor sensitivity with the factors 0,5, 1,0 and 1,5 shall account for all runs. A result with a negative inventory, concentration or surface-bound number is counted as an error even if the sum were correct.

What counts as an error: `NaN`/`Inf`, negative inventory, broken volume balance, failure to report external flow, or placebo 6 passing the criterion 1. A falling criterion is a result, not something that can be hidden or reinterpreted after the run.

## The frozen core of the model

The model has two fluid domains with concentrations `c_L` and `c_R`, a moving membrane mode `x_m`, and separate substance inventories `N_s^L`, `N_s^R`, `N_s^M` and ion inventories `N_+^L`, `N_+^R`, `N_+^M`. `M` is a surface bound layer. For a rigid outer domain, `V_L=A_r x_m` and `V_R=A_r(L-x_m)` are

When `x_m` increases with `dx`, the geometric volume flow is `Q_delta=A_r dx`. For each solute, `Q_delta c_R` is then transported from right to left, as `V_L` increases as `x_m` increases. Without this volume flow record, changing control volumes alone would produce an artificial concentration change. The neutral matter flux is `J_s=k_s A_m(c_L-c_R)`. For a positive monovalent ion is `J_+=k_+ A_m[c_L exp(-alpha F V_m/RT)-c_R exp(alpha F V_m/RT)]`.

The surface mass `A_m` follows a frozen deformation `epsilon_m(t)`. The capacitance is `C_m=c_A A_m`, where `c_A` is specific capacitance. The surface bound charge is used as `V_m=(q_b A_m+F N_+^M)/C_m`, with `q_b` as an assumed fixed surface charge bias. In this way, the voltage changes when the surface is deformed or the number of bound ions changes. The total charge of the ions is `F(N_+^L+N_+^R+N_+^M)`; a fixed countercharge in the bulk is not included in this first model and is an explicit limitation.

Binding to the membrane is a neat Langmuir-like process `r_s=k_on,s A_m c_s^L(1-N_s^M/N_max,s)-k_off,s N_s^M`, and the corresponding process for ions. It is positive when the substance binds to the surface and negative when it is released. Fringe fluxes are donor cell fluxes with specified inflow concenration; total inflow and outflow must balance the volume.

## Reference value

Primary source: Pérez-Mitta, G. & MacKinnon, R. (2023), “Freestanding lipid bilayer tensiometer for the study of mechanosensitive ion channels”, *PNAS* 120(12):e2221541120, DOI `10.1073/pnas.2221541120`. In Fig. 3B (specific capacitance versus pressure, six bilamellae) the interval `0.3–0.5 microF/cm^2` is reported; the article also states that capacitance per area is constant within about 10 % during printing. The value is thus **ACQUIRED**, not a made-up measurement. The model uses the center of the interval `c_A=0.4 microF/cm^2 = 0.004 F/m^2`; this is an assumed pre-representative value within published range.

## Parameters that are not measurements

`L`, `A_r`, `x_m(0)`, `x_m(T)`, `k_s`, `k_+`, `alpha`, `T`, `q_b`, binding rates, binding capacities, and initial concentrations are explicit model assumptions for a controlled mechanistic test run. They are not claimed metrics. `F=96485.33212 C/mol` and `R=8.314462618 J/(mol K)` are used as physical constants.

## Building on

- `inputs/QUESTION.md`: Q082's functionality surface with balance residual, positivity, concentration, surface bound number, deformation, edge flux and capacitance.
- `inputs/NIGHT_PREAMBLE.md`: requirements for pre-registration, unit control, conservative diagnoses and that sensitivity is reported.
- Published primary source: Pérez-Mitta & MacKinnon (2023), DOI `10.1073/pnas.2221541120`, Fig. 3B.
- Working directory's previous break: no `model.py`, `test_model.py`, `PREREG.md`, or `results.json` existed to build on.

## Not redone

No external musculoskeletal solver run, no body model, no invented measurement, no batch/LOSO-svep and no hidden data directory is used. The full Poisson–Nernst–Planck model with counterions, ion channels, three-dimensional deformation, and measurable surface tension is not claimed to be solved in this first executable step.

## Frysning

The criteria, parameter limits and source are established in this text before running `model.py`. The sensitivity factors are exactly 0,5, 1,0 and 1,5.
