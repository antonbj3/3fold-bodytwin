# BT-HX-Q031 — PREREG

## Status and scope

This is a first executable mechanistic model, not a validated measurement model. It does not use internal BodyTwin, restricted model data or personal data and makes no parameter adaptation to new data. All measurement values below are literature anchors; all other values are explicit assumptions.

## Hypotes

When an osmotic perturbation changes intracellular osmolarity, membrane potential and ion transport must be solved together with water volume if any of the following is true: the electrical potential drives a measurable net ion flux, ion flux changes total intracellular osmolarity over the relevant timescale, or geometry changes membrane capacitance and conductance. A pure Nernst/water model is therefore insufficient when electrical and osmotic feedback exceed the frozen error scenario.

## Proposed model

The model is a spherical approximation with effective area/volume and three permeable monovalent ions (Na+, K+, Cl−), an impermeable negative proton anion, an electrogenic Na/K pump and water flow. The state is Vm, V, NNa, NK, NCl, water amount, three cumulative outgoing ion fluxes and a small uncontrolled countercharge reserve NΔ that maintains bulk electroneutrality. Internal charge and the capacitance voltage law are treated separately; the pump's net ion current is 3 Na out and 2 K in per pump cycle. Electrical and chemical potential enter the Nernst fluxes, and water flow is calculated from ideal osmotic pressure balance.

The model is a reduced first-principles core. It does not include ion activity coefficients, pump regulation, membrane elasticity, hydrostatic pressure or the cell type's full ion channels. These omissions are intentional and reported as limitations.

Before the first complete run, a unit check of the water law was performed. The published Pf unit cm/s is not treated as m/(s Pa). It is therefore refrozen to Lhyd = 6,1×10^-6 m/s / (150 mol m^-3 · R T) = 1,65×10^-11 m s^-1 Pa^-1, with Pf interpreted as a coefficient per 150 mM step. To represent a cell layer, a fixed hydraulic surface efficiency fw = 0,025 of the electrical geometric area is also used. This is an explicit layer geometry assumption, not a parameter that may be changed after seeing the result. Electrical area, capacitance and ion fluxes use the full area; only water flow uses fw.

## Referensankare

Primary source: Farinas J, Verkman AS. Cell volume and plasma membrane osmotic water permeability in epithelial cell layers measured by interferometry. Biophysical Journal. 1996;71(6):3511–3522. DOI: 10.1016/S0006-3495(96)79546-2. PMID 8968620, PMCID PMC1233838.

Verified values in the article:

- Figure 6 and text on p. 3517–3518: after a step from 300 to 150 mOsm, MDCK cells' relative volume increases approximately twofold; the time curve indicates a half-time of approximately 6 s at 24 °C.
- Results text on p. 3518: calculated plasma membrane water permeability Pf = 6,1 ± 2 × 10^-4 cm/s at 24 °C for MDCK. Figure 6's microscopic estimate is given as 9 × 10^-4 cm/s; this is a separate estimation method, not a mean to call it.
- Text on p. 3516: mean volume 1,7 pL and area/volume 6750 cm^-1 for the measured MDCK layer at 300 mOsm.

The reference value frozen here is the ideal osmotic limit V∞/V0 = 300/150 = 2,0 and the experimental order of magnitude t50 ≈ 6 s for this reference protocol. The values are not claims that the generic model describes MDCK exactly.

Structural support: Jakobsson E. Interactions of cell volume, membrane potential, and membrane transport parameters. American Journal of Physiology-Cell Physiology. 1980;238(5):C196–C206. DOI: 10.1152/ajpcell.1980.238.5.c196. PMID 7377338. The abstract describes a system with n transport equations plus bulk electroneutrality and isotonicity. Jagobsson is used here as structural support, not as a source for a new measurement.

## Frozen computational conditions

1. Initially, electrical potential is extrapolated to zero net cell current with the extrinsic 300 mOsm concentrations. At t = 0, extrinsic ion concentrations are set to 50 % of their original values, giving 150 mOsm. Internal osmolarity and all internal amounts are kept as the initial state.
2. The full model is compared with a placebo where electrogeneration, conductivity and the pump are set to zero and total intracellular osmolarity is held constant. The placebo must give V∞/V0 = 2,0 within 1 % in the ideal limit.
3. The electrical placebo/countertest model holds Vm at its initial value but retains the same ion and water fluxes. The primary prediction is χV = max_t |Vfull(t) − VVm-clamp(t)| / V0.
4. Coupling is considered necessary for this frozen protocol only if χV > 0,02. If χV ≤ 0,02, this operational placebo distinction is insufficient to require full dynamic coupling in this particular case.
5. The reference gate is that the full model lies within 1,8 ≤ V∞/V0 ≤ 2,2 and 3 s ≤ t50 ≤ 10 s. Passing the gate is a reproducibility requirement, not a claim of external validation.
6. The full model's numerical requirement is that the relative final closed amount and water balance is < 10^-6. The charge check is F[(NNa+NK−NCl−NX)+∫(JNa+JK−JCl)dt] against the initial value; it must be < 10^-6. Bulk electroneutrality including NΔ is reported separately. Capacitance charge is not incorrectly treated as an extra net exchange of ionic charge; its KCL is C_m A dVm/dt = −(I_p+I_ion).

## Countertests and placebo

- Zero-flux limit: each Nernst flux must be zero when Vm = Eion.
- Osmotic limit: with 300 mOsm inside and 150 mOsm outside, ideal water volume must become 2V0.
- Voltage clamp: the same initial state and the same external step, but Vm is held at the initial value.
- Pump placebo: pump flux is set to zero; any change in Vm must not be reported as an error of the full model.
- No data matching, no training split and no assumption that a successful run is a measurement finding.

## Frozen sensitivity parameters

For each parameter, exactly ±50 % is applied while all other values are locked:

- Lhyd, hydraulic water permeability; ±50 % changes only the hydraulic coefficient, not fw.
- gK, K+ conductance density.
- p0, maximum pump cycle per area.

For each run, V∞/V0, t50, χV and the full charge/amount balance are reported. Sensitivity is defined as (Y+50% − Y−50%)/Y0 for Y = V∞/V0, t50 and χV. No ranking or causal conclusion is made before running.

## Builds on

- The local working directory contained only BRIEF.md, NIGHT_PREAMBLE.md, QUESTION.md, ALLOW_WEB and agent.log; no previous model files, results.json or RESULTS.md existed to build on.
- question capability: Q031 in inputs/QUESTION.md, particularly the K01/K03/K06 requirements: simultaneous potential, ion amount, volume and balance data.
- Jakobsson 1980 above: n transport equations + electroneutrality + isotonicity.
- Farinas–Verkman 1996 above: geometric/osmotic reference and water permeability.

## Not redone

- No access to or use of ~/projects/bodytwin, restricted model data, the collaborative data or other internal result trees; no internal node-ids or files could therefore be reused.
- No external musculoskeletal solver run or comparison with external musculoskeletal solver: the task is a cell mechanism and lacks relevant external solver runtime or measurement material.
- No full channel model, pump regulation, active mechanics, pH, ATP coupling or heat transfer.
- No interpretation of sensitivity as measured data or as a universal bound across other cells.

## Run limit

One thread, only a short local simulation, numpy/scipy/pandas, no batch/LOSO sweep and no new files outside this directory.
