# First principles — common research rule

Anton, 2026-09-26: prioritize first principles in every domain. Recursively break a question into its smallest relevant constituents, then reconstruct the observable from their interactions.

For each seed, deliver `DECOMPOSITION.json` containing:

1. **Decision and observable:** what changes if this mechanism matters; measurable output, units, time window and regime.
2. **Constituents and state:** molecules/receptors/cells/control volumes, or mathematical primitives and computational operators. Name locations, inventories, geometry and boundary conditions.
3. **Interactions:** stoichiometry, binding/unbinding, transport, signaling, force, energy or information flow. Distinguish a conserved quantity from a created/destroyed population or a statistical score.
4. **Governing relations:** equations, dimensions, physical constraints and timescales. For stochastic models specify events, propensities and observation process. For algorithms specify the operator, representation and its assumptions.
5. **Leaf status:** `DERIVED_UNDER_ASSUMPTIONS`, `EXTERNALLY_MEASURED`, `CONSTITUTIVE_CLOSURE`, or `UNKNOWN`. Include provenance and applicability. A named pathway, fitted Hill curve or literature citation alone is not a derivation.
6. **Reconstruction:** dimensioned ports, cancellation of internal fluxes, coupling direction and feedback/delays. Identify every closure required to go from microstates to the requested output.
7. **Discriminating test:** strongest matched existing baseline, adverse case, intervention/measurement that distinguishes alternatives, error budget and a result that would reject the hypothesis.

Decompose recursively while a finer level changes a prediction, validity condition, identifiability or error bound. Record the stopping argument for every leaf and its unresolved deeper mechanism. Molecular/electronic/quantum detail is a legitimate next layer when reaction energetics or binding physics controls the answer; an unevaluated atomistic model is not evidence that a tissue prediction became accurate. Retain a detailed specification even when execution needs a reduced model; quantify the reduction error where possible.

Inventories, positivity, charge, energy and thermodynamic consistency must be checked at the appropriate level. Activation signals are not substance flows. A causal biological network is not automatically a statistically independent evidence graph. Numerical/Lean certificates cover their explicitly stated mathematical model and arithmetic; empirical biological closures still need independent observations.

Start from existing code, corrections and failed results. Every new branch must identify its missing mechanism and downstream consumer. Return useful negative results or explicit knowledge gaps when the foundation is insufficient. No invented parameter values presented as measurements, no automatic scientific admission.
