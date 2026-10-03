# PREREG — BT-HX-Q115

## Scope and hypothesis

This is a first-principles, one-thread model of one mouse small-intestinal villus unit coupled to eight crypts. The predicted primary quantity is the mean epithelial residence time from crypt exit to villus-tip shedding, in hours. The mechanistic hypothesis is that crypt production, migration distance, stochastic transit amplification, age-dependent cell area and age-dependent transport/barrier properties jointly create a moving absorptive area; a constant-area null cannot represent the same uptake and barrier predictions after a geometry perturbation.


## Frozen reference

Primary reference: Kai, Y. (2021), “Intestinal villus structure contributes to even shedding of epithelial cells”, *Biophysical Journal* 120(4), 699–710, DOI `10.1016/j.bpj.2021.01.003`. Table 1, p. 704, reports for the finger-like villus model: total epithelial cells `N = 1885 cells`, cell supply/shedding `n = 26 h^-1`, and shortest crypt-to-tip migration time `tau_f = 65 h`. The derived reference mean shedding age is

`T_ref = N/n = 1885/26 h = 72.5 h = 3.02 d`.

The same article’s Introduction (pp. 699–700) states epithelial renewal every 3–5 days. The table value and the derived value are kept separate: the latter is the numerical target, the former is the source check.

## Frozen approval criteria

1. **Primary turnover criterion:** the healthy model mean shedding age must satisfy `0.70*T_ref <= T_pred <= 1.30*T_ref`, i.e. 50.75–94.25 h.
2. **Production/census criterion:** the crypt output must be within a factor of two of 26 cells h^-1 and the healthy epithelial census must be within a factor of two of 1885 cells.
3. **Mechanistic criterion:** with the prespecified villus-path shortening factor 0.55, the age-resolved model must change effective uptake area and normalized barrier conductance relative to the constant-area null by at least 10%, while conserving the crypt-to-shed flux at steady state.
4. **Scientific answer to “better than constant area”:** UNKNOWN unless matched measured uptake and barrier data are supplied. The model comparison can establish non-identifiability of the mechanisms, not empirical superiority.

A numerical criterion is marked `PASS` only if all inequalities in that criterion hold. A failed primary criterion is reported as `FAIL`; it is not repaired by post-run tuning. Conservation, non-negativity and dimensional checks are mandatory; a violation is a model error.

## Countertests and falsifiers

- Constant-area null: hold healthy epithelial area fixed, remove age-dependent area and transport, and compare the same geometry perturbation.
- No-shedding limit: set the shedding hazard to zero and verify that the cell census grows by the crypt input rather than reaching a spurious steady state.
- No-migration limit: set migration velocity to zero and verify that the minimum residence time diverges rather than silently changing the hazard.
- Zero-production placebo: set crypt production to zero and verify zero steady-state epithelial cells.
- Conservation check: for every discrete update, `N(t+dt)-N(t) = J_crypt*dt - J_shed*dt` within numerical tolerance.
- The proposal is falsified for a dataset if a held-out pulse-chase/3D-histology dataset cannot reproduce both the area-time course and the age distribution under the same parameter set.

## Model decisions frozen before the first run

- One villus unit, eight crypts, four stem cells per crypt.
- Stem divisions every 16 h; symmetric expected renewal gives one transit-amplifying daughter per stem division.
- Each transit-amplifying lineage undergoes a six-round stochastic division process; each round divides with probability 0.5, giving expected terminal amplification `(1+0.5)^6 = 11.390625`.
- A cell cannot shed before the shortest crypt-to-tip path is traversed; after that, the constant tip extrusion hazard is 0.12 h^-1.
- Cell apical area and transport maturity increase from birth toward a mature plateau; barrier conductance density is higher for immature cells. These are explicit assumptions, not measured values.
- Baseline shortest path is 1.95 mm and net migration speed is 0.03 mm h^-1, giving a 65 h minimum transit. The geometry perturbation scales path length to 55% of baseline.

## Builds on

- `inputs/QUESTION.md` — task mechanism, mathematical starting identity and requested outputs.
- `inputs/NIGHT_PREAMBLE.md` — precision, preregistration, countertest and reporting constraints.
- `BRIEF.md` — bounded deliverable and time/resource limits.
- No existing model, result node or direct code anchor was present in this directory when the run began. Node IDs: none available in the supplied files.

## Not redone


## Reproduction

Run `python3 model.py --out results.json`, then `python3 test_model.py`. The first run is after this preregistration and its SHA-256 record are frozen.
