# BT-HX-Q005 — preregistration

## Scope and frozen target

This is a first, assumption-explicit, renal-transporter model. The primary organ
outcome is creatinine renal clearance (`CL_Cr`, mL/min), a transporter-sensitive
renal readout rather than a claim of glomerular injury. The secondary independent
readout is metformin renal clearance (`CL_Met`, L/h). Both are computed from
one common parameter chain, while independent residual errors are allowed to be
correlated. No patient, cohort, or internal measurement is used.

The value question is operationalized as: which uncertain molecular parameter
should be measured/computed more accurately if the only objective is to reduce
prediction error in `CL_Cr` per unit of calculation cost? The model reports the
same ranking for `CL_Met` and for a joint normalized error, but the primary
decision target is `CL_Cr`.

## Mechanistic hypothesis and falsification

A shared plasma unbound concentration enters OCT2-mediated basolateral uptake.
OCT2 flux saturates with the substrate affinity `K_OCT2`; the secretory flux
then passes through an MATE-mediated apical step. At steady state, the
creatinine and metformin readouts are `f_u * GFR` plus the resulting tubular
secretion clearance. Therefore, a better `K_OCT2` or `f_u` should reduce the
renal-outcome error budget, while a parameter that is downstream of a
saturating step may have little value. This is a mechanistic hypothesis, not a
claim that a particular drug causes kidney injury.

The null model is filtration only: `CL = f_u * GFR`, with no active secretion.
A parameter is counted as useful only if it produces a positive, numerically
resolved change in the primary outcome and beats the filtration-only null in
the frozen plausibility check. The model is rejected if any output is
non-finite, a flux is negative at nonnegative concentrations, the unit check
fails, or the primary prediction is outside 0.5–2.0 times the cited human
reference range. A parameter is not declared a winner if its score is within
10% of the next score; that case is reported as a tie.

## Public reference anchors

These anchors set scale and are not silently treated as validation data for the
new model.

1. Severance, Sandoval & Wright, *J Pharmacol Exp Ther* 2017,
   DOI `10.1124/jpet.117.242552`, Table 1 (`PMC5539587`): metformin OCT2
   `Jmax = 1046 +/- 171 pmol cm^-2 min^-1` and apparent
   `Ktapp = 518 +/- 45 microM`. These are in-vitro CHO-cell transport values.
2. Topletz-Erickson et al., *J Clin Pharmacol* 2021,
   DOI `10.1002/jcph.1750`, Table 3 (`PMC7984390`): iohexol-derived
   `GFR = 94.99 mL/min/1.73 m^2` in the metformin-alone arm. Figure 3A and
   Results: mean baseline `CrCl = 152.3 mL/min` and `117.2 mL/min` after
   transporter inhibition. Table 2: metformin renal clearance
   `30.0 L/h` and 24-hour urinary fraction `28.5%` in the control arm.
3. Mathialagan et al., *J Pharmacol Exp Ther* 2024,
   DOI `10.1124/jpet.123.001890`, abstract: active tubular secretion is
   approximately 30% of creatinine renal clearance and OCT2/MATE transport is
   implicated. This supports the mechanism but is not used to invent a
   parameter value.

If a cited value cannot be checked, mark it `OVERIFIERAD`; no unverified value
is used as a measured input. All other model parameters are explicitly marked
as assumptions or derived scenario values.

## Frozen model choices

- A 70-kg adult, total kidney, constant unbound substrate concentrations during
  the local clearances: `C_Cr = 100 microM` and `C_Met = 7.5 microM`. The
  metformin concentration is a scenario scale consistent with the cited
  control Cmax (1.3 microg/mL, using 129.9 g/mol); it is not a new observation.
- The nominal `f_u = 1.0` for both probes is a model assumption. Creatinine is
  treated as freely filtered; metformin is treated as effectively unbound in
  this bounded scenario. The uncertainty prior is lognormal with mean 0.95 and
  CV 0.20, then clipped below one to respect the physical bound.
- The OCT2 kinetic scale uses the cited metformin values. A fixed substrate
  factor `kappa = 0.40` maps the common OCT2 affinity to creatinine; this is a
  derived assumption, not a measured creatinine `K` value. A fixed factor
  `0.067` maps the cited metformin `Jmax` scale to creatinine. Both factors are
  frozen before the run.
- Effective renal tubular area is `2.0e5 cm^2` per kidney, a geometric
  assumption. Apical efflux is linear MATE permeability; its parameter is
  included for the tissue-exposure sensitivity, but a linear efflux step does
  not change steady-state secretion clearance when uptake is rate-limiting.
- Candidate molecular parameters and prior relative standard deviations (CV) are
  frozen: `K_OCT2_ref` 0.30, `Jmax_OCT2_ref` 0.25, `f_u` 0.20, and
  `P_MATE` 0.30. `K_OCT2_ref` and `Jmax_OCT2_ref` have lognormal correlation
  0.70, motivated by the substrate dependence reported in the cited OCT2 paper.
  A common `f_u` draw couples the two probes; independent residual readout
  errors have correlation 0.35, with residual CV 0.05 for `CL_Cr` and 0.08 for
  `CL_Met`.
- Costs are frozen relative units: refine `K_OCT2_ref` = 2, `Jmax_OCT2_ref` =
  1, `f_u` = 1, `P_MATE` = 2, and `GFR` = 3. These are planning units, not
  measured compute times. The primary score is
  `(RMSE_before - RMSE_after) / cost`, where `after` halves only the named
  parameter's log-scale CV. The error budget share is the one-at-a-time
  variance contribution divided by the total latent variance.

## Frozen acceptance criterion

A candidate is called the primary parameter worth improving only if all of the
following hold in the deterministic run and 20,000 seeded draws:

1. its primary-outcome score is positive;
2. halving its prior CV lowers primary `CL_Cr` RMSE by at least 10%;
3. its score exceeds every other candidate by at least 10% (otherwise report a
   tie); and
4. the ±50% one-at-a-time perturbation changes the primary outcome in the
   direction implied by the equation and remains finite.

A failure of any criterion is reported as `not established`; no post-hoc
parameter, cost, or threshold changes are allowed.

## Data provenance and boundaries

`Builds on`: the two input files in this job, the public primary studies listed
above, and the mechanistic OCT2/MATE/filtration equations. K02, K08, K10 and
K11 could not be located in the job directory or its accessible workspace; no
node content is invented.

claims, unmeasured organ injury, and a general literature review. The model
cannot identify a universal molecular parameter: its result is conditional on
the stated renal readouts, dose/concentration scale, priors, and cost units.

## Reproduction

Run `python3 model.py` to create `results.json`, then run
`python3 test_model.py`. The model is single-threaded, uses only NumPy, and
writes no data outside this job directory.
