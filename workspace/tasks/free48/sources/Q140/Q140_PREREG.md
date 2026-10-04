# BT-HX-Q140 — Preregistration

## Scope and node

This node is a first runnable mechanistic model of renal elimination of metformin. The model predicts the plasma amount/concentration trajectory after an intravenous 500 mg dose, total renal clearance in mL/min, and cumulative urinary mass in mg over 24 h. Filtration, segmental secretion, and segmental reabsorption remain separate mass flows. The model does not estimate oral absorption, hepatic metabolism, interindividual distribution, or a patient-specific clinical dose.

The bounded workspace contains no earlier implementation. `agent.log` records only the interrupted inspection of `BRIEF.md`; `model.py`, `test_model.py`, `PREREG.md`, and result files were absent before this run.

## Builds on

- Node `BT-HX-Q140` and the renal-flow question in `inputs/QUESTION.md`.
- `inputs/NIGHT_PREAMBLE.md` for the required node discipline, source separation, and output files.
- `agent.log` for the interrupted-session check; no prior numerical result was available to reuse.
- Public transporter and clearance sources are listed under the external reference below.

## Not redone

- No external solver runtime, fitting sweep, oral absorption model, or patient-specific inference was added.
- No measured urine or plasma time-course data were invented; unsupported reabsorption is kept UNKNOWN.

## Frozen hypothesis and endpoint

Primary endpoint: the effective renal clearance at the initial model concentration, defined as the total drug mass flux from the systemic plasma compartment into the kidney divided by total plasma concentration, in mL/min.

Secondary endpoints: plasma concentration at 5 h and 24 h, cumulative urinary mass at 24 h, and the filtration, secretion, and reabsorption contributions to renal clearance.

Prediction rule:

`CL_R = f_u,p GFR + CL_sec − CL_reabs`

The implementation uses a separate plasma amount, segmental lumen amounts and volumes, and a segmental epithelial-cell amount. Secretion is represented by basolateral OCT2 uptake and apical MATE1/MATE2-K efflux. Reabsorption is represented by an apical uptake branch and basolateral return. Water reabsorption changes luminal concentration and residence time. pH enters the ionization and apical H+-coupled transport factor.

## External reference and acceptance criterion

Reference: Graham GG, Punt J, Arora M, Day RO, Doogue MP, Duong JK, et al. (2011), “Clinical pharmacokinetics of metformin,” *Clinical Pharmacokinetics* 50:81–98, DOI `10.2165/11534750-000000000-00000`. The PubMed record reports a population mean renal clearance of `510 ± 130 mL/min` in healthy subjects and diabetic patients with good renal function and renal clearance `4.3 ± 1.5` times creatinine clearance. Source page: `https://pubmed.ncbi.nlm.nih.gov/21241070/`. The transporter assignment is cross-checked against Gong L, Goswami S, Giacomini KM, Altman RB, Klein TE (2012), *Pharmacogenetics and Genomics* 22:820–827, DOI `10.1097/FPC.0b013e3283559b22`, source page `https://pmc.ncbi.nlm.nih.gov/articles/PMC3651676/`.

Frozen acceptance criterion: the predicted primary endpoint must lie between `357` and `663 mL/min` (`510 × 0.70` to `510 × 1.30`). This is a calibration/validation target for a first mechanistic run, not a claim of patient-level accuracy. The filtration-only null model is rejected if its clearance is within the acceptance interval or if its 24 h urinary mass is within 30% of the full model. A negative mass, negative concentration, or failure of the analytic transport test is a model failure. The external 24 h urinary-mass value is UNKNOWN; the only urine-mass check in this run is conservation against the modeled dose.

## Data and assumptions

No internal data, patient data, fitted concentration-time series, or urine measurements are used. Published values are used only for the reference and transporter/ionization context. GFR, plasma flow, water reabsorption fractions, segment volumes, pH, transporter parameters, and effective distribution volume are transparent assumptions or derived quantities. The model is therefore exploratory and not a substitute for measured GFR, free concentration, transporter abundance, urine volume, or urine drug concentration.

The reabsorption capacity is set to a small effective baseline because no in-vivo metformin reabsorption component is identified by the cited transporter review. It is not interpreted as a measured zero. A ±50% perturbation of that capacity is reported as an identifiability/sensitivity result.

## Counterfactual and failure tests

1. `CL_filt = f_u,p GFR` is the filtration-only null model.
2. Set all segmental secretion capacities to zero while retaining the filtration term.
3. Set reabsorption capacity to zero and to twice its baseline value.
4. Reduce each of the three prespecified dominant candidates (`GFR`, epithelial secretory capacity, and water-reabsorption fraction) by 50% and increase each by 50%.
5. The analytic no-transport limit must return the filtered load without adding a secretory or reabsorptive mass flow.

The claim that filtration, secretion, and reabsorption are all required is not accepted for metformin from clearance alone: filtration and secretion must be present to explain a clearance above GFR, while a reabsorption term is needed only if the observed filtered load exceeds the observed urinary excretion. With the present reference, reabsorption remains UNKNOWN rather than being reported as observed.

## Frozen outputs

- `PREREG.sha256`: SHA-256 of this file before the first model execution.
- `model.py`: equations, parameter table, dimensional checks, deterministic run entry point.
- `test_model.py`: analytic limit and conservation tests.
- `results.json`: every numeric endpoint and sensitivity result.
- `RESULTS.md`: first line exactly `BT-HX-Q140`, endpoint comparison, acceptance status, limitations, and next data/geometry/measurement requirements.
