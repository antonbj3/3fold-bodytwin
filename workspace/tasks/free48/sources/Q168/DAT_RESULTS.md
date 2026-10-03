BT-DAT-Q168

## What this pass did

Located and downloaded **public measured** datasets that constrain the parameter set of the
synthetic tracer model in `inputs/Q168_model.py` (BT-HX-Q168), checked form, units and frame
for each, and recorded them in `DATA_SOURCES.json`. No model parameter was changed, the model
was not re-run, and no biological effect is claimed.

`PREREG.md` (sha256 `88fbc086…b5cf`, in `PREREG.sha256`) was frozen **before the first
download**; only endpoint liveness, HTTP status, licence/unit/size metadata and directory
listings were read beforehand. All numbers below come from `results.json`, which is produced
by `verify_samples.py` from the files in `samples/`.

## Datasets found (5 accepted, 3 rejected, 5 gaps)

| ID | Repository / accession | Licence | Downloaded | Subjects | Constrains |
|---|---|---|---|---|---|
| MW-ST000488 | Metabolomics Workbench `ST000488`/`AN000754` | CC BY 4.0 | 22 731 B | 36 human, plasma | endogenous pools, `sigma` frame |
| MW-ST005015 | Metabolomics Workbench `ST005015`/`AN008521` | CC BY 4.0 | 128 962 B | 148 samples, mouse+human tissue | host-conversion observability |
| MTBLS15673 | MetaboLights `MTBLS15673` | CC BY 4.0 (repo default; **not** verified inside the sample) | 27 269 675 B | 166 samples: 128 mouse, 32 reference, 6 *E. coli* | design precedent only |
| doi:10.7910/DVN/GK8J5W | Harvard Dataverse, PAMPA | CC0 1.0 | 153 920 B | compounds, not subjects | **nothing numeric** |
| CHEMBL4739856 | ChEMBL, bacterial TMA-lyase | CC BY-SA 3.0 | 66 332 B | recombinant assays | perturbation range only |

Total samples 27 650 002 B, largest file 25.94 MB, cap 50 MB/file — **C1 UPFYLLT**.
`DATA_SOURCES.json` is generated from `results.json`, so licence/size/format/subject counts
cannot drift apart.

## What the measured values actually say

**MW-ST000488, human plasma, micromolar, 36 subjects** (`results.json`
`datasets[MW-ST000488].stats_uM`): median TMAO **2.646 µM** (IQR 1.6195–3.8645), median TMA
**13.955 µM**, median carnitine **38.42 µM**, median betaine 38.27 µM, median choline 8.099 µM.
TMA exceeds TMAO in this cohort: median TMA/TMAO **5.2703**, median TMAO/TMA **0.19**. The
median TMA/carnitine ratio is **0.3677** with CV **0.035** over 36 subjects — a tight
stoichiometric invariant. Between-subject CV is 0.2166 (TMA) and 2.0211 (TMAO, driven by one
54.59 µM value; robust IQR-based SD ≈ 1.665 µM).

Two precision notes, kept separate on purpose:

- The BT-HX-Q168 anchor of 4.6 µM is a *general-population fasting* median; this cohort is
  sleep-apnoea/cardiovascular (`factors` field for all 36 samples). 2.646 µM is **not** a
  refutation and **not** validation of the anchor.
- The model's `sigma_plasma = 0.10 µM` is an *analytical* noise assumption. ST000488 gives
  *biological* between-subject spread. The two are not the same quantity, and the measured
  spread does not license replacing the noise model.

**MW-ST005015, 148 samples**: measured `d9-TMAO-1`, `d9-TMAO-2`, `TMA`, `TMAO` with
30 min / 120 min resolution across BAT, blood, eWAT, sWAT, gastrocnemius, kidney, liver,
pancreas, islet, hepatocyte, MIN6/INS-1E lysate, `human_liver` and `human_islet`; treatments
include AAV-FMO3, AAV-P153L, AAV-GFP, a 4-week antibiotic cocktail, d9-TMA 200 µM and
vehicle. Two hard defects, both reported not repaired:

1. The `Units` field lists **three** units in one column (`pmol/mg tissue; pmol/µg protein;
   µM`), so a cell cannot be mapped to a unit without the method section.
2. The **labelled substrate d9-TMA is absent** (non-empty cells 56/48/10/34 of 148). The
   preregistration anticipated a measured d9-TMA/d9-TMAO pair; it does not exist, so the
   ratio d9-TMAO-1/TMA has **n = 0** usable pairs and is reported as empty. This is logged as
   a preregistration/data mismatch in `results.json`, not silently fixed.

**doi:10.7910/DVN/GK8J5W (PAMPA)**: the two TSV files parse as `ID / SMILES / Y` with a
**binary class label only** (1740/295 and 50/92). There is no Papp, no cm/s and no exchange
area, so **no number for `P_epithelial` can be read from it**. Logged as `FRAME-GAP` under C3
— neither a pass nor a silent conversion.

**CHEMBL4739856**: 35 activities, 33 numeric, IC50 median **200 µM**, range **2.4–10 000 µM**
(conversion `1 µM = 1000 nM`, factor 0.001). This is inhibitor potency on a bacterial
enzyme, **not** a production flux, so it does not bound `Vmax_microbe` (criterion C6). The
legitimate use: every value is 1–3 orders of magnitude above the measured plasma carnitine
(38.42 µM median), so a microbial-production perturbation through this enzyme class is only
plausible as a **high-luminal, pre-systemic** intervention — which supports the lumen-site
design the model already implies.

**MTBLS15673**: 38 096 MAF rows, 166 samples, abundance columns are relative intensity with
**no calibration curve and no unit column**, and no TMAO/TMA analyte row is nameable in the
retrieved files (the compound appears only in the study description). Design precedent only.

## Counter-checks and negative controls

- The PAMPA dataset was the pre-registered route to `P_epithelial` and it failed the
  magnitude test; the permeability-versus-area non-identifiability from BT-HX-Q168 is
  therefore **not** resolved.
- ChEMBL IC50s were **not** converted into a `Vmax_microbe` bound.
- MAF relative intensities were **not** converted into µM.
- No taxonomy or gene-presence quantity was used as a flux bound anywhere.
- Plasma-matrix data was not used to bound lumen or portal parameters; the one transfer made
  (IC50 vs plasma carnitine) is named as an assumption.

## Coverage gaps (declared in `PREREG.md`, unchanged)

`Vmax_FMO3`, `K_FMO3`, `Vmax_microbe`, `K_microbe`, `P_epithelial`, `A_epithelial`,
`K_barrier`, `k_renal`, `V_s`, `V_l`, `V_p`, `sigma_lumen`, `sigma_portal`, `sigma_urine`
have **no public measured counterpart** in this pass. FMO3 kinetics: SABIO-RK returned
`302 → /ui/404` and `404`; ChEMBL has no human FMO3 record. Portal operator: no public
portal/liver-accessible TMA-TMAO dataset found. Rejected candidates are listed with the
observed reason in `results.json` (`rejected_candidates`), including the NHANES creatinine
route, which was dropped at the freeze rather than repaired after it.

## Status

Biological identification stays **UNKNOWN**. This pass added measured human plasma scales, a
measured labelled-TMAO product pool with tissue/time resolution, and a measured
perturbation-potency range for the microbial link. It did not make any causal estimand
identified, and it did not change a single model parameter.

## Next step

The single highest-value missing artefact is a public measured **d9-TMA → d9-TMAO** pair with
labelled substrate, stated units and a portal or liver-accessible matrix. MW-ST005015 already
carries the FMO3 and antibiotic arms for such a pair; the barrier magnitude needs an assay
that publishes Papp, not a class label.

## Files

`PREREG.md`, `PREREG.sha256`, `DATA_SOURCES.json`, `results.json`, `verify_samples.py`,
`make_data_sources.py`, `samples/` (22 files, sha256 manifest `samples/SHA256SUMS.txt`).
Reproduce with `python3 verify_samples.py && python3 make_data_sources.py`.
