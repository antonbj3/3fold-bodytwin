# Supply external inputs

See [CITATIONS.md](../CITATIONS.md) for exact upstream-required citations and current primary-source licence/access checks (4 October 2026). Historical implementation notes remain source disclosures; this register records the current public terms separately.

The repository contains no scans, labels, anatomical meshes, patient-derived arrays, clinical records or model weights. Download from the source under its terms into a directory outside the repository. Retain the version, licence, archive hash, selected member hashes, units and coordinate frames. Extract only the members required by the chosen experiment.

| Input family | Source / acquisition | Licence scope |
|---|---|---|
| Bits2Bites | [dataset page](https://ditto.ing.unimore.it/bits2bites/) | CC BY-NC-SA 4.0, linked by current dataset page; registration required |
| Bite2Text | [dataset page](https://ditto.ing.unimore.it/bite2text/) | CC BY-NC-SA 4.0, linked by current dataset page; registration required; bind the selected archive version |
| ToothFairy2 | [dataset page](https://ditto.ing.unimore.it/toothfairy2/) | CC BY-SA 4.0 in official benchmark dataset.json; registration required |
| ToothFairy / IAN | [ToothFairy](https://ditto.ing.unimore.it/toothfairy/), [Maxillo](https://ditto.ing.unimore.it/maxillo/); bind the selected release/provider manifest | ToothFairy CC BY-SA per official challenge (version unspecified there); Maxillo data permission UNKNOWN; registration required |
| STS-Tooth3D | [Zenodo record 10597292](https://doi.org/10.5281/zenodo.10597292) | CC BY 4.0 in the source delivery |
| Mandibular defects | [Figshare 28052240.v2](https://doi.org/10.6084/m9.figshare.28052240.v2); selected member manifest remains an external prerequisite for X4 | CC BY 4.0 in the supplied release brief; no archive is bundled |
| Teeth3DS / Teeth3DS+ | [dataset site](https://crns-smartvision.github.io/teeth3ds/) and its linked access instructions | CC BY-NC-ND 4.0 for data in official challenge README; selected version binding required; website licence differs |
| Open-Full-Jaw | [upstream repository](https://github.com/diku-dk/Open-Full-Jaw), including its LFS download instructions | CC BY-NC-SA 4.0 in the checked source delivery |
| MMDental | [Figshare record 28505276](https://doi.org/10.6084/m9.figshare.28505276) | CC BY 4.0 in Figshare data metadata; article is separately CC BY-NC-ND 4.0 |
| Pulpy3D | [dataset page](https://ditto.ing.unimore.it/pulpy3d/) | Exact redistribution permission UNKNOWN in the source delivery |
| teethPreparationData | [authors' repository](https://github.com/intellident-ai/teethPreparationData) | MIT; Copyright (c) 2024 intellident-ai; retain licence notice with upstream assets |
| AlignerMovement | [Zenodo record 11280343](https://doi.org/10.5281/zenodo.11280343) | CC BY 4.0; source-labelled crown-motion data, excluded here |
| Published laboratory tables | Cunali DOI 10.1590/0103-6440201601531 and Sagheb DOI 10.1186/s40729-023-00473-3 | CC BY 4.0; only the specified numerical group summaries are bundled |

These source links are acquisition locators. This release checked public citation/licence metadata without downloading any dataset payload; archive-version binding remains the consumer's responsibility. A source title or shared case key does not establish patient matching across datasets.

## Configure paths

The [catalogue](CATALOG.md) lists the input families named by each source README. `demos.json` lists input-root variables found in the extracted code. The full original pipelines are SOURCE_ONLY_EXTERNAL_INPUTS; setting a variable supplies a path and does not supply missing fixtures, packages or empirical calibration.

| Variable | Meaning |
|---|---|
| DENTAL_CORPUS_ROOT | Literature corpus root, with the `europepmc/fulltext` subtree expected by applicable sources |
| DENTAL_DATA_ROOT | Anatomical dataset root; preserve the layout referenced by the chosen source |
| DENTAL_INPUT_ROOT | Other excluded workspace artefacts and input contracts |
| DENTAL_EXTERNAL_ROOT | Other excluded engines, tools, archives or input snapshots |
| DENTAL_WORK_ROOT | Writable scratch outside the repository; default is the system temporary directory |
| DENTAL_CASE_ID | A user-selected source case key; the shipped default is `synthetic` |
| DENTAL_PYTHON | Python interpreter for the entry scripts; default is `python3` |
| DENTAL_IMPLEMENTATIONS | Extracted implementation root; default is this repository's `implementations` directory |

Source-specific environment names already present in the original code remain visible there. For optional source scripts use `PYTHONPATH=.` from the repository root. Do not run the source acquisition, freeze or graph-registration helpers as a substitute for a ported pipeline; the default 23 profiles do not call them. Full pipeline ports must preserve their original criteria and source versions.

NumPy, SciPy, mpmath, h5py and Shapely suffice for the 23 profiles. Optional full implementations also reference packages such as matplotlib, lxml, trimesh, scikit-image, numba and torch, and external FE/CAM executables. Their interpreter/runtime snapshots and executable binaries are excluded; installing those packages alone does not establish reproduction.

## Four-observation laboratory input

X59 calibration and validation CSVs require `specimen_id`, `protocol_id`, `region`, `point_id`, `batch_id`, `point_pair_verified`, `ct_state_bias_calibrated`, `single_observation_bound_um`, `dry_ct_um`, `in_situ_replica_ct_um`, `sectioned_replica_ct_um`, and `sectioned_replica_optical_um`. Supply measured values and measured bounds. The future-design CSV contains specimen/protocol/region/batch identity and `replica_mean_um`; it must contain no CT or validation outcomes. The API refuses training specimens in validation predictions and reports conditional method-scale predictions rather than absolute true gap.

## Added public numerical fixtures

| Fixture | Source quantity and locator | Terms / scope |
|---|---|---|
| pulpotomy_publication.json | DOI 10.1111/iej.14144 Table 4 counts; DOI 10.4103/jcd.jcd_118_22 Chart 1 comparator allocation | Primary CC BY 4.0; comparator CC BY-NC-SA 4.0. Only numerical facts and locators are bundled, not article text or figures. The software licence does not relicense third-party expression. Comparator denominator ambiguity 32/33 is retained. |
| milling_publication.json | DOI 10.3390/healthcare9080983 Table 1 four regional milling RMS means | CC BY 4.0; group means, not a pure beam-deflection calibration |
| guide_publication.json | DOI 10.1111/clr.13578 published Table 4, MANDIBLE (accepted-manuscript Table 3); DOI 10.1186/s40729-020-00272-0 Comparison of accuracy/Figure 4; profile-specific locators retained | Numerical population group summaries only; primary article rights remain with publisher. Sample-to-population and implant-to-drill transfer uncalibrated. |
| drill_protocols.json | Three manufacturer depth/datum statements with document version and printed-page locator | Manufacturer document rights retained; no PDF/article redistribution, no commercial tool shape inferred |
| crown_diagnosis.json | DOI 10.1016/j.jds.2025.03.016 Table 1; DOI 10.4047/jap.2025.17.1.1 Table 1 and methods | Numerical facts/locators only. CAD-file RMS and manufactured/scanned cusp variation remain different quantities. |

The reviewed Teeth3DS licence locator is https://osf.io/download/9dutn/ (CC BY-NC-ND 4.0 in the source review). Dataset licence differs from article licence. Bite2Text's current public page links CC BY-NC-SA 4.0; the older local archive's exact version binding remains an external prerequisite. Maxillo and Pulpy3D public data terms remain UNKNOWN; establish them before redistribution. None of these datasets is bundled.
