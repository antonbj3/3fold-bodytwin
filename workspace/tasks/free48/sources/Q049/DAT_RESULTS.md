BT-DAT-Q049

## What was built on
`inputs/Q049_MODEL.py` + `inputs/Q049_results.json` (1D layered poroelastic column; computational core FAIL: `psi_s(1 s)=0.0532` against frozen interval [0.90, 1.00]) and `inputs/Q049_PREREG.md`. PREREG.md was frozen with sha256 `23ee6865…` before download. `agent.log` from the interrupted run contained only a directory listing — nothing to resume.

## Five verified public datasets (DATA_SOURCES.json)
| source | licence | size (B) | format | people | constrains |
|---|---|---|---|---|---|
| OAIZIB-CM, Zenodo `10.5281/zenodo.14934086` | CC BY-NC 4.0 | 12 652 466 616 | ZIP/NIfTI-1 | 103 test volumes (measured) | `superficial.thickness`, `calcified_barrier.thickness`, `normalization_strain` |
| 3D Knee MRI Cartilage Segmentation, `10.17632/py8tp69jjh.1` | CC BY 4.0 | 214 300 249 | ZIP (not submitted) | UNKNOWN | total column thickness 0.20+1.80 mm |
| Coronal-PD knee, `10.5281/zenodo.10673359` | CC BY 4.0 | 8 021 423 926 | ZIP/patient, PD-slices | 20 patients | lateral geom that the 1D column neglects |
| Knee adduction moment, `10.5281/zenodo.7215806` | CC BY 4.0 | 27 103 679 | .mat/.m | 1 person, 3 speeds ± brace | `REFERENCE_STRESS_PA`=36 500 Pa, `drainage_half_time_s`=0.528 s |
| NLM Visible Human | public domain | UNKNOWN | .raw/.rgb/.png @0.33/0.17 mm | 2 | mm units and coordinate frame for the thin layers |

Nine candidates rejected with the stated rule (MRNet lacks a licence number, OAI/OpenCap/fastMRI could not be verified here, Cartigram is a PDF, BodyParts3D/AddBiomechanics JS shells).

## Sample (samples/OAIZIB-CM_labelsTs.zip)
14 122 969 B ≤ 50 MB, sha256 `72a81672…` = the file's own. ZIP: 103 NIfTI-1 masks. File read: **160×384×384 uint8**, pixdim **0.7×0.364583×0.364583 mm** (unit code 1 = mm), sform code 1 with the affine `[0.698288,−0.025432,0,2.080258 / 0.048829,0.363695,0,−72.379135 / 0,0,0.364583,−74.478035]`, FOV 112×140×140 mm, five labels 1–5 + background. ROI names missing in the file and acquired metadata → label–anatomy = UNKNOWN.

## What failed
- **C1 FAIL**: 5 verified sources against frozen requirement 6. Reported, not padded.
- **The barrier lacks external data**: no public tidemark dataset — but the barrier effect is the largest in the model (0.978 relative in `J_bone(10 s)`).
- **The film lacks measurement entirely**: no published `film_conductance` (assumed 1e-8), `film_thickness`, `film_modulus` or `viscosity`. The film's "necessity" (Δ`psi_s`=0.0527) rests entirely on assumptions. This is the most important finding: the data gap is exactly where the model fails.
- Plausibility check on the sample: **not executable** (no ROI names; all bboxes span the entire FOV in two axes, so bbox ≠ thickness).

## Next steps
Acquire `info.zip` (51 680 B, same record) for the ROI names, measure a real thickness distribution against 0.20/1.80/0.10 mm, and search specifically for a deposited tidemark dataset. K01/K04/K05/K10 lack definitions in that bounded folder — UNKNOWN, not claimed.
