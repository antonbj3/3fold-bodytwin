BT-HX-Q115

## Results

Builds on `inputs/QUESTION.md`, `inputs/NIGHT_PREAMBLE.md` and `BRIEF.md`; no previous model files or result nodes were in the directory. `PREREG.md` was frozen before running and is included in `PREREG.sha256`.

The model is an age-structured crypt–villus model with stem cells, transit-amplifying divisions, migration, tip rejection, age-dependent cell area, uptake and barrier conductance. Reference: Kai, Y. (2021), *Biophysical Journal* 120(4):699–710, DOI `10.1016/j.bpj.2021.01.003`, Table 1, p. 704: 1 885 cells, 26 cells h⁻¹ and 65 h. Derived reference lifetime is 72,5 h.

## Prediction and criteria

- Baseline: 73,33 h (3,056 day), 1,15 % from reference; **PASS** according to ±30 %.
- Kryptproduktion: 22,78 cells h⁻¹ versus 26; **PASS** inom faktor 2.
- Cellantal: 1 670,6 versus 1 885; **PASS** inom faktor 2.
- At crypt-villus length factor 0,55: villus area 4,455 → 2,476 mm²; effective absorption area 3,614 → 1,763 mm²; uptake proxy 0,361 → 0,176 mm³/day. It is −60,43 % against the constant area zero.
- Normalized conductance: 0,01581 → 0,01049 mS/cm², only +4,85 % towards zero. The pre-registered mechanistic criterion (at least 10 % for both) is therefore **FAIL**.
- OIf uptake and barrier are measured simultaneously, empirical superiority is **UNKNOWN**; no measurement data has been found.

## Sensitivity and next steps

Biggest impact at ±50 %: TA-delningssannolikhet (max 152,16 %), stem cycle (100,00 %) and crypt migration power (62,59 %). `python3 -m unittest -v test_model.py` passes all tests; dimensional_analysis is PASS. The next resolution steps are paired 3D geometry, pulse-chase line tracking, uptake measurement and transepithelial conductance under the same perturbation. All numbers and provenance can be found in `results.json`.

---

## CORRECTION 2026-10-03: the PASS above are no longer valid

The text above is unchanged and not touched, as it is reviewed. But the numbers it rests on have fallen.

The migration rate of the cell was never a measurement. It fell out of `drag_pN_h_per_mm = 100.0` (`Q115_model.py:26`),
whose own parameter table on line 114 said **"chosen to give 0.03 mm/h; assumption"**. The
speed was the target and the drag was the dial: (2,1 + 0,9) / 100 = 0,03 mm/h = 30,0 µm/h.

With the sustained measurement 9,00 ± 0,465 µm/h (DOI 10.1096/fj.201601002) as declared input:

| | before | now |
|---|---|---|
| turnover time | 73,33 h | **223,16 h** |
| wrong versus 72,5 h-referensen | 1,15 % | **207,81 %** |
| cellantal | 1670 | 5103 |
| upptag | — | **+282,68 %** |
| `primary_turnover` | PASS | **FAIL** |
| `cell_census` | PASS | **FAIL** |
| `mechanistic_perturbation` | FAIL | **PASS** |

The 1,15 % compliance was produced **twice**: also `baseline_path_mm = 1.95` was chosen
"with v=0.03 to represent 65 h", and 65 h is the figure the reference's 72,5 h is derived
from. Thus, both the speed and the distance were set to meet the reference.

`--migration-speed-um-h 30.0` reproduces the old blocks identically, so nothing is lost.

Checked separately: the model age ceiling of 240 h does not explain the
outcome — the ceiling flux is exactly 0,0 at both 30,0 and 9,00 µm/h.

Status PENDING_INDEPENDENT_REVIEW.
