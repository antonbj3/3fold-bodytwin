BT-HX-Q013

## Results

A workable first-principles model was built in `model.py` with separate PK and cell branches:

1. Oral first-order absorption and one-compartment clearance: `dA_g/dt=-ka A_g`, `dA_p/dt=F ka A_g-CL A_p/V`.
2. Lokal exponering: `dC_local/dt=k_in C_p-k_out C_local`.
3. Target occupancy: `O=C_local^h/(EC50^h+C_local^h)`.
4. Power: `E=E0+gamma Emax O` .

`gamma=1` is unchanged cell sensitivity. `delta_E_PD=E_kombination-E_PK(C_local)` tests whether the effect deviates after local exposure is controlled. An independent local measurement can be specified to `decompose`; plasma alone is not sufficient for cell sensitivity grading.

## Source and derivation

**Source:** Tachibana M et al., *Clinical and Translational Science* 2025;18(11):e70330, DOI `10.1111/cts.70330`, PMID `41150706`, Table 3. Digoxin 0,25 mg + valemetostat 200 mg QD versus digoxin ensamt gav `AUC_last` GMR **1,27** (90% CI **1,06–1,52**, `n=16`). The source is systemic plasma exposure, not local cellular concentration.

**HDerivation:** The F branch was frozen to `F_mult=1,27`, `CL_mult=1`; the model then gives `R_plasma=R_local=1,27`. This is an identification assumption, not an observed partition between absorption and clearance.

**Hypothesis/antagande:** `gamma=1,5` is a synthetic sensitivity scenario. Other PK, EC50, Emax and transfer values ​​are normalized assumptions. No cell or patient measurement data has been created.

## Driving

`python3 model.py --output results.json` succeeded. Calibration passed: plasma AUC-ratio **1,2700000000000002**, local AUC-ratio **1,27**. All five tests in `test_model.py` passed, including analytical limit case `C_local=EC50 -> O=0,5`.

| scenario | AUC plasma | AUC local | max power | residual/Emax | model class |
|---|---:|---:|---:|---:|---|
| medicine only | 1,999991 | 2,857142 | 0,318790 | — | — |
| PK-only | 2,539989 | 3,628571 | 0,372778 | 0,000000 | `PK_ONLY` |
| cell-only | 1,999991 | 2,857142 | 0,478186 | 0,159395 | `CELL_SENSITIVITY` |
| PK + cell | 2,539989 | 3,628571 | 0,559166 | 0,186389 | `CELL_SENSITIVITY` |

Sensitivity, ±50 %:

- `clearance` : AUC local **7,254680 / 3,628571 / 2,419048**; max power **0,653313 / 0,559166 / 0,492354**.
- `EC50`: AUC local unchanged **3,628571**; max power **0,814650 / 0,559166 / 0,425671**.
- `gamma`: AUC local unchanged **3,628571**; max power **0,279583 / 0,559166 / 0,838749**.

## What remains

The next resolution step requires a factorial measurement with drug-only, PK-only, cell-only and combination. At the same time, plasma exposure, independent local target-site exposure and the same cellular response are required. Then `F_mult` and `CL_mult` can be identified separately and `EC50`/`gamma` estimated with a frozen model selection. Until this exists, the `F` verse is the `CL` breakdown, the cell type, and the true sensitivity factor is `UNKNOWN`.
