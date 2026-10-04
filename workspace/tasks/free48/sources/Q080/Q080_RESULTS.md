BT-HX-Q080

## Result

The frozen core model provides **support for the mechanistic hypothesis**, not a biological measurement result. For an identical total number of 880 receptors, four Gaussian patches give a 2.7723 times greater peak activation at the end of the pulse than uniform distribution. The CV difference increases by 1.03494. The total activation at the end of the pulse is 1.1723 times larger and the time-integrated activation 1.0167 times larger. Both frozen requirements (`peak_ratio >= 1.20`, `CV increase >= 0.10`) are fulfilled. Source: `results.json` → `primary`.

| Distribution | Peak at 0.10 s (receptor/µm²) | Total activation at 0.10 s (receptor) | CV at 0.10 s |
|---|---:|---:|---:|
| Uniform | 10.2680 | 58.4106 | 0.502716 |
| Patched | 28.4659 | 68.4739 | 1.537658 |

At 1.2 s a clear peak remains (101.0767 vs 29.6127 receptor/µm²), while total bound receptors are almost equal (462.0351 vs 460.1279). This is a deterministic consequence of local mass action, the ligand field, and the slow receptor/dissociation dynamics, not an observed effect.

## Model and reference

`model.py` uses the frozen first-principles equations
`dc/dt = kon*L*r - koff*c + DR*Laplacian(c)`,
`dr/dt = -kon*L*r + koff*c + DR*Laplacian(r)` and
`dL/dt = DL*Laplacian(L) - kclear*L + source(t)`.
The receptor density is `Rbar*(r+c)`, the activation proxy is the bound receptor density `A=Rbar*c`, and both patterns are renormalized to exactly the same integrated number. The unit check is consistent; largest relative receptor mass error over main, control, and sensitivity runs is `4.3924e-15`.

Looked up primary sources:

- Briddon et al., PNAS 101 (2004), Fig. 6b and text: `Kd = 33 nM` and approximately `55 receptor/µm²`; DOI `10.1073/pnas.0400420101`.
- Suzuki et al., Biophysical Journal 88 (2005), Fig. 4d and the conclusion: `210 nm` compartment scale and `45 ms` residence; DOI `10.1529/biophysj.104.048538`.
- Calebiro et al., PNAS 110 (2013), Fig. 2D/I: median `D = 0.052 µm²/s` for β1AR and `0.039 µm²/s` for β2AR; DOI `10.1073/pnas.1205798110`.

Patch width, pulse, domain, effective ligand diffusion and clearance are assumptions. No internal data or measurement data is used. Building on K01–K03 according to `inputs/QUESTION.md`; their underlying scripts/data are not included in the package.

## Sensitivity

| Parameter | Value | Factor | Peak Ratio | CV Increase | Total ratio at pulse end | Peak time, both patterns |
|---|---:|---:|---:|---:|---:|---:|
| `D_L` | 12.5 µm²/s | −50 % | 2.0583 | 0.740956 | 1.2668 | 0.10 s |
| `D_L` | 37.5 µm²/s | +50 % | 3.2698 | 1.16654 | 1.1242 | 0.10 s |
| `koff` | 0.33 s⁻¹ | −50 % | 2.7670 | 1.03344 | 1.1727 | 0.10 s |
| `koff` | 0.99 s⁻¹ | +50 % | 2.7775 | 1.03643 | 1.1719 | 0.10 s |
| `sigma_patch` | 0.15 µm | −50 % | 9.1382 | 2.64964 | 1.1564 | 0.10 s |
| `sigma_patch` | 0.45 µm | +50 % | 1.3839 | 0.468285 | 1.1834 | 0.10 s |

The direction and frozen passport requirement remain in all six cases. This is structural deterministic sensitivity, not a confidence interval. Source: `results.json` → `sensitivity`.

## Checks, limit and next step

The analytical zero diffusion limit deviates at most `1.0845e-11 receptor/µm²`; the zero pulse gives exactly zero activation; fields are non-negative; all validity checks pass. The PREREG hash is `6da63a31da2a47e4dbcc012d6e7bcab06c85a2c5810699ee6071082f3a70d6df` and matches `PREREG.sha256`.

The next resolution step is same-cell measured receptor coordinates and a calibrated ligand dose/dose field, followed by pre-registered fit of patch width, effective ligand diffusion and kinetics on holdout cells. Next, extracellular layer geometry and modeling of receptor–receptor cooperativity, desensitization and endocytosis are needed before quantitative biological interpretation.

Reproduce with `python3 model.py --output results.json` and `python3 -m unittest -v`.
