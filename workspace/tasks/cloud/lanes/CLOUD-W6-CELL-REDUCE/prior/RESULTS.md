# CLOUD-W3-CELL-COUPLE: Results

**Verdict:** the numerical criterion is **PASS** (synthetic check only). Empirical calibration is **UNKNOWN** because no matched public data could be reached. One planned counterexample did not finish and is reported as **not run** (see §6).

This result comes only from a model. It does not show that real cells behave this way.

## 1. What was built
The model is an original toy cell, written in stdlib-only Python 3.11 with no third-party packages. It is in `code/cellcouple.py`.

- **Geometry:** a spherically symmetric cell with nondimensional radius R = 1.
  - An impermeable nucleus spans r < 0.30.
  - A mitochondrial shell at r ∈ [0.35, 0.55] makes ATP from ADP (Michaelis–Menten in ADP).
  - A sub-membrane shell at r ∈ [0.85, 1.0] consumes ATP back to ADP (first order in ATP).
  - The plasma membrane has no adenylate flux.
- **Equations:** all quantities are nondimensional. Length is scaled by R, concentration by the initial total adenylate A0, and time by R²/D_ATP.
  - ∂T/∂t = ∇²T + Da_p·χ_m·P/(κ+P) − Da_c·χ_c·T
  - ∂P/∂t = δ∇²P − Da_p·χ_m·P/(κ+P) + Da_c·χ_c·T
- **Dimensionless groups:**
  - Da_c = k_c R²/D_T
  - Da_p = V_p R²/(D_T A0)
  - κ = K_m/A0
  - δ = D_ADP/D_ATP
- **Baseline values (assumed, not calibrated):** Da_c = 5, Da_p = 50, κ = 0.1, δ = 1.2. Initial state: T = 0.9, P = 0.1.
- **Numerics:**
  - Finite volumes on spherical shells, with exact shell volumes and face areas r_f².
  - Each region is represented by its exact volume fraction in each cell.
  - Backward Euler in time.
  - Picard iteration on the Michaelis–Menten term, with a lagged denominator, run to 1e-12. This makes the system matrix an M-matrix, which keeps concentrations non-negative.
  - A 2×2 block-tridiagonal direct solve.
  - Steady states use pseudo-transient continuation: dt starts at 0.01, doubles each step up to 100, and stops when the change is below 1e-9.

**Metrics reported:**
- `gradient_index` = 1 − ⟨T⟩_consumption shell / ⟨T⟩_mitochondrial shell
- `effectiveness` = ATP flux ÷ the flux of the same cell if it were well mixed
- `ATP_ADP_cons`, the ATP/ADP ratio in the consumption shell

## 2. Commands (exact)
```
cd code
python3 run_all.py 240     # budgeted run, writes ../results.json; took 25.6 s
```
- Python 3.11.15. No random elements except the prior sampling, which uses seed 20260924.

## 3. Preregistered criteria
| Criterion (frozen, from PREREG.md) | Result | Status |
|---|---|---|
| Coupled mass/energy balance < 1e-6 | Worst relative error 1.36e-11 across 8 transient runs: baseline and stress parameters, N = 112 and 448, dt = 0.001, 0.01 and 0.1. This covers both adenylate drift \|ΣT+ΣP − A0\|/A0 and the ATP budget residual \|ΔΣT − ∫(production − consumption)dt\|/A0. | **PASS** |
| Mesh convergence < 2 % | Relative change between N = 224 and N = 448 was at most 5.8e-5 at baseline and 2.9e-4 at stress (Da_c = 100, Da_p = 1000). The observed order was 1.99–2.00. | **PASS** |
| Empirical calibration | No matched data were available (see §5). | **UNKNOWN** |

## 4. Baseline and nondimensional sensitivity
Steady state at N = 112 with baseline parameters:
- T_mito = 1.0125
- T_cons = 0.6540
- ATP/ADP in the consumption shell = 1.996
- gradient_index = 0.354
- effectiveness = 0.683 (the well-mixed T would be 0.957)
- Steady production/consumption mismatch < 1e-11

Numerical uncertainty: the Richardson-type discretisation error is below 1e-4 relative for every metric.

**Elasticities** d ln(metric)/d ln(parameter), central ±1 %, N = 112:

| parameter | gradient_index | effectiveness | ATP/ADP (cons.) |
|---|---|---|---|
| Da_c | 0.657 | −0.298 | −1.035 |
| Da_p | −0.042 | −0.005 | 0.166 |
| κ | 0.025 | 0.005 | −0.111 |
| δ | 0.029 | 0.195 | 0.876 |
| m_out (mito shell outer radius) | −1.511 | 0.605 | 2.538 |
| c_in (consumption shell inner radius) | −2.644 | 1.159 | 4.257 |
| r_nuc | 0.000 | 0.004 | 0.000 |

In this toy, **organelle geometry matters most**, specifically the distance from the mitochondria to the consumption sites. The diffusion Damköhler number Da_c comes next. Da_p and κ matter little because production is near saturation at baseline.

**Da_c sweep** at a fixed ratio Da_p/Da_c = 10, so the well-mixed state stays at T = 0.957:

| Da_c | 0.01 | 0.1 | 1 | 3 | 10 | 30 | 100 |
|---|---|---|---|---|---|---|---|
| gradient_index | 0.001 | 0.012 | 0.104 | 0.252 | 0.514 | 0.750 | 0.905 |
| effectiveness | 0.999 | 0.990 | 0.912 | 0.779 | 0.527 | 0.282 | 0.110 |

The well-mixed assumption is within about 1 % only when Da_c ≲ 0.1. It is off by more than 10 % once Da_c ≳ 1.

**Prior propagation (synthetic):** the four kinetic groups were each drawn log-uniformly within ×/÷3 of baseline, 200 draws at N = 56.
- gradient_index: 5/50/95 % = 0.175 / 0.369 / 0.595
- effectiveness: 5/50/95 % = 0.404 / 0.707 / 0.881
- Minimum concentration over all draws: 2.1e-7, so all draws stayed non-negative.

This spread reflects the assumed prior. It is not an empirical uncertainty.

## 5. Data and sources
- I searched for a public structure/function reference matched to one species or cell line, for example organelle geometry together with ATP diffusivity and turnover. I could not reach any source. Every host I tried was blocked by the session egress proxy: bionumbers.hms.harvard.edu, api.crossref.org, doi.org, europepmc.org, ncbi.nlm.nih.gov, pmc.ncbi.nlm.nih.gov, ebi.ac.uk, openorganelle.janelia.org, allencell.org, zenodo.org, wikipedia.org, arxiv.org, biorxiv.org and api.openalex.org.
- **No URLs or DOIs were read, so none are cited.**
- No parameter value here comes from the literature, and none should be read as a physiological estimate.
- Mixing species or cell lines would be a further risk: for example, geometry from one cell line with diffusivities or kinetics from another. That was not attempted, because no data were obtained.

## 6. Counterexamples and failure modes
1. **Non-conservative discretisation is detected.** I deliberately used cell-centre areas in place of face areas. This gave an adenylate drift of 8.6e-3 and an ATP budget residual of 3.25, so the balance check can tell a correct scheme from a broken one.
2. **Conservation and mesh checks do not guarantee a physical answer.** With Newton linearisation and a single huge first pseudo-time step (dt = 100), one parameter draw converged to a spurious steady state:
   - Draw: Da_c = 2.46, Da_p = 136, κ = 0.051, δ = 1.076.
   - Spurious state: min concentration −22.4, ATP/ADP = −1.05, effectiveness = 17.8.
   - It still balanced production and consumption to round-off. In an exploratory run whose output was not saved, it also looked grid-independent (N = 56, 224, 448).
   - The shipped Picard solver with a ramped time step gives a physical state for the same draw: effectiveness = 0.811 and min concentration ≥ 0.
   - A positivity check is therefore needed in addition to the preregistered criteria.
3. **The well-mixed approximation fails at Da_c ≳ 1.** At Da_c = 100 the effectiveness is 0.11.
4. **Stiff consumption layers need resolution.** At Da_c = 100 the coarsest grid (N = 28) is off by 1.97 % in effectiveness (0.1079 vs 0.1100). It passes the 2 % criterion only from N = 56 upward.
5. **Not run / incomplete:** a planned under-resolution counterexample at Da_c = 1e4 (with grids up to N = 896) did not finish in about 13 minutes with the Picard solver. I stopped it and removed it from the budgeted run, so its result is **UNKNOWN**. The 2 % mesh result has therefore been shown only for Da_c ≤ 100.

## 7. Limitations
- The model is 1-D spherically symmetric. Real mitochondrial networks are 3-D and heterogeneous.
- It tracks only two species. There is no phosphocreatine shuttle, AMP, Pi, pH or buffering.
- Diffusion is free, with no crowding or anisotropy, and kinetics are simple Michaelis–Menten and first order.
- All values are nondimensional and uncalibrated.
- The PASS applies to numerical consistency only, not to biological validity.
