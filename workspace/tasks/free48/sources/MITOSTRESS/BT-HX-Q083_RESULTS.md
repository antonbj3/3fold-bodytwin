BT-HX-Q083

# Resultat

Status: FIRST RUNNABLE MECHANISTIC MODEL. The model is built on the previously interrupted work in `PREREG.md`, `model.py` and `test_model.py`; no internal BodyTwin data was used. `results.json` is regenerated from current code.

## Core mechanism

Each compartment is represented by `S_i = (V_i, n_i, psi_i)`. During fusion, mass-conserving exchange is used,

`J_s = K_s (c_i,s - c_j,s)`, `K_s = 1 / (tau_s (1/V_i + 1/V_j))`,

which gives each species its own time constant. At division, soluble species are pooled by the volume fraction `q`; `enzyme`, `damage`, `opa1` and `activation` use a lobe-first partition. The potential is divided by

`psi_1 = psi_mean - 2(1-q) delta_psi`, `psi_2 = psi_mean + 2q delta_psi`.

The mitochondrion's expected response is

`J_ATP = V k_ATP c_enzyme c_substrate/(K_M+c_substrate) exp((psi-psi_ref)/psi_scale) exp(-c_damage) - V k_leak`.

For vesicles, a slow activation/priming state is combined with a Ca2+-dependent release hazard. This is a reused topology model; synaptic exocytosis is not assumed to give a literal daughter organelle.

## Canonical mitochondrial run

| Storhet | Modellresultat | Publicerad referens |
|---|---:|---|
| `tau_post` | 77 s | Twig et al. 2008, Fig. 2B: 77 ± 71 s in INS1 |
| `gap_psi` | 12 mV | Twig et al. 2008, Fig. 3B: +5,9 ± 3,7 and −6,5 ± 3,8 mV |
| hyperpolarized/depolarized risk ratio | 4,1689 | Twig et al. 2008, Fig. 4: about 6 times within 3 min |
| response gap after 0,5 substrate perturbation | 0,0129093 nmol/s | model derivation, not measured data |
| max relative mass error | 1,85 × 10^-16 | frozen criterion ≤ 1 × 10^-10 |

The daughters have `psi = -168,5` and `-156,5 mV` respectively; the probability of the next event is 0,8688 and 0,2084 respectively. The volume-weighted mean potential is preserved exactly (`-162,5 mV`). The total-amount model makes the daughters equal (`gap = 0`; `total_only_gap_fraction_of_full = 0`). Placebo gap is 0, and the analytical complete-mixing limit has maximum concentration error 0.

All frozen criteria pass: mass conservation, time window, potential gap, risk ratio 4–8, placebo, total-amount null, analytical limit, unit check and finiteness. The unit check gives `K = 0,0416667 um^3 s^-1`, sample flow `0,0291667 nmol s^-1` and the dimensional check exactly 1.

## Vesikelgren

After a 30 s mixing/partition window, the two vesicles' activation concentrations are 0,88293 and 0,21707 nmol/um^3. With the frozen Ca2+ perturbation 0,5, expected load becomes 7,4045 × 10^-6 and 1,8227 × 10^-6 nmol respectively; the gap is 5,5819 × 10^-6 nmol. The total-amount model gives equal daughter load (`gap = 0`), while potential-null does not change the vesicle response. Reference values are Luo et al. 2015, Fig. 5C (5–6 % per open channel) and Kim et al. 2012, Fig. 3E–F (`3,7 ± 0,9 × 10^-3 s^-1` with Syt1; `1,16 ± 0,25 × 10^-4 s^-1` without Syt1).

## Sensitivity (±50 %)

| Parametergrupp | responsgap 0,5x / 1x / 1,5x (nmol/s) | riskkvot 0,5x / 1x / 1,5x | potentialgap 0,5x / 1x / 1,5x (mV) |
|---|---|---|---|
| `tau_enzyme`, `tau_damage`, `tau_opa1` gemensamt | 0,009998 / 0,012909 / 0,013967 | 3,024 / 4,169 / 4,743 | 12 / 12 / 12 |
| `delta_psi` | 0,018655 / 0,012909 / 0,007350 | 4,801 / 4,169 / 3,604 | 6 / 12 / 18 |
| `q` | 0,010230 / 0,012909 / 0,017529 | 1,381 / 4,169 / 3,396 | 12 / 12 / 12 |

No recalibration was done after the sensitivity run. `delta_psi` and `q` govern different parts of memory: potential gap and volume/lobe partition respectively. Full values and mass errors are in `results.json`.

## What does not follow

No frozen control failed. However, `J_ATP` is a normalized first-pass response, not a fitted metabolic network. No mtDNA copies, protein structure or cell-type-specific measured data are derived. The model's `tau_*`, initial amounts and response coefficients without a source are explicitly assumptions.

## Next resolution step

1. Measure species-specific exchange and partition efficiencies in a controlled fusion/fission.
2. Calibrate potential, OPA1 and damage simultaneously with time-lapse playback and perturbation.
3. Distinguish pre-fusion priming from post-fission memory in the vesicle case.

## Reproducerbarhet

Run `python3 model.py`, `python3 -m unittest -v` and `sha256sum -c PREREG.sha256`. Exact test results: 6 tests passed. Preregistration hash: `62ab12d1b2c882f7189db8cf9132d2c8939cab3299b941ed6e6dec7f390ab70e`.
