# PREREG — BT-HX-Q009

Status: frozen before the first model run. This is a first executable mechanistic model, not a validated clinical prediction.

## Hypothesis

A similar total plasma curve gives different free exposure in the target tissue when the free plasma fraction differs, or when transport, tissue binding or local clearance differs. Total plasma is therefore not a sufficient description of the free exposure trajectory.

## Predicted quantity and frozen criterion

Primary quantity:

`Q = AUC_t,u(f_u,p=0.40) / AUC_t,u(f_u,p=0.20)`

under an identical total plasma curve and identical other parameters.

Predicted model quantity: `Q = 2.00`. This follows from linearity in the free plasma input; it is an analytical structure check and not an empirical validation.

The technical PASS criterion is frozen at:

- the maximum relative deviation between the imposed total plasma curves is at most `1e-12`;
- the observed ratio lies between `1.90` and `2.10`;
- the same-binding control gives ratio `1.00` within `1e-10`;
- all unit checks and the analytical limiting case pass.

Physical PASS means only that the mechanism is implemented numerically and dimensionally. Physical VALID means that a matched tissue dataset reproduces the curve; this is `UNKNOWN` in this run.

## Frozen numerical experiment

- Time axis: `0–24 h`.
- Controlled total plasma curve for all cases: `C_p,total(t) = 10 exp(-0.08 t) + 1 exp(-0.8 t) mg/L`.
- Plasma cases: `f_u,p=0.40`, `f_u,p=0.20`, `f_u,p=0.05`.
- Tissue cases: same `k_transport=0.80 h^-1`, `k_elimination=0.05 h^-1`, `P_t=1.0 mg/L`, `K_d,t=4.0 mg/L`, `k_on,t=0.25 L/(mg h)`, `k_off,t=1.0 h^-1`.
- The plasma fraction's `K_d,p` is derived from `P_p=0.8 mg/L` and the chosen `f_u,p`; no claimed measured data are used.
- AUC is computed by the trapezoidal method on the same time grid. `C_t,u` is a microdialysis-like extracellular free concentration; `C_t,total=C_t,u+C_t,b` is an idealized tissue sum.
- Sensitivity: `f_u,p`, `k_transport` and `k_elimination` are varied by `−50 %` and `+50 %`.

## The model's core

Plasma:

`f_u,p = 1/(1 + P_p/K_d,p)`

`C_p,u = f_u,p C_p,total`

`C_p,b = (1 − f_u,p) C_p,total`

Tissue:

`P_t,free = P_t,total − C_t,b`

`dC_t,u/dt = (PS/V_t)(C_p,u − C_t,u) − k_elimination C_t,u − k_on,t P_t,free C_t,u + k_off,t C_t,b`

`dC_t,b/dt = k_on,t P_t,free C_t,u − k_off,t C_t,b`

`C_t,total = C_t,u + C_t,b`

`PS/V_t` is expressed as `k_transport` in the code. The model assumes fast plasma-binding equilibrium, passive membrane exchange and irreversible clearance of free tissue molecule. It makes no claim about intracellular concentration.

## Reference value

Verified primary source: Gill CM, Fratoni AJ, Shepard AK, Kuti JL, Nicolau DP. *Omadacycline pharmacokinetics and soft-tissue penetration in diabetic patients with wound infections and healthy volunteers using in vivo microdialysis.* Journal of Antimicrobial Chemotherapy 2022;77:1372–1378. DOI: `10.1093/jac/dkac055`, Table 2.

The source reports for infected patients and healthy volunteers, respectively:

- plasma fraction `0.21 (0.03)` and `0.20 (0.02)`, dimensionless;
- total plasma AUC `6.27` and `14.06 mg·h/L`;
- tissue AUC `0.82` and `1.37 mg·h/L`;
- tissue penetration `0.66` and `0.54`, dimensionless, defined as `AUC_tissue/(f_u,p AUC_plasma)`.

These values are an external reference and measurement-method anchor, not training data. The abstract's free plasma AUC for patients (`1.13 mg·h/L`) differs from Table 2 (`1.30 mg·h/L`); therefore that number is not used in the model criterion.

## Builds on

- `inputs/QUESTION.md`, K01/K11/K07 as the question's inputs.
- `source_repository/data/MECHANISM_ANCHOR_GRAPH.json`, node `MET-PHARMACOKINETICS-ADME`: one-compartment PK is seed design without executable plasma–tissue coupling.
- `scripts/msk/hepatic_clearance.py` read as background on the free fraction's significance; no existing code couples it to tissue.
- `inputs/NIGHT_PREAMBLE.md` and `BRIEF.md` for resource, provenance and writing rules.

## Not redone

No internal personal data, external solver runtime, patient-specific parameters, invented measurements or training on the literature curve. No clinical cutoff is made from tissue penetration. Time-dependent binding kinetics, pH partition, transporters, cell subdivision and intracellular measurement are not identified here.

## Countertests and error definitions

- Null model: the same `f_u,p` in both cases must not create an exposure difference.
- Transport countertest: `k_transport=0` must give zero tissue exposure from a zero initial state.
- Analytical countertest: without tissue binding or clearance, the free-tissue ODE must reproduce the solution for exponential plasma.
- `FAIL`: non-finite values, dimensional errors, failed limiting case or broken matching.
- `UNKNOWN`: empirical validity without matched free plasma and tissue data.
