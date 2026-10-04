BT-HX-Q019

## Outcome

The bounded first runnable model is complete and its preregistered acceptance checks pass (`results.json:2-10`). It is a reduced, single-superficial-nephron model scaled to one kidney; it is not a clinical measurement or a fitted patient model.

## Mechanistic answer

The segment-specific variation that changes the independent organ response is not one universal coefficient: collecting-system water permeability controls the flow/osmolarity response, while thick-ascending-limb active NaCl reabsorption controls sodium excretion. Proximal water permeability is comparatively buffered by downstream transport. The dynamic equations are in `model.py:203`; the independent terminal observables are generated in `model.py:164` and the preregistered one-factor perturbations in `model.py:462`.

The baseline terminal outputs are urine flow `13.239098354243048 ml/min/kidney`, NaCl-equivalent `3400.469579125694 micromol/min/kidney`, Cl-equivalent `2843.772718163543 micromol/min/kidney`, urea `149.52092615391996 micromol/min/kidney`, and osmolality `522.9455200827852 mM` (`results.json:30-35`). These are derived model outputs. The checked primary reference is Layton & Layton (2019), DOI `10.1371/journal.pcbi.1006108`, whose Table 3 reports `0.62 ml/min` flow and `280 mM` urea; the discrepancy is reported, not hidden or converted into a measurement (`results.json:369-375`).

## Sensitivity

For the preregistered plus/minus 50% perturbations, normalized half-range sensitivity is:

- urine flow: `CD_Lp = -0.14977271705864037`, `TAL_vmax = -0.09746762631313273`, `PCT_Lp = -0.0015763033305374657`; ranking `CD_Lp > TAL_vmax > PCT_Lp` (`results.json:278-281`, `results.json:307-312`, `results.json:333-338`, `results.json:359-364`).
- NaCl-equivalent excretion: `TAL_vmax = -0.17877366595600433`, `PCT_Lp = -0.0014647147300674273`, `CD_Lp = -0.0003290450654906647`; ranking `TAL_vmax > PCT_Lp > CD_Lp` (`results.json:283-286` and the perturbation entries above).

## Controls and provenance

The filtration/unit check passes (`model.py:423`); the no-transport analytical limit returns `100.0 ml/min/kidney` against target `100.0` with relative error `1.4210854715202004e-16` (`results.json:133-137`). The null perturbation difference is `0.0`, and the integrated water-balance residual is `-1.906388819953679e-15` (`results.json:18-19`, `results.json:132-137`). Segment geometry and interstitial anchors are source values from Layton & Layton Tables 1–2; reduced residence times, water permeabilities, Michaelis–Menten `Vmax` values, and the lumped NaCl representation are explicit assumptions in `model.py:21-55` and `model.py:340-420`. No internal data or invented observations were used.

## Next resolution step

Calibrate the segment conductances against time-resolved independent urine volume, sodium, and osmolality measurements; then replace the well-mixed compartments with an axial countercurrent model and add nephron-population heterogeneity. This requirement is recorded in `results.json:142`.
