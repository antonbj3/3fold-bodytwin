BT-HX-Q127

## Results
Builds on `PREREG.md` (SHA-256: `0443a207051c2e1c10d1b2f90322530154aea9e9a64be5d377a7cb02249ac0a7`), `inputs/QUESTION.md` and a public longitudinal human study. The model in `model.py` has nuclear/fiber identity, satellite cell stock, domain-driven fusion, age/atrophy-dependent loss, protein-area balance, load signal, epigenetic priming and spatial repulsion; `test_model.py` passes four tests including zero flow, exponential loss and dimensional audit.

Reference: Cumming et al. (2024), *The Journal of Physiology* 602(17):4171–4193, DOI `10.1113/JP285675` (verified page: https://www.ovid.com/journals/jphy/fulltext/10.1113/jp285675~muscle-memory-in-humans-evidence-for-myonuclear-permanence). Type-II means: `M0=2.4`, `M1=3.3`, `M_D=3.2`, control `M_C,D=2.4` nuclei/fiber; frozen `R_ref=(3.2−2.4)/(3.3−2.4)=0.8888889`.

Deterministic prediction: `M1=2.9345089`, `M_D=3.0765813`, `M_C,D=2.3527240` nuclei/fiber; `R_excess=1.3542475`. Relative error is `52.35%`, so the ±30% criterion is **NOT MET**. The raw ratio `M_D/M1=1.0484144` is not the primary criterion. `R_excess>1` breaks the expected 0–1 retention bound; the bounded diagnostic value is therefore not used for approval. The 2 000-fiber stochastic run gave `R_excess=1.4027553`. The model's fusion flow was `0.0085089`, `0.0022092`, `0.0081554` events/fiber/day during first training, detraining and retraining; density at phase boundaries was `685.57`, `624.72`, `657.01` nuclei/mm².

`r_history=0` changed the final nuclear count by only `0.008909` nuclei/fiber; this is not motor memory. The model overpredicts excess retention through basal fusion during detraining and simultaneous control loss. Conclusion: the mechanism is runnable but quantitatively failed; biological retention and retraining effect are **UNKNOWN**.

Next step: longitudinal fiber-identified PCM1/DAPI measurement with registered biopsy position and injury site; measure fusion/loss, satellite cells, fCSA and fiber type per fiber. Without this data, nuclear retention should not be interpreted as functioning memory.

All numbers: `results.json`.
