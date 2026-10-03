BT-HX-Q043

## Result

The first runnable model is a synthetic, fixed-budget benchmark, not a measurement fit. It simulates one force channel and one differential surface-EMG channel at `1000 Hz` for `4 s`, with `120` motor units. The full unit-level model and the pooled one-factor model receive the same observed budget.

The frozen rule was `R_pooled >= 0.10` and `G_MU >= 0.20`, where `G_MU = 1 - R_MU/R_pooled`. The nominal run gives:

- pooled held-out joint NRMSE: `0.45674171680830905`
- full-model oracle NRMSE: `0.027003081745450925`
- error reduction: `0.9408788802254647`
- decision: `detail_needed = true`

The modeled force/EMG correlation is `0.7570554523923134`; the code-defined lag is `-80.0 ms` (negative means force leads EMG). The heterogeneity diagnostics are geometry-visibility CV `0.6028301028297379` and force-capacity CV `0.5407935377846519`. These are generated-model diagnostics, not anatomical measurements.

## Mechanistic answer

A pooled scalar drive is sufficient only in the matched limit: homogeneous units, proportional force/EMG temporal responses, and no observation noise. The homogeneous control has pooled NRMSE `0.07567907187800736` and returns `detail_needed = false`; the zero-drive control returns zero force and zero EMG.

Motor-unit detail is needed when the same spikes are observed through different weights and delays: size-ordered recruitment changes which units contribute to force, electrode geometry changes which units are visible and can cancel, and a finite action potential blurs spike timing. In this run, the nominal geometry and recruitment heterogeneity make the two channels non-rank-one even though they originate from one drive. Thus the useful detail is latent unit state plus geometry, not simply adding another observed channel.

The `20 ms` action-potential control increases pooled NRMSE to `0.5962449509566248` and lowers correlation to `0.5287597868329311`. This is directionally consistent with the cited source, but is not a direct reproduction of its correlation statistic.

## Reference and sensitivity

Kutch, Kuo, and Rymer, DOI `10.1152/jn.00956.2009`, Table 1 supplies the anchored `120`-unit, `5/20 ms` action-potential, `8 Hz` minimum-rate, `45 Hz` first-unit peak-rate, `10 Hz` peak-rate difference, `0.2` ISI-CV, and `90 ms` longest-contraction values. Figure 5B reports correlations `0.91`, `0.73`, `0.53`, and `0.30` at `5`, `10`, `15`, and `20 ms`; these values and locations are also stored in `results.json`.

All preregistered ±50% cases retain `detail_needed = true`:

| controlling parameter | −50% pooled NRMSE | +50% pooled NRMSE |
|---|---:|---:|
| geometry conduction length | 0.4660883355487717 | 0.4644831389193208 |
| action-potential duration | 0.45423300946665346 | 0.49363990855867645 |
| recruitment-threshold spread | 0.4814988576568068 | 0.5176120062699385 |

The aligned nominal pooled error is `0.45674171680830905`; phase-shuffled EMG gives `0.6594664708947369`, a worsening of `0.20272475408642782`. The force-kernel integral check is `1.0000000000000002`; all analytical tests pass.

## Limits and next resolution step

No measured data were fitted, and no claim is made that the synthetic force scale, electrode layout, geometry, or heterogeneity applies to a particular person or muscle. The full-model error is a forward-model oracle score, not empirical predictive validation. The gain-formula and null-control corrections made before accepting the run are recorded in `PREREG.md`; the corrected file is hashed in `PREREG.sha256`.

The next step is simultaneous, held-out human surface EMG and force recordings from the same muscle with calibrated electrode geometry, plus an independent motor-unit constraint such as spike-triggered averaging or validated decomposition. Re-estimate geometry, recruitment, and noise on training data, then test whether the frozen `R_pooled`/`G_MU` rule predicts improvement on unseen subjects. Until then, the defensible conclusion is conditional: motor-unit detail is required when electrode-weighted EMG and force cease to share a proportional latent response.

Run reproduction:

```text
python3 model.py
python3 -m unittest -v test_model.py
```
