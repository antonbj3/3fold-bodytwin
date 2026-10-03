BT-HX-Q058

## Resultat

I built a normalized first-principles model with fine detail (activation, filter, central adaptation, fast and slow fatigue) and a coarse model retaining only `a` and `m_s`. Calibration is sourced to Potvin & Fuglevand (2017), DOI 10.1371/journal.pcbi.1005581: `c0 = 0,050 s`, `beta_ct = 0,379` and `tau_h = 22 s`. No measurement data have been generated or fitted.

At `t_s = 79,5 s` (`t_p = 80,0 s`), naive switching failed: the force jump was `0,0381 F_max` and maximum error after perturbation `0,1240 F_max`, against frozen bounds `0,02` and `0,05` respectively. It lost `h`, `m_f` and `x` although `m_s` was transferred.

Closure switching retained `a`, `h`, `m_f` and `m_s`: the jump was `3,17e-6 F_max` and post-perturbation error `0,0491 F_max`; it thus met the criterion, but close to the bound and not universally. Placebo switching at zero history gave `0,1227 F_max`, showing that the coarse dynamics also need validation. Relative sensitivity for `tau_f`, `k_f`, `tau_s` was `0,00414`, `0,00426`, `-0,00314`; the error was dominated by omitted adaptation.

Conclusion: switching can preserve memory within this frozen tolerance window when all slow states are transferred and fast states are initialized or observed. Exact losslessness also requires the fast filter to be represented. The next step is paired measurements of the same excitation/preload history before and after switching, identification of `tau_f`/`tau_s`, and a test of state reconstruction. Unit checks and tests passed; the result is specific to this normalized surrogate, not a clinical threshold.

See `PREREG.md`, `model.py`, `test_model.py` and `results.json`.
