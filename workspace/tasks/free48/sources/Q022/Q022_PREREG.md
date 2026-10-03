# BT-HX-Q022 — preregistration

## Scope and question

This is a bounded, first-principles lumped model of one cardiac excitation–contraction episode. It represents a normalized electrical activation signal, releasable sarcoplasmic-reticulum calcium, cytosolic calcium, coronary perfusion pressure/flow, oxygen supply, strain, and active/total force. It is a mechanistic hypothesis generator, not a patient model and not a fit to internal data.

**Predicted quantity.** The primary predicted quantity is the perfusion-to-force gain

\[
G_F(P)=F_{\mathrm{active,peak}}(P=70\ {\rm cmH_2O}) /
F_{\mathrm{active,peak}}(P=0\ {\rm cmH_2O}),
\]

together with the temporal ordering and lags electrical peak → activation peak → active-force peak and the strain minimum. The implementation also exposes normalized coronary flow and total wall force.

## Hypothesis

A perfusion increase can change active force without requiring a change in the electrical trace through two explicit routes: (1) pressure-driven flow increases oxygen/ATP availability and calcium reloading, and (2) a saturating perfusion-dependent contractility gain (the Gregg-like route) multiplies calcium-activated force. Electrical activation is upstream and must precede calcium release and force. At very low flow, the model includes an explicit hypoxia penalty; at normal autoregulated flow, the perfusion gain is bounded rather than extrapolated indefinitely.

## Reference values and source status

Values below were looked up in public primary literature before implementation. “Verified” means the cited figure or an exact result in the accessible primary abstract was checked. No value below is a measurement from the present job.

| Quantity | Primary source and location | Value | Unit | Status |
|---|---|---:|---|---|
| Maximum twitch cytosolic calcium anchor | Fabiato, A. (1981), *J Gen Physiol* 78:457–497, DOI `10.1085/jgp.78.5.457`, Fig. 10 and Results text | pCa 5.30–5.40 | pCa (dimensionless) | VERIFIED |
| Full myofilament activation anchor | Same primary source, Fig. 10 discussion | pCa approximately 4.90 | pCa (dimensionless) | VERIFIED |
| pCa 5.40 converted for the model | Derived from `10^-pCa` | 3.9810717e-6 | M | VERIFIED DERIVATION |
| pCa 4.90 converted for the model | Derived from `10^-pCa` | 1.2589254e-5 | M | VERIFIED DERIVATION |
| Perfusion-induced peak-force change | Schouten, V.J., Allaart, C.P., Westerhof, N. (1992), *J Physiol* 451:585–604, DOI `10.1113/jphysiol.1992.sp019180`, PubMed abstract Results paragraph 3 | 74 ± 20 (n=11) for 0 → 70 cmH2O | % peak force | OVERIFIERAD for figure/table; exact abstract result verified |
| Flow/contractility counterexample | Schulz, R., Guth, B.D., Heusch, G. (1991), *Circulation* 83:1390–1403, DOI `10.1161/01.cir.83.4.1390`, abstract Results | no significant change across 88–186 mmHg in the autoregulatory range; fall only at 57 ± 13 mmHg | mmHg, % wall thickening | VERIFIED abstract result |
| Whole-heart corroboration | Goto, Y., Slinker, B.K., LeWinter, M.M. (1991), *Circ Res* 68:482–492, DOI `10.1161/01.res.68.2.482`, abstract Results | flow +99 ± 76%, Emax +18 ± 15% at 93 ± 11 mmHg | %, mmHg | VERIFIED abstract result |

The Schouten paper’s numerical result is used as a calibration target, not treated as a universal law. The Schulz result is a preregistered falsification boundary: the model must not predict unbounded force gain in the autoregulatory range.

## Frozen model specification

The implementation will use the following planned parameter names and values; only a documented correction before the first run would require a new preregistration hash.

- Time step `dt = 0.001 s`; electrical control pulse is synthetic, dimensionless, Gaussian, peak time 0.160 s, width 0.025 s.
- Activation gate: `k_on = 100 s^-1`, `k_off = 50 s^-1`, with `a_dot = k_on e (1-a) - k_off a`.
- Calcium: `ca_rest = 1.0e-7 M`, `s_max = 1.0e-5 M`, `k_release = 20 s^-1`, `tau_ca = 0.12 s`, `tau_load = 0.25 s`; release and reloading are explicit mass balances.
- Myofilament activation: Hill coefficient 4, `ca_ec50 = 1.9952623e-6 M` (pCa 5.70 sensitivity assumption below the Fabiato 5.30–5.40 twitch anchor; not measured), clipped at the full-activation anchor.
- Perfusion: `p_ref = 70 cmH2O`, linear effective flow `q = clip((P-p_ext)/p_ref, 0, 1)`, `p_ext = 0 cmH2O`; `hypoxia_penalty = 0.20`; `gregg_gain = 0.392` derived so the static gain ratio at q=1/q=0 is 1.74; flow saturation is an explicit autoregulation assumption.
- Force: normalized length factor with Gaussian width 0.20, Hill force–velocity half-velocity 0.50 s^-1, normalized mechanical mass 1, passive stiffness 40 s^-2, damping 10 s^-1, and afterload 0.10 normalized force.
- Compression feedback: positive strain reduces effective flow by `compression_coef = 0.25` per normalized strain, capped to keep flow nonnegative.

The model must report both the static perfusion gain ratio and dynamic simulated peak-force ratio; they are not conflated.

## Frozen acceptance and error criteria

The run passes the preregistered model checks if all of the following hold:

1. `abs(static_gain_ratio - 1.74) <= 0.02`.
2. At the synthetic stimulation, `electrical_peak_time < activation_peak_time < active_force_peak_time`.
3. Electrical-to-activation peak lag is in `[0, 0.080] s`; activation-to-force lag is in `[0, 0.150] s`.
4. Increasing pressure from 0 to 70 cmH2O increases modeled peak active force and peak strain shortening; effective flow is finite, nonnegative, and nondecreasing with pressure at fixed strain.
5. Electrical-off placebo (`e(t)=0`) produces `max(active_force) <= 1e-8` normalized force; perfusion-only cannot create active force.
6. All states and outputs are finite, calcium and releasable calcium remain nonnegative, and the unit ledger passes.

A failed numerical or dimensional check is an error, not a scientific null. A failed sign/timing check is a failed model specification. Failure of the source anchor is not hidden by refitting to generated data.

## Counterfactuals and controls

- Electrical-off placebo: hold perfusion at 70 cmH2O while zeroing the electrical trace.
- Perfusion-off placebo: hold the electrical pulse while setting pressure to 0 cmH2O; the residual force is not interpreted as a measured ischemic value.
- Autoregulation control: saturate the flow/gain above `p_ref`; no additional force gain is permitted.
- Causal order control: an identical electrical trace is used for all pressure conditions; force cannot precede the causal activation chain.
- No internal data, no external solver runtime, and no clinical or joint-force inference are used.

## Builds on

- `BRIEF.md`; `inputs/NIGHT_PREAMBLE.md`; `inputs/QUESTION.md` (Q022; K01, K05, K06, K09 are named as inputs but are not present in this bounded job).
- Public primary anchors listed above: Fabiato 1981, Schouten 1992, Schulz 1991, and Goto 1991.
- The referenced `~/projects/bodytwin` anchor graph, `scripts/msk`, `docs/MECHANISM_*`, `bt_memory`, and prior `results/<id>` directories were not mounted/available in this job; no node ID is therefore invented.

## Not redone

