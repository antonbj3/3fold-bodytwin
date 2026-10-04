#!/usr/bin/env python3
"""Score the twin's two decisions the way a weather service scores a forecast.

Where this came from. The improvement batch ran twelve angles on the free tier; the operator insisted on
lenses that start outside my own measurements rather than lines extending them, and the NEIGHBOUR lens --
"what does another field do routinely that we do not do at all?" -- returned this: forecast verification.
Weather services do not publish a point number. They publish a probability and then score it against a
reference forecast, so a user knows both the answer and how much to trust it. We publish point numbers.

Why it is worth building rather than noting. Three separate measurements tonight pointed at the same
thing from different directions: the richer corneal representation changed 0 of 20 power choices, the
instrument's own repeat range exceeds the clinical threshold in 43 of 300 eyes, and the flag work showed
unreliability is partly predictable in advance. More accuracy is not the binding constraint; the output
form is. A probability is a different output, not a better one, which is also what the USER lens asked
for independently.

The design, and the two circularities it has to avoid.

1. The forecast may not see the outcome it forecasts. So every probability here is leave-one-out: the
   error scale is fitted on the other patients only, and the patient's own error is never in its own fit.

2. The sharpness may not come from hindsight. A forecast that varies per patient is only honest if what
   makes it vary is known before surgery. For the implant power that is the sensitivity dR/dP, which
   comes from measured biometry: the same uncertainty in chosen power lands as more dioptres of
   refraction in an eye with a steeper conversion. For the toric decision it is the magnitude of corneal
   cylinder being corrected, also measured preoperatively. Neither is derived from the result.

The reference forecast is climatology -- the leave-one-out base rate, the standard against which skill is
defined in forecast verification. Beating it means the per-patient probability carries information the
cohort rate does not. Not beating it is a real outcome and is reported as one, because the decomposition
then says which half failed: calibration or discrimination.
"""
from __future__ import annotations

import json
import math
import statistics
from pathlib import Path

W = Path('')
OUT = W / 'results/ASSEMBLY_FORECAST_VERIFICATION'
TOL_D = 0.5          # the event: the decision lands within half a dioptre
BINS = 5


def p_within(sigma: float, tol: float = TOL_D) -> float:
    """P(|e| <= tol) for a zero-mean normal error of scale sigma."""
    if sigma <= 0:
        return 1.0
    return math.erf(tol / (sigma * math.sqrt(2.0)))


def scale_from(errors: list[float], covars: list[float]) -> float:
    """One parameter, fitted on OTHER patients: sigma_i = c * x_i.

    For a half-normal, E|e| = sigma*sqrt(2/pi), so c = mean(|e|/x) * sqrt(pi/2). One parameter is
    deliberate: twenty patients cannot support more, and a richer fit would buy apparent sharpness
    out of the same data that scores it.
    """
    ratios = [abs(e) / x for e, x in zip(errors, covars) if x > 1e-9]
    if not ratios:
        return 0.0
    return statistics.mean(ratios) * math.sqrt(math.pi / 2.0)


def verify(name: str, errors: list[float], covars: list[float], covar_name: str) -> dict:
    n = len(errors)
    obs = [1.0 if abs(e) <= TOL_D else 0.0 for e in errors]
    fc, clim = [], []
    for i in range(n):
        others_e = [errors[j] for j in range(n) if j != i]
        others_x = [covars[j] for j in range(n) if j != i]
        c = scale_from(others_e, others_x)
        fc.append(p_within(c * covars[i]))
        clim.append(statistics.mean(1.0 if abs(e) <= TOL_D else 0.0 for e in others_e))

    bs = statistics.mean((p - o) ** 2 for p, o in zip(fc, obs))
    bs_ref = statistics.mean((p - o) ** 2 for p, o in zip(clim, obs))
    skill = 1.0 - bs / bs_ref if bs_ref > 0 else None

    # Murphy decomposition: BS = reliability - resolution + uncertainty, over equal-width bins.
    base = statistics.mean(obs)
    unc = base * (1.0 - base)
    rel = res = 0.0
    diagram = []
    for b in range(BINS):
        lo, hi = b / BINS, (b + 1) / BINS
        idx = [i for i in range(n) if (lo <= fc[i] < hi or (b == BINS - 1 and fc[i] == 1.0))]
        if not idx:
            diagram.append({'bin': f'{lo:.1f}-{hi:.1f}', 'n': 0})
            continue
        pbar = statistics.mean(fc[i] for i in idx)
        obar = statistics.mean(obs[i] for i in idx)
        rel += len(idx) / n * (pbar - obar) ** 2
        res += len(idx) / n * (obar - base) ** 2
        diagram.append({'bin': f'{lo:.1f}-{hi:.1f}', 'n': len(idx),
                        'mean_forecast': round(pbar, 4), 'observed_frequency': round(obar, 4)})

    return {
        'decision': name,
        'n': n,
        'event': f'|error| <= {TOL_D} D',
        'observed_base_rate': round(base, 4),
        'sharpness_covariate': covar_name,
        'brier_score': round(bs, 5),
        'brier_score_climatology_reference': round(bs_ref, 5),
        'brier_skill_score_vs_climatology': round(skill, 4) if skill is not None else None,
        'murphy_decomposition': {'reliability_lower_is_better': round(rel, 5),
                                 'resolution_higher_is_better': round(res, 5),
                                 'uncertainty': round(unc, 5),
                                 'identity_check': round(rel - res + unc - bs, 9)},
        'reliability_diagram': diagram,
        'forecast_range': [round(min(fc), 4), round(max(fc), 4)],
        'leave_one_out': 'the error scale is fitted on the other patients only; no patient is in its own fit',
    }


def main() -> int:
    iol = json.loads((W / 'results/ASSEMBLY_IOL_DECISION/DECISION_V1.json').read_text())['rows']
    toric = json.loads((W / 'results/ASSEMBLY_TORIC_DECISION/TORIC_V1.json').read_text())['rows']

    # Implant power: the twin's miss expressed in dioptres of refraction at the spectacle plane, which
    # is what the tolerance is stated in. Covariate is the measured conversion |dR/dP|.
    iol_err = [r['twin_power_miss_D'] * abs(r['sensitivity_dR_dP']) for r in iol]
    iol_cov = [abs(r['sensitivity_dR_dP']) for r in iol]

    # Toric: residual cylinder, scored exactly as the toric cell scores it. Covariate is the measured
    # corneal cylinder magnitude -- the amount being corrected.
    tor_err = [abs(r['measured_posterior']['predicted_cyl_D'] - r['measured_cyl_D']) for r in toric]
    tor_cov = [abs(r['measured_posterior']['corneal_cyl_D']) for r in toric]

    results = [verify('implant_power', iol_err, iol_cov, 'measured refraction-per-power sensitivity |dR/dP|'),
               verify('toric_cylinder', tor_err, tor_cov, 'measured corneal cylinder magnitude, D')]

    summary = {
        'practice_adopted': 'forecast verification: probability output scored against a reference forecast',
        'field_it_comes_from': 'operational weather forecasting',
        'what_we_did_before': 'both decisions emitted a point number with no statement of confidence',
        'reference_forecast': 'climatology, i.e. the leave-one-out cohort base rate',
        'claim_type': 'capability',
        'control': 'the climatology forecast, which is the standard reference in this practice; a positive '
                   'skill score means the per-patient probability carries information the cohort rate does not',
        'scope': 'the same 20 and 69 patients the decisions were scored on; no new data',
        'results': results,
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'FORECAST_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    for r in results:
        print(f"  {r['decision']:15s} n={r['n']:3d} base={r['observed_base_rate']:.3f}  "
              f"BS={r['brier_score']:.4f}  ref={r['brier_score_climatology_reference']:.4f}  "
              f"skill={r['brier_skill_score_vs_climatology']}")
        m = r['murphy_decomposition']
        print(f"    reliability {m['reliability_lower_is_better']:.5f}  "
              f"resolution {m['resolution_higher_is_better']:.5f}  "
              f"identity {m['identity_check']:+.2e}  forecasts {r['forecast_range']}")
    print(f"  written: {OUT / 'FORECAST_V1.json'}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
