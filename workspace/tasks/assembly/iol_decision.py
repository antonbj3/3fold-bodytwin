#!/usr/bin/env python3
"""Use the assembled eye chain to make the surgical decision, and score it against what happened.

Why this exists. The eye chain already runs end to end: our molecular hydration model gives the stromal
refractive index, the patient's measured biometry gives the geometry, a modelled intraocular lens gives
the implant, and the result is a predicted spectacle refraction. R15 scored that PREDICTION against 20
held-out patients and landed inside the clinical 0.25 D tolerance in 7 of them, with a bias of -0.3865 D
under a thin lens and +0.1455 D under a thick one.

But predicting a refraction after the fact is not what a surgeon needs. The decision in cataract surgery
is WHICH LENS POWER TO IMPLANT. That is an inverse question, and the chain can answer it: the predicted
refraction is monotone in implant power, so the power that would have produced a target refraction can be
found by bisection. This script does that and then asks the only question that matters — would the
twin's choice have left the patient closer to the target than the choice that was actually made?

Honesty of the test. A model with a known bias would recommend a biased power, so the bias has to be
removed. It is removed by LEAVE-ONE-OUT: for each patient the bias is estimated from the other 19 and
never from that patient's own outcome. The patient's own measured refraction enters only in the final
scoring. The raw, uncorrected recommendation is reported alongside, so the contribution of the
correction is visible rather than hidden.

What this is not. The chain carries a mean-curvature scalar cornea with the preoperative geometry frozen
after surgery, no astigmatic vector and no patient-specific wavelength, and the modelled lens is a
published design rather than the implanted serial-number lens. Those limits are inherited from R15 and
are stated in its own scope line. So this is a decision-grade EXERCISE on a real cohort, not a clinical
recommendation.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

R15 = Path('./results/LANE_EYE_OPTICAL_TWIN/r15')
OUT = Path('./results/ASSEMBLY_IOL_DECISION')
sys.path.insert(0, str(R15))
sys.path.insert(0, str(R15.parent))

from clinical_anchor import predict          # the assembled chain, unchanged
TARGET_D = 0.0                               # emmetropia at infinity, the usual surgical target
TOL_D = 0.25                                 # clinical autorefractor tolerance
POWER_STEP = 0.5                             # lenses are manufactured in half-dioptre steps


def refraction(row: dict, power: float, model: str = 'thick') -> float:
    """Predicted spectacle refraction at infinity for a hypothetical implant power."""
    r = dict(row)
    r['IOLP'] = power
    return predict(r, model)[0]['refraction_infinity_D']


def invert(row: dict, target: float, model: str = 'thick') -> float:
    """Implant power whose predicted refraction equals the target. Monotone, so bisect.

    The bracket is checked rather than assumed: a silent non-bracketing bisection would return an
    endpoint and look like an answer.
    """
    lo, hi = 5.0, 40.0
    f_lo, f_hi = refraction(row, lo, model) - target, refraction(row, hi, model) - target
    if f_lo * f_hi > 0:
        raise ValueError(f'target {target} not bracketed on [{lo}, {hi}]: {f_lo}, {f_hi}')
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f_mid = refraction(row, mid, model) - target
        if f_lo * f_mid <= 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)


def main() -> int:
    geometry = json.loads((R15 / 'GEOMETRY_INPUTS_R15_V1.json').read_text())
    held = json.loads((R15 / 'HELD_REFRACTION_R15_V1.json').read_text())
    outcome = {r['sheet_row']: r for r in held['rows']}
    rows = [r for r in geometry['rows'] if r['IOLT'] == 0]

    # Measured spherical equivalent, the quantity the chain predicts.
    mrse = {r['sheet_row']: outcome[r['sheet_row']]['sph'] + outcome[r['sheet_row']]['cyl'] / 2
            for r in rows}
    # Prediction error at the power actually implanted, per patient.
    err = {r['sheet_row']: refraction(r, r['IOLP']) - mrse[r['sheet_row']] for r in rows}

    results = []
    for r in rows:
        key = r['sheet_row']
        others = [err[k] for k in err if k != key]
        bias = sum(others) / len(others)          # leave-one-out: never this patient's own outcome

        raw_power = invert(r, TARGET_D)
        corrected_power = invert(r, TARGET_D + bias)
        chosen = float(r['IOLP'])

        # A counterfactual outcome cannot be scored against measurement: we do not know what this
        # patient would have seen with another lens. Scoring our recommendation against our own
        # prediction was circular and gave an impossible 20 of 20 -- that version is withdrawn.
        #
        # What CAN be tested without circularity: the HINDSIGHT-CORRECT power. The measured refraction
        # at the implanted power says how far off the target the eye actually landed, and the local
        # sensitivity dR/dP converts that miss into the power that would have hit the target. Then the
        # question is whether our recommendation sits closer to that hindsight power than the
        # implanted one did. The measured outcome is the anchor; only the sensitivity comes from the
        # model, and a sensitivity is a far milder quantity than an absolute prediction.
        dRdP = (refraction(r, chosen + 0.5) - refraction(r, chosen - 0.5)) / 1.0

        rec = round(corrected_power / POWER_STEP) * POWER_STEP
        hindsight = chosen - (mrse[key] - TARGET_D) / dRdP
        results.append(dict(
            sheet_row=key,
            implanted_power_D=chosen,
            recommended_power_raw_D=round(raw_power, 4),
            recommended_power_bias_corrected_D=round(corrected_power, 4),
            recommended_power_on_manufacturing_grid_D=rec,
            power_difference_D=round(rec - chosen, 4),
            leave_one_out_bias_D=round(bias, 4),
            measured_refraction_D=mrse[key],
            sensitivity_dR_dP=round(dRdP, 4),
            hindsight_correct_power_D=round(hindsight, 4),
            twin_power_miss_D=round(abs(rec - hindsight), 4),
            implanted_power_miss_D=round(abs(chosen - hindsight), 4),
        ))

    n = len(results)
    same = sum(1 for x in results if abs(x['power_difference_D']) < 1e-9)
    differ = n - same
    twin_closer = sum(1 for x in results if x['twin_power_miss_D'] < x['implanted_power_miss_D'] - 1e-9)
    surgeon_closer = sum(1 for x in results if x['implanted_power_miss_D'] < x['twin_power_miss_D'] - 1e-9)
    within_twin = sum(1 for x in results if x['twin_power_miss_D'] <= 0.5)
    within_real = sum(1 for x in results if x['implanted_power_miss_D'] <= 0.5)
    mae_twin = sum(x['twin_power_miss_D'] for x in results) / n
    mae_real = sum(x['implanted_power_miss_D'] for x in results) / n

    summary = dict(
        patients=n,
        target_refraction_D=TARGET_D,
        tolerance_D=TOL_D,
        same_power_as_surgeon=same,
        different_power=differ,
        max_abs_power_difference_D=max(abs(x['power_difference_D']) for x in results),
        twin_closer_to_target=twin_closer,
        surgeon_closer_to_target=surgeon_closer,
        twin_within_half_dioptre_of_hindsight=within_twin,
        implanted_within_half_dioptre_of_hindsight=within_real,
        mean_power_miss_twin_D=round(mae_twin, 4),
        mean_power_miss_implanted_D=round(mae_real, 4),
        bias_removal='leave-one-out over the other 19 patients; a patient never sees its own outcome',
        control=('the power actually implanted, scored against the same hindsight-correct power '
                 'derived from the same measured refraction'),
        inherited_scope=('mean-curvature scalar cornea with preoperative geometry frozen after surgery, '
                         'no astigmatic vector, modelled lens design rather than the implanted serial '
                         'number lens; see results/LANE_EYE_OPTICAL_TWIN/r15'),
        claim_type='capability',
        review_state='PENDING_INDEPENDENT_REVIEW',
    )

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'DECISION_V1.json').write_text(json.dumps(
        dict(summary=summary, rows=results), indent=1, ensure_ascii=False))
    for k, v in summary.items():
        if isinstance(v, (int, float)):
            print(f'  {k} = {v}')
    print(f'  written: {OUT / "DECISION_V1.json"}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
