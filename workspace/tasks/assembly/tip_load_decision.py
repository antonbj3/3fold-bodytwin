#!/usr/bin/env python3
"""A decision outside the eye: what tip load may a robot-held instrument assume from shaft sensing?

Why this exists. A surgical robot has sensors in the shaft or the joint, never in the tip. Everything
the twin wants to say about tissue at the tip -- adhesion, cutting, damage -- has to be inferred from
the wrong end of the instrument. PROOF_LANE_TIP_FROM_SHAFT settled what is and is not recoverable, and this
file is that theorem turned into a decision that returns a number and refuses when it cannot.

The theorem, and the part of it this coordinator re-derived. For a shaft force channel

    y = A(T, h) + b*v + S

the shaft identifies only the SUM. Exchanging adhesion against other resistance preserves the signal
while changing the tip load, so no amount of force data separates them: 53 248 indistinguishable pair
cases at exact observation error zero, and one sharp example admitting 0.2 to 0.8 N of tip load at
identical shaft force, a factor of four. Temperature identifies adhesion CHANGES but not the absolute
offset. Speed, phase and depth each fail on their own, because adhesive memory reproduces a damping
phase exactly.

What does work is a second mechanical channel. With tip load L at projected lever arm l and the
nuisance resultant S at projected lever arm r, the shaft reads

    y = L + S,     M = l*L + r*S

and solving the pair gives L = (M - r*y) / (l - r). This coordinator re-derived that inversion
symbolically rather than taking it: solving the two equations returns exactly (M - r*y)/(l - r).

The refusal is in the same algebra. The inversion is singular when l = r, which is the axially
collinear case -- the tip load and the nuisance acting on the same projected lever arm. No calibration
fixes that; it is a geometry the instrument must not be in when a tip number is needed. This decision
returns no tip load when the lever arms are closer than a stated separation.

Why it does not depend on the tissue model. The inversion is mechanical: it recovers the total
resultant whatever the adhesion law turns out to be. That matters because the adhesion-dominance
conclusion our own lane tested does not transfer to our contact pressures -- 1273.2395447351628 kPa at
the default 0.7853981633974483 mm^2 patch under 1 N, which is 36.378273x the 35 kPa gel where the
conclusion was measured. A decision that needed the adhesion law would inherit that gap. This one does
not.

The control is Coulomb, and it does not escape either: y = mu*P + S carries the same ambiguity unless
the normal load or the nuisance resistance is independently known. The coefficient usually quoted is
static, and proof_lane notes its reported dispersion is not a rigorous uncertainty bound -- Urrea et al.,
J Mech Behav Biomed Mater 56:98-105 (2016), PMID 26700572, doi 10.1016/j.jmbbm.2015.11.024. An earlier
night-log row of mine treated 0.295 +/- 0.056 as if the dispersion bounded uncertainty; it does not,
and the sigma distances computed from it are descriptive only.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('')
OUT = W / 'results/ASSEMBLY_TIP_LOAD'

# Minimum lever-arm separation, as a fraction of the tip arm, below which the inversion is refused.
# Set from the conditioning of the inverse rather than from taste: at separation f the tip load error
# is amplified by 1/f relative to the moment error, so f = 0.1 means a tenfold amplification and
# anything smaller turns a millinewton-metre of moment noise into newtons of tip load.
MIN_ARM_SEPARATION_FRACTION = 0.1


def decide(shaft_force_N: float, shaft_moment_Nm: float, tip_arm_m: float, nuisance_arm_m: float,
           moment_noise_Nm: float = 1e-3, force_noise_N: float = 1e-2) -> dict:
    out = {'shaft_force_N': shaft_force_N, 'shaft_moment_Nm': shaft_moment_Nm,
           'tip_arm_m': tip_arm_m, 'nuisance_arm_m': nuisance_arm_m}
    sep = tip_arm_m - nuisance_arm_m
    if tip_arm_m <= 0:
        out['decision'] = 'NO_TIP_LOAD_RETURNED'
        out['reason'] = 'tip lever arm must be positive'
        return out
    frac = abs(sep) / tip_arm_m
    out['arm_separation_fraction'] = frac
    if frac < MIN_ARM_SEPARATION_FRACTION:
        out['decision'] = 'NO_TIP_LOAD_RETURNED'
        out['reason'] = (f'lever arms separated by {frac:.4f} of the tip arm, below '
                         f'{MIN_ARM_SEPARATION_FRACTION}; the inversion L = (M - r*y)/(l - r) is '
                         f'singular at equal arms, which is the axially collinear case, and no '
                         f'calibration removes it')
        return out
    L = (shaft_moment_Nm - nuisance_arm_m * shaft_force_N) / sep
    # Propagated as an interval, both endpoints, not as a width: the endpoint rule.
    amp_M, amp_y = 1.0 / abs(sep), abs(nuisance_arm_m) / abs(sep)
    half = amp_M * moment_noise_Nm + amp_y * force_noise_N
    out.update({
        'decision': 'TIP_LOAD_N',
        'decision_value_N': L,
        'tip_load_enclosure_N': [L - half, L + half],
        'enclosure_half_width_N': half,
        'moment_noise_amplification_per_m': amp_M,
        'force_noise_amplification': amp_y,
        'nuisance_resultant_N': shaft_force_N - L,
        'control_coulomb': ('y = mu*P + S carries the same ambiguity unless the normal load or the '
                            'nuisance is independently known; the quoted coefficient is static and '
                            'its dispersion is not a rigorous bound (PMID 26700572)'),
        'independent_of_adhesion_law': True,
    })
    return out


def main() -> int:
    cases = [
        # tip arm, nuisance arm: a grasper with the nuisance contact near the trocar, then progressively
        # closer to the tip, ending in the collinear case the theorem excludes.
        (1.0, 0.020, 0.150, 0.010),
        (1.0, 0.020, 0.150, 0.075),
        (1.0, 0.020, 0.150, 0.135),
        (1.0, 0.020, 0.150, 0.145),
    ]
    rows = [decide(*c) for c in cases]
    summary = {
        'question': 'what tip load may be assumed from shaft force and shaft bending moment',
        'claim_type': 'capability',
        'theorem': ('y = L + S and M = l*L + r*S give L = (M - r*y)/(l - r), re-derived symbolically '
                    'by this coordinator; singular at l = r, the axially collinear case'),
        'why_force_alone_fails': ('the shaft identifies only the sum A(T,h) + b*v + S; 53248 '
                                  'indistinguishable pair cases at exact observation error zero, and '
                                  'one sharp example admits 0.2 to 0.8 N of tip load at identical '
                                  'shaft force, a factor of 4'),
        'what_temperature_buys': 'adhesion changes, not the absolute offset',
        'what_fails_alone': 'speed, phase and depth; adhesive memory reproduces a damping phase exactly',
        'refusal': {'min_arm_separation_fraction': MIN_ARM_SEPARATION_FRACTION,
                    'reason': 'below it the moment noise amplification 1/(l-r) exceeds tenfold'},
        'control': 'Coulomb y = mu*P + S, same ambiguity; PMID 26700572, doi 10.1016/j.jmbbm.2015.11.024',
        # Quoted by the edge, so emitted here. The dispersion is the source's reported scatter and
        # NOT a rigorous uncertainty bound, which is why no sigma distance is computed from it.
        'coulomb_control_coefficient': {
            'mu_static': 0.295, 'mu_static_reported_scatter': 0.056,
            'mu_dynamic': 0.255, 'mu_dynamic_reported_scatter': 0.086,
            'source': 'Urrea et al., J Mech Behav Biomed Mater 56:98-105 (2016), PMID 26700572',
            'dispersion_is_not_an_uncertainty_bound': True},
        'falsifier': ('if a shaft-only observable is found that separates adhesion from other '
                      'resistance without a second mechanical channel, the theorem is wrong and this '
                      'decision is unnecessary'),
        'geometry_swing_at_identical_shaft_reading': {
            'note': ('the two admissible cases below carry the SAME shaft force and the SAME shaft '
                     'moment and differ only in where the nuisance contact sits'),
            'nuisance_arm_m': [0.010, 0.075],
            'tip_load_N': [rows[0]['decision_value_N'], rows[1]['decision_value_N']],
            'swing_N': rows[0]['decision_value_N'] - rows[1]['decision_value_N'],
            'reading': ('a swing of this size at one shaft reading is the same order as the 0.2 to '
                        '0.8 N ambiguity proof_lane exhibited, so the lever arms are not a detail of the '
                        'calibration -- they are an input the decision cannot do without'),
            'sign_changes': (rows[0]['decision_value_N'] > 0) != (rows[1]['decision_value_N'] > 0)},
        'boundary_note': ('the third case has separation fraction 0.09999999999999999 in floating '
                          'point and is refused at exactly the stated threshold; the refusal is '
                          'reported rather than smoothed'),
        'rows': rows,
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
        'no_claim_of_biological_validation': True,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'DECISION_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f'{"y N":>6} {"M Nm":>7} {"l m":>7} {"r m":>7} {"sep":>7} {"beslut":>20} {"L N":>9} {"holje +/- N":>12}')
    for r in rows:
        if r['decision'] == 'NO_TIP_LOAD_RETURNED':
            print(f'{r["shaft_force_N"]:6} {r["shaft_moment_Nm"]:7} {r["tip_arm_m"]:7} '
                  f'{r["nuisance_arm_m"]:7} {r.get("arm_separation_fraction", 0):7.4f} '
                  f'{"INGET SVAR":>20}')
        else:
            print(f'{r["shaft_force_N"]:6} {r["shaft_moment_Nm"]:7} {r["tip_arm_m"]:7} '
                  f'{r["nuisance_arm_m"]:7} {r["arm_separation_fraction"]:7.4f} '
                  f'{"spetslast":>20} {r["decision_value_N"]:9.4f} {r["enclosure_half_width_N"]:12.4f}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
