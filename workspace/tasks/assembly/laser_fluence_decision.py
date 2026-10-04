"""A decision outside the eye, scored against a published law: what fluence removes the mass you want.

Why this and why now. Every decision the twin makes was in the eye, and I had taken that for a property
of the project when it is a property of where I looked: the constraint net holds 12 edges outside the
eye that carry an external reference against 3 inside it. The laser line is the sharpest of them,
because its external anchor is a LAW and not a single value -- a published mass-loss regression with
both an intercept and a slope.

  published (held out)   threshold 1.15 J/cm^2, slope 267 ug/J
  our chain              threshold 1.4548 J/cm^2, slope 200.33 ug/J

So the chain is 26.51 % high on the threshold and 24.97 % LOW on the slope, and the two errors have
opposite signs. That is the interesting part and it is not visible from either number alone: a
threshold set too high removes too little, a slope set too low asks for too much fluence to compensate,
and whether they cancel depends on how much mass you are trying to remove. Nobody has asked what that
does to a decision, because the comparison was reported as two relative errors rather than inverted.

The decision. Given a target mass removed per unit area, choose the fluence. Ours chooses
F = our_threshold + target / our_slope. The published law then says what that fluence actually removes.
The error is in micrograms per square centimetre of tissue, which is the unit the decision is about,
rather than in percent of a regression coefficient.

Controls, both equally informed.
  1. The published law used as the decision rule, which is the best anyone could do with the outside
     literature and no model of ours. This is the control that matters.
  2. A fixed fluence at the middle of the range, which is what a setting on a device does. A decision
     that cannot beat a constant is not a decision.

A convention is in the way and is reported rather than resolved: the same lane computes the threshold as
1.4548 J/cm^2 for a Gaussian beam proxy and 1.9081 for a top-hat, a 31 % spread from the beam profile
alone. Both are carried through below, because picking one silently would hide that the beam convention
matters more than the slope error.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('')
SRC = W / 'results/LANE_LASER_SURGERY/r17/EXTERNAL_COMPARISON_V1.json'
OUT = W / 'results/ASSEMBLY_LASER_FLUENCE_DECISION'

HELD_THRESHOLD, HELD_SLOPE = 1.15, 267.0          # J/cm^2 and ug/J, published, held out
TARGETS_UG_CM2 = [25.0, 50.0, 100.0, 200.0, 400.0, 800.0]


def main() -> int:
    comps = json.loads(SRC.read_text())['comparisons']
    # Key on criterion AND profile. The file holds six comparisons, three criteria x two beam
    # profiles, and keying on profile alone silently kept the LAST criterion -- a boiling-onset
    # variant with a NEGATIVE threshold of -0.765 J/cm^2. That produced a chosen fluence below the
    # published threshold and a flat -100 % error, which is what sent me back to the file. A negative
    # ablation threshold is physically impossible and was the tell.
    ours = {(c['criterion'], c['profile']): (c['predicted_intercept_J_cm2'],
                                             c['predicted_slope_microgram_per_J'])
            for c in comps if 'predicted_slope_microgram_per_J' in c}
    ours = {prof: v for (crit, prof), v in ours.items()
            if crit == 'complete_vaporization_primary'}
    assert all(thr > 0 for thr, _ in ours.values()), 'a threshold must be positive'

    rows = []
    for profile, (thr, slope) in sorted(ours.items()):
        for target in TARGETS_UG_CM2:
            f_ours = thr + target / slope
            f_held = HELD_THRESHOLD + target / HELD_SLOPE
            # what the published law says our chosen fluence actually removes
            delivered = HELD_SLOPE * max(f_ours - HELD_THRESHOLD, 0.0)
            rows.append({
                'profile': profile,
                'target_ug_cm2': target,
                'our_threshold_J_cm2': round(thr, 4),
                'our_slope_ug_per_J': round(slope, 2),
                'fluence_we_would_choose_J_cm2': round(f_ours, 4),
                'fluence_the_published_law_would_choose_J_cm2': round(f_held, 4),
                'mass_actually_removed_under_published_law_ug_cm2': round(delivered, 2),
                'mass_error_ug_cm2': round(delivered - target, 2),
                'mass_error_percent': round(100.0 * (delivered - target) / target, 2),
            })

    # Control 2: one fixed fluence for every target, chosen as the middle of what the published law
    # would need across the target range.
    fixed = HELD_THRESHOLD + (min(TARGETS_UG_CM2) + max(TARGETS_UG_CM2)) / 2 / HELD_SLOPE
    fixed_rows = [{'target_ug_cm2': t,
                   'mass_removed_ug_cm2': round(HELD_SLOPE * (fixed - HELD_THRESHOLD), 2),
                   'mass_error_percent': round(100.0 * (HELD_SLOPE * (fixed - HELD_THRESHOLD) - t) / t, 2)}
                  for t in TARGETS_UG_CM2]

    gauss = [r for r in rows if r['profile'] == 'gaussian_proxy']
    summary = {
        'question': 'what fluence removes the mass you want, and what does our calibration error cost',
        'external_anchor': {'threshold_J_cm2': HELD_THRESHOLD, 'slope_ug_per_J': HELD_SLOPE,
                            'kind': 'published mass-loss regression, held out, both coefficients'},
        'our_calibration_error_as_reported': {'threshold_relative': 0.2651, 'slope_relative': -0.2497,
                                              'note': 'opposite signs, so the effect on a decision '
                                                      'depends on the target and was never inverted'},
        'beam_profile_convention': {
            'gaussian_proxy_threshold_J_cm2': ours.get('gaussian_proxy', (None,))[0],
            'top_hat_threshold_J_cm2': ours.get('top_hat_diagnostic', (None,))[0],
            'note': 'the profile convention moves the threshold more than the slope error does, and is '
                    'reported rather than resolved here'},
        'decision_rows': rows,
        'control_fixed_fluence': {'fluence_J_cm2': round(fixed, 4), 'rows': fixed_rows},
        'control_published_law_as_rule': 'exact by construction, it is the facit',
        'claim_type': 'capability',
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
        'scope': ('a decision against a published law, not against measured tissue outcomes; no claim '
                  'of biological or clinical validity'),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'LASER_DECISION_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    print(f'  published law: threshold {HELD_THRESHOLD} J/cm2, lutning {HELD_SLOPE} ug/J\n')
    print(f"  {'target ug/cm2':>11s} {'our fluence':>11s} {'publicerad':>11s} {'faktiskt bort':>13s} {'fel %':>8s}   fast fluens fel %")
    for r, fr in zip(gauss, fixed_rows):
        print(f"  {r['target_ug_cm2']:>11.0f} {r['fluence_we_would_choose_J_cm2']:>11.4f} "
              f"{r['fluence_the_published_law_would_choose_J_cm2']:>11.4f} "
              f"{r['mass_actually_removed_under_published_law_ug_cm2']:>13.1f} "
              f"{r['mass_error_percent']:>8.1f}   {fr['mass_error_percent']:>+8.1f}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
