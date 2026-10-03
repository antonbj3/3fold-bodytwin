#!/usr/bin/env python3
"""Can the twin say in advance WHICH eyes it cannot decide precisely for?

Why this is the useful question now. Measured 2026-10-03 from 900 exams over 300 eyes, three per eye:
the instrument's own repeat range exceeds the 0.25 D clinical threshold in 43 of 300 eyes and 0.50 D in
6. For roughly a seventh of patients the target precision is finer than the measurement can resolve, so
no model can win there. The richer corneal representation was separately measured to give zero clinical
gain, which is consistent with part of the error living in the input rather than the model.

So the remaining value is not more precision. It is telling a surgeon, before the operation, that THIS
eye is one where the decision cannot be trusted to 0.25 D. That is a capability nobody has from a single
exam, and it is testable with the data on disk.

The circularity that had to be avoided. Predicting an eye's KrSEQ repeat range from its own KrSEQ
spread is not prediction: both come from the same three exams. So the features here are the OTHER
measures of the same exam set and the zone difference, never the target's own spread. Stated plainly
because the first version of the surgical decision work was circular in exactly this way and had to be
withdrawn.

The control is the base rate. Always answering "reliable" is right 257 times in 300, so accuracy is a
worthless score here. What matters is whether a flagged subgroup carries a materially higher rate of
exceeding the threshold than the cohort does.
"""
from __future__ import annotations

import json
import statistics
from pathlib import Path

W = Path('.')
SRC = W / 'results/LANE_CORNEA_SHAPE/r3/REPEATABILITY_V1.json'
OUT = W / 'results/ASSEMBLY_RELIABILITY_FLAG'
THRESHOLD_D = 0.25


def main() -> int:
    data = json.loads(SRC.read_text())
    rows = data['records']

    # Target: does this eye's primary corneal power measure repeat worse than the clinical threshold?
    # Features: only quantities that do NOT come from the target's own three-exam spread.
    eyes = []
    for r in rows:
        kr = r['KrSEQ']
        f3, f6, zone = r['Kf3rSEQ'], r['Kf6rSEQ'], r['zone6_minus_zone3_D']
        eyes.append({
            'eye': r['anonymous_eye'],
            'target_range_D': kr['range_D'],
            'unreliable': kr['range_D'] > THRESHOLD_D,
            # features available from the same exam set but not from KrSEQ's own spread
            'f3_sd_D': f3['sample_sd_D'],
            'f6_sd_D': f6['sample_sd_D'],
            'zone_sd_D': zone['sample_sd_D'],
            'zone_mean_D': zone['mean_D'],
            'power_mean_D': kr['mean_D'],
            'f3_minus_f6_D': f3['mean_D'] - f6['mean_D'],
        })

    n = len(eyes)
    positives = sum(1 for e in eyes if e['unreliable'])
    base_rate = positives / n

    # One interpretable rule, not a fitted model: flag an eye when the OTHER measures of the same exam
    # set are themselves unstable. The threshold is chosen on the feature's own distribution (its upper
    # quartile), not tuned against the target, so this is a declared rule rather than a fit.
    results = {}
    for feat in ('f3_sd_D', 'f6_sd_D', 'zone_sd_D'):
        values = sorted(e[feat] for e in eyes)
        cut = values[int(0.75 * n)]
        flagged = [e for e in eyes if e[feat] >= cut]
        hits = sum(1 for e in flagged if e['unreliable'])
        unflagged = [e for e in eyes if e[feat] < cut]
        miss = sum(1 for e in unflagged if e['unreliable'])
        results[feat] = {
            'cut_at_upper_quartile_D': round(cut, 6),
            'flagged': len(flagged),
            'unreliable_among_flagged': hits,
            'rate_in_flagged': round(hits / len(flagged), 4) if flagged else None,
            'rate_in_unflagged': round(miss / len(unflagged), 4) if unflagged else None,
            'lift_over_base_rate': round((hits / len(flagged)) / base_rate, 3) if flagged else None,
            'sensitivity': round(hits / positives, 4) if positives else None,
        }

    # And the combined rule: flag when ANY of the three other-measure spreads is in its upper quartile.
    cuts = {f: sorted(e[f] for e in eyes)[int(0.75 * n)] for f in ('f3_sd_D', 'f6_sd_D', 'zone_sd_D')}
    flagged = [e for e in eyes if any(e[f] >= cuts[f] for f in cuts)]
    hits = sum(1 for e in flagged if e['unreliable'])
    unflagged = [e for e in eyes if e not in flagged]
    combined = {
        'flagged': len(flagged),
        'unreliable_among_flagged': hits,
        'rate_in_flagged': round(hits / len(flagged), 4) if flagged else None,
        'rate_in_unflagged': round(sum(1 for e in unflagged if e['unreliable']) / len(unflagged), 4)
                             if unflagged else None,
        'lift_over_base_rate': round((hits / len(flagged)) / base_rate, 3) if flagged else None,
        'sensitivity': round(hits / positives, 4) if positives else None,
    }

    summary = {
        'eyes': n,
        'unreliable_eyes': positives,
        'base_rate': round(base_rate, 4),
        'threshold_D': THRESHOLD_D,
        'single_feature_rules': results,
        'combined_any_upper_quartile': combined,
        'circularity_avoided': ("the target is KrSEQ's repeat range; no feature uses KrSEQ's own spread, "
                                "only the other measures of the same exam set"),
        'control': ("the base rate. Always answering reliable is correct in "
                    f"{n - positives} of {n} eyes, so accuracy is not a usable score; the question is "
                    "whether a flagged subgroup carries a materially higher rate"),
        'scope': ("a different instrument and cohort from the 89 surgical eyes; this establishes whether "
                  "unreliability is PREDICTABLE at all, not a flag for those patients"),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'RELIABILITY_V1.json').write_text(json.dumps(
        {'summary': summary, 'eyes': eyes}, indent=1, ensure_ascii=False))

    print(f"  eyes {n}  unreliable {positives}  base rate {base_rate:.4f}")
    for f, r in results.items():
        print(f"  {f:12s} flagged {r['flagged']:3d}  rate {r['rate_in_flagged']}  "
              f"lift {r['lift_over_base_rate']}  sens {r['sensitivity']}")
    print(f"  combined     flagged {combined['flagged']:3d}  rate {combined['rate_in_flagged']}  "
          f"lift {combined['lift_over_base_rate']}  sens {combined['sensitivity']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
