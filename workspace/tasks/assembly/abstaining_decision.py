"""A decision that declines: which eyes cannot be decided to the clinical threshold, said beforehand.

CORRECTION, from auditing this cell against its own headline. What it composes is NOT a floor, and I
reported it as one. A floor bounds every unit; 50 of 89 eyes have a realised miss below their own
composed value, averaging 0.2035 D against 0.4508, which is what a composed standard deviation does.
The median was the tell and I read past it: realised 0.3675 D against composed 0.4794 D. The honest
claim is that input uncertainty is the same magnitude as the residual error, so there is little
headroom left -- not that 0.08 D of model error remains.

This is the capability the failed probability work was reaching for and could not get. Scoring the
decision as a probability gave a Brier skill of +0.038 and -0.085 against the cohort rate, with a
resolution of 0.0034, because the uncertainty was fitted from outcomes and an additive split of a
once-observed total is not identifiable. Three measurements since then changed that: the floor terms are
now known per eye and measured independently of the outcome.

  label         ISO 11979-2 permits +-0.30 to +-1.00 D by power band, and measured deviations run
                0.08 D below 20 D and 0.24 D at or above it
  lens position predicted to 0.1517 mm of spread, and the chain converts position to refraction at a
                per-eye sensitivity measured from that eye's own biometry
  cornea        the two instruments in the workbook disagree on mean corneal power with a spread of
                0.3610 D on the same 89 eyes

A convention had to be separated from a disagreement first, for the third time tonight. The raw
difference in mean corneal power between the devices is +0.8511 D systematic, and that is the known gap
between keratometric power at an assumed index and true net corneal power including the posterior
surface, not a measurement conflict. Using the offset would have inflated every eye's floor by the
better part of a dioptre. The SPREAD, 0.3610 D, is the disagreement term.

What is new here. Every term is known before the operation and none comes from the result, so the
decision can refuse in advance: this eye's inputs do not support a 0.25 D answer. It is then scored
against what actually happened, which is the test a confidence claim has to pass. The control is the
cohort: deciding every eye.
"""
from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path

W = Path('.')
ENGINE = Path('the public staging tree/3fold-graph-engine/src')
OUT = W / 'results/ASSEMBLY_ABSTAINING_DECISION'
THRESHOLD_D = 0.25

POSITION_SPREAD_MM = 0.1517     # measured: leave-one-out prediction spread
CORNEA_SPREAD_D = 0.3610        # measured: between-device spread in mean corneal power, offset removed
LABEL_LOW_D, LABEL_HIGH_D = 0.08, 0.24      # measured deviations below and at/above 20 D

sys.path.insert(0, str(ENGINE))
from graph_engine.stage_decorrelation_verifier import covariance_aware_margin   # noqa: E402


def main() -> int:
    rows = json.loads((W / 'results/ASSEMBLY_PROSPECTIVE_POWER/PROSPECTIVE_V1.json').read_text())['rows']

    eyes, stages = [], []
    for r in rows:
        label = LABEL_HIGH_D if abs(r['implanted_power_D']) >= 20 else LABEL_LOW_D
        position = POSITION_SPREAD_MM * abs(r['dR_dAQD_D_per_mm'])
        cornea = CORNEA_SPREAD_D
        stages.append([label, position, cornea])
        floor_rss = math.sqrt(label ** 2 + position ** 2 + cornea ** 2)
        floor_sum = label + position + cornea
        miss = abs(r['power_with_predicted_position_D'] - r['hindsight_correct_power_D'])
        eyes.append({
            'sheet_row': r['sheet_row'],
            'label_term_D': round(label, 4),
            'position_term_D': round(position, 4),
            'cornea_term_D': cornea,
            'floor_if_independent_D': round(floor_rss, 4),
            'floor_if_fully_shared_D': round(floor_sum, 4),
            'decidable_to_threshold': floor_rss <= THRESHOLD_D,
            'realised_miss_D': round(miss, 4),
        })

    import numpy as np
    bracket = covariance_aware_margin(np.asarray(stages), THRESHOLD_D, ci_pctl=99.0,
                                      n_boot=300, seed=0, min_units=20)

    n = len(eyes)
    abstain = [e for e in eyes if not e['decidable_to_threshold']]
    decide = [e for e in eyes if e['decidable_to_threshold']]
    miss_all = [e['realised_miss_D'] for e in eyes]

    summary = {
        'question': 'can the chain say in advance which eyes its inputs cannot decide to 0.25 D',
        'eyes': n,
        'terms_and_where_each_was_measured': {
            'label': f'{LABEL_LOW_D} D below 20 D, {LABEL_HIGH_D} D at or above, from the external '
                     f'search; ISO 11979-2 permits more',
            'position': f'{POSITION_SPREAD_MM} mm of prediction spread times this eye own dR/dAQD',
            'cornea': f'{CORNEA_SPREAD_D} D, the between-device SPREAD in mean corneal power on these '
                      f'same 89 eyes; the +0.8511 D systematic offset is a keratometric index '
                      f'convention and is deliberately excluded',
        },
        'floor_D': {
            'median_if_independent': round(statistics.median(
                e['floor_if_independent_D'] for e in eyes), 4),
            'median_if_fully_shared': round(statistics.median(
                e['floor_if_fully_shared_D'] for e in eyes), 4),
            'module_bracket': {k: bracket.get(k) for k in
                               ('verdict', 'margin', 'rss_margin', 'sum_margin')},
        },
        'decision': {
            'eyes_declared_undecidable_to_0_25_D': len(abstain),
            'eyes_declared_decidable': len(decide),
            'of': n,
        },
        'scored_against_what_happened': {
            'mean_miss_all_eyes_D': round(statistics.mean(miss_all), 4),
            'mean_miss_among_declared_decidable_D': round(
                statistics.mean(e['realised_miss_D'] for e in decide), 4) if decide else None,
            'mean_miss_among_abstained_D': round(
                statistics.mean(e['realised_miss_D'] for e in abstain), 4) if abstain else None,
        },
        'control': ('the cohort: deciding every eye. An abstention is only worth anything if the eyes it '
                    'declines carry a materially larger realised miss than the ones it keeps'),
        'what_is_new': ('every term is measured before the operation and none comes from the outcome, so '
                        'the refusal is prospective; the earlier probability attempt fitted its spread '
                        'from outcomes and reached a resolution of 0.0034'),
        'claim_type': 'capability',
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
        'eyes_detail': eyes,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'ABSTAIN_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=str))

    s = summary['scored_against_what_happened']
    print(f"  eyes {n}")
    print(f"  floor: median {summary['floor_D']['median_if_independent']} D if independent, "
          f"{summary['floor_D']['median_if_fully_shared']} D if fully shared")
    print(f"  declared undecidable to {THRESHOLD_D} D: {len(abstain)} of {n}")
    print(f"  realised miss: all {s['mean_miss_all_eyes_D']} | kept "
          f"{s['mean_miss_among_declared_decidable_D']} | declined {s['mean_miss_among_abstained_D']} D")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
