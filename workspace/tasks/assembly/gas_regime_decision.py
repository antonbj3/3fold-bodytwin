"""A decision outside the eye: does a measured gas output report production, or only the overflow?

Where the numbers come from. A creative brief harvested a complete measured gas balance for the colon
into `results/LANE_READ_CREATIVE/CONSUMABLE.json`, every row marked `measurement` with a locator: total
flatus 705 mL/24h median, hydrogen 361, carbon dioxide 68, an unidentified residual of 213, day and
sleep rates of 34 and 16 mL/h, a fibre-free total of 214, a hydrogen absorption fraction spanning 90 to
20 per cent, and the amount below which essentially all hydrogen is absorbed, 76 mL/6h.

Why that last number makes a decision possible. An absorption fraction that swings by 4.5x is not noise
around a mean -- it is two regimes. Below the threshold the gut absorbs essentially everything, so a
measured output reports almost nothing about production. Above it the excess passes through, and the
output tracks production with a known loss. Any inference from a measured output therefore depends on
which regime the subject is in, and that question is answerable from the measurement itself.

The decision: given a measured 6-hour hydrogen output, state the regime, and give the interval of
production consistent with it. The honest answer in the absorbing regime is a LOWER BOUND and a wide
interval, not a number, and saying so is the point.

The control is the regime-blind reading: taking the measured output as proportional to production with
one fixed absorption fraction, which is what a single-coefficient model does. The facit is the measured
pair of fractions and the measured threshold, none of which this cell fits.

Mass closure is checked rather than assumed, because the residual is 30 per cent of the total and an
unexamined 30 per cent is where a conservation error hides.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('.')
OUT = W / 'results/ASSEMBLY_GAS_REGIME_DECISION'

TOTAL_24H = 705.0          # mL/24h, median
H2_24H = 361.0
CO2_24H = 68.0
RESIDUAL_24H = 213.0
METHANE_OBSERVED = [3.0, 26.0, 120.0]      # mL/24h, in three of ten individuals
RATE_DAY, RATE_SLEEP = 34.0, 16.0          # mL/h
FIBRE_FREE_TOTAL = 214.0
ABSORB_HIGH, ABSORB_LOW = 90.0, 20.0       # per cent, the two ends of the measured range
THRESHOLD_6H = 76.0                        # mL/6h below which essentially all hydrogen is absorbed


def regime(measured_6h: float) -> str:
    if measured_6h < THRESHOLD_6H:
        return 'ABSORBING'
    return 'OVERFLOW'


def production_interval(measured_6h: float) -> tuple[float, float | None, str]:
    """Production consistent with a measured output, given the two measured absorption fractions."""
    if regime(measured_6h) == 'ABSORBING':
        # Essentially all of it is absorbed, so the output is a floor on production and the ceiling is
        # whatever the threshold permits before overflow would have been seen.
        lo = measured_6h / (1.0 - ABSORB_LOW / 100.0)
        return (lo, None, 'lower bound only: the output does not report production in this regime')
    lo = measured_6h / (1.0 - ABSORB_LOW / 100.0)
    hi = measured_6h / (1.0 - ABSORB_HIGH / 100.0)
    return (lo, hi, 'bracketed by the two measured absorption fractions')


def main() -> int:
    closure = TOTAL_24H - (H2_24H + CO2_24H + RESIDUAL_24H)
    rows = []
    for measured in (20.0, 60.0, 76.0, 100.0, 200.0, 400.0):
        lo, hi, note = production_interval(measured)
        # the regime-blind control: one fixed fraction, the midpoint of the measured range
        mid = (ABSORB_HIGH + ABSORB_LOW) / 2.0
        blind = measured / (1.0 - mid / 100.0)
        rows.append({
            'measured_6h_mL': measured,
            'regime': regime(measured),
            'production_low_mL_6h': round(lo, 1),
            'production_high_mL_6h': round(hi, 1) if hi is not None else None,
            'interval_width_mL': round(hi - lo, 1) if hi is not None else None,
            'regime_blind_single_fraction_mL': round(blind, 1),
            'blind_error_against_low_end_percent': round(100.0 * (blind - lo) / lo, 1),
            'note': note,
        })

    summary = {
        'question': ('does a measured gas output report production or only the overflow, and what '
                     'production interval is consistent with it'),
        'why_a_decision_exists': ('the measured absorption fraction spans 90 to 20 per cent, a 4.5-fold '
                                  'swing, and a measured threshold of 76 mL/6h separates the two '
                                  'regimes. That is two regimes, not scatter about a mean'),
        'measured_inputs': {
            'total_flatus_mL_24h': TOTAL_24H, 'hydrogen_mL_24h': H2_24H,
            'carbon_dioxide_mL_24h': CO2_24H, 'unidentified_residual_mL_24h': RESIDUAL_24H,
            'methane_in_three_of_ten_mL_24h': METHANE_OBSERVED,
            'rate_day_mL_h': RATE_DAY, 'rate_sleep_mL_h': RATE_SLEEP,
            'fibre_free_total_mL_24h': FIBRE_FREE_TOTAL,
            'absorption_fraction_percent_range': [ABSORB_HIGH, ABSORB_LOW],
            'threshold_mL_6h': THRESHOLD_6H,
            'all_marked_measurement_in': 'results/LANE_READ_CREATIVE/CONSUMABLE.json',
        },
        'mass_closure': {
            'total_minus_named_components_mL_24h': round(closure, 1),
            'residual_as_fraction_of_total_percent': round(100.0 * RESIDUAL_24H / TOTAL_24H, 1),
            'reading': ('the named components leave ' + f'{closure:.0f}' + ' mL/24h unaccounted beyond '
                        'the residual already listed, and the residual is itself 30 per cent of the '
                        'total. Nitrogen and oxygen are the expected remainder and are not in the '
                        'harvested rows, so closure is stated as open rather than claimed'),
        },
        'rate_consistency_check': {
            'day_plus_sleep_over_24h_mL': round(RATE_DAY * 16 + RATE_SLEEP * 8, 1),
            'measured_total_mL_24h': TOTAL_24H,
            'reading': ('a 16-hour waking day at the day rate plus 8 hours at the sleep rate gives a '
                        'figure to compare against the measured total; the two come from the same '
                        'source, so agreement is a consistency check and not evidence'),
        },
        'decision_rows': rows,
        'control': ('the regime-blind reading: one fixed absorption fraction at the midpoint of the '
                    'measured range, which is what a single-coefficient model does'),
        'claim_type': 'capability',
        'scope': ('arithmetic on published measured values; no biological or clinical validity is '
                  'claimed and no subject is modelled'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'GAS_REGIME_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    print(f'  mass balance: named components leave {closure:.0f} mL/24h unaccounted for, the residual is {100 * RESIDUAL_24H / TOTAL_24H:.0f} % av totalen')
    rc = summary['rate_consistency_check']
    print(f"  takt: {rc['day_plus_sleep_over_24h_mL']} mL/24h from day and night on a measured basis {rc['measured_total_mL_24h']}")
    print(f"\n  {'measured 6h':>8s} {'regim':>10s} {'produktion lo':>14s} {'hi':>8s} {'regimblind':>11s} {'fel %':>7s}")
    for r in rows:
        print(f"  {r['measured_6h_mL']:>8.0f} {r['regime']:>10s} {r['production_low_mL_6h']:>14.1f} "
              f"{str(r['production_high_mL_6h']):>8s} {r['regime_blind_single_fraction_mL']:>11.1f} "
              f"{r['blind_error_against_low_end_percent']:>7.1f}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
