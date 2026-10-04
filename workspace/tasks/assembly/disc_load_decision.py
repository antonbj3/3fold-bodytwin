"""A decision in the spine: what load keeps nucleus pressure under a limit, and what does the facit say?

Why the spine and why now. The detail-layer harvest put 99 edges into the net and raised its external
references from 15 to 27, and the strongest of them was already a documented failure rather than an
agreement. From `MECHANISM_INTERVERTEBRAL_DISC.md`, verified line by line in the source:

  disc area                  1800 mm^2
  nominal stress F/A         0.9045 MPa
  predicted nucleus pressure 1.176 MPa at k=1.3, 1.357 at k=1.5, 1.378 at the empirical k=1.524
  Wilke in-vivo walking      0.53 to 0.65 MPa
  verdict                    FAIL, both k land above the measured range, by 1.81x and 2.09x

So the model overpredicts the pressure inside a living human disc by roughly a factor of two, and the
source says so itself. That is a calibration gap with an external facit, which is exactly the shape the
laser fluence had, and there the two errors compounded once the chain was inverted.

The decision. A forward model that overpredicts pressure is not harmless: inverted, it prescribes a
load limit. Ask what load keeps the nucleus under a stated pressure, answer it with our multiplier and
with the multiplier the measurement implies, and report the difference in newtons rather than in percent
of a coefficient. Percent of a coefficient is what hid the laser compounding.

The implied multiplier is the quantity the facit actually pins: k_implied = P_measured / (F/A), which
from 0.53 and 0.65 MPa against 0.9045 MPa gives 0.586 to 0.719 -- less than half the empirical 1.3 to
1.5. Nothing is fitted here; both are read off.

The control is the regime our own model would prescribe, which is the equally informed comparison: same
geometry, same nominal stress, only the multiplier changes.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('')
OUT = W / 'results/ASSEMBLY_DISC_LOAD_DECISION'
SRC = 'source_documents/INTERVERTEBRAL_DISC.md'

A_DISC_MM2 = 1800.0
NOMINAL_STRESS_MPA = 0.9045          # F/A at the source's reference load
K_MODEL = (1.3, 1.5)                 # our empirical pressure multiplier range
K_EMPIRICAL_SINGLE = 1.524
WILKE_WALKING_MPA = (0.53, 0.65)     # in-vivo, held out


def implied_k() -> tuple[float, float]:
    return (WILKE_WALKING_MPA[0] / NOMINAL_STRESS_MPA,
            WILKE_WALKING_MPA[1] / NOMINAL_STRESS_MPA)


def load_for_pressure(limit_mpa: float, k: float) -> float:
    """The axial load in newtons that puts the nucleus at the limit, for a given multiplier."""
    # P = k * F / A, with A in mm^2 and P in MPa gives F in newtons directly
    return limit_mpa * A_DISC_MM2 / k


def main() -> int:
    ki_lo, ki_hi = implied_k()
    reference_load_n = NOMINAL_STRESS_MPA * A_DISC_MM2

    rows = []
    for limit in (0.5, 0.65, 1.0, 1.5, 2.0):
        ours_strict = load_for_pressure(limit, K_MODEL[1])      # our most conservative multiplier
        ours_loose = load_for_pressure(limit, K_MODEL[0])
        facit_strict = load_for_pressure(limit, ki_hi)
        facit_loose = load_for_pressure(limit, ki_lo)
        rows.append({
            'pressure_limit_MPa': limit,
            'load_our_model_N': [round(ours_strict, 0), round(ours_loose, 0)],
            'load_facit_implied_N': [round(facit_strict, 0), round(facit_loose, 0)],
            'our_limit_is_lower_by_N': round(facit_strict - ours_strict, 0),
            'our_limit_is_lower_by_factor': round(facit_strict / ours_strict, 3),
        })

    summary = {
        'question': ('what axial load keeps nucleus pressure under a limit, and how far is our answer '
                     'from the one the in-vivo measurement implies'),
        'source': {'path': SRC, 'all_values_read_off_not_fitted': True},
        'measured_inputs': {
            'disc_area_mm2': A_DISC_MM2,
            'nominal_stress_MPa': NOMINAL_STRESS_MPA,
            'reference_load_N': round(reference_load_n, 0),
            'our_multiplier_range': list(K_MODEL),
            'our_single_empirical_multiplier': K_EMPIRICAL_SINGLE,
            'wilke_in_vivo_walking_MPa': list(WILKE_WALKING_MPA),
        },
        'the_multiplier_the_facit_pins': {
            'k_implied_low': round(ki_lo, 4),
            'k_implied_high': round(ki_hi, 4),
            'our_k_over_implied_k': [round(K_MODEL[0] / ki_hi, 3), round(K_MODEL[1] / ki_lo, 3)],
            'reading': ('the measurement implies a multiplier of 0.59 to 0.72 where our model uses 1.3 '
                        'to 1.5, so ours is high by a factor of 1.8 to 2.6 on the coefficient itself'),
        },
        'decision_rows': rows,
        'what_the_overprediction_costs_as_a_decision': (
            'inverted, a multiplier that is twice too high prescribes a load limit that is about half '
            'what the measurement supports. Our model is therefore CONSERVATIVE in this direction: it '
            'under-permits load rather than over-permitting it, which is the safe side of a clinical '
            'error but makes the model useless for asking how much load is tolerable'),
        'control': ('our own multiplier on the same geometry and the same nominal stress; only the '
                    'coefficient changes, which is the equally informed comparison'),
        'falsifier': ('if the in-vivo range and our nominal stress are not measured at the same posture '
                      'and load, the implied multiplier is not comparable and this whole inversion is '
                      'void. The source states walking for the in-vivo range; our nominal stress is at '
                      'its own reference load, and whether those match is NOT established here'),
        'claim_type': 'information_link',
        'scope': 'arithmetic on values read from one document; no clinical validity is claimed',
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'DISC_LOAD_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    t = summary['the_multiplier_the_facit_pins']
    print(f"  multiplier: spring {K_MODEL[0]}-{K_MODEL[1]}, facit implicerar {t['k_implied_low']}-{t['k_implied_high']}, that is {t['our_k_over_implied_k']}x too high")
    print(f"\n  {'limit MPa':>10s} {'our load N':>22s} {'facit N':>22s} {'faktor':>7s}")
    for r in rows:
        print(f"  {r['pressure_limit_MPa']:>10.2f} {str(r['load_our_model_N']):>22s} "
              f"{str(r['load_facit_implied_N']):>22s} {r['our_limit_is_lower_by_factor']:>7.2f}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
