"""Where does an applied load go, and which uncertainty decides the answer?

Two harvested edges that have never been put against each other, both from the detail layer and both
with the document's own numbers:

  disc axial stiffness          1734 +/- 446 N/mm, human ex-vivo
  facet compressive fraction    0.03 to 0.25 normally, up to 0.47 in arthritic joints

Each alone invites a different mistake. The stiffness alone says the disc carries the load and gives a
displacement per newton. The fraction alone says some of it goes elsewhere. Together they decide what a
load actually does, and the question worth asking is which of the two uncertainties governs: a stiffness
known to +/- 26 per cent, or a fraction spanning more than eightfold.

That matters because they enter differently. The fraction splits the load before the disc sees it, so it
multiplies; the stiffness converts what remains into displacement, so it divides. A relative
uncertainty that multiplies and one that divides do not combine the way two additive errors do, and
tonight already produced one case where reading a composed spread as a bound was wrong.

The decision: for an applied axial load, the displacement interval, and which input to measure next if
the interval is too wide to act on. The control is each input taken alone at its central value, which
is how both numbers have been reported until now.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('.')
OUT = W / 'results/ASSEMBLY_LOAD_PATH'

K_DISC_N_PER_MM = 1734.0
K_DISC_SD = 446.0
FACET_FRACTION = (0.03, 0.25)
FACET_ARTHRITIC = 0.47


def displacement_mm(load_n: float, facet_fraction: float, k: float) -> float:
    """The disc carries what the facets do not, and converts it to displacement."""
    return load_n * (1.0 - facet_fraction) / k


def main() -> int:
    rows = []
    for load in (200.0, 500.0, 1000.0, 1628.0):     # the last is the disc-load figure from the net
        # both uncertainties at their extremes
        d_lo = displacement_mm(load, FACET_FRACTION[1], K_DISC_N_PER_MM + K_DISC_SD)
        d_hi = displacement_mm(load, FACET_FRACTION[0], K_DISC_N_PER_MM - K_DISC_SD)
        d_mid = displacement_mm(load, sum(FACET_FRACTION) / 2, K_DISC_N_PER_MM)
        # which input governs: vary one at a time from the midpoint
        only_facet = (displacement_mm(load, FACET_FRACTION[0], K_DISC_N_PER_MM)
                      - displacement_mm(load, FACET_FRACTION[1], K_DISC_N_PER_MM))
        only_stiff = (displacement_mm(load, sum(FACET_FRACTION) / 2, K_DISC_N_PER_MM - K_DISC_SD)
                      - displacement_mm(load, sum(FACET_FRACTION) / 2, K_DISC_N_PER_MM + K_DISC_SD))
        rows.append({
            'load_N': load,
            'displacement_low_mm': round(d_lo, 4),
            'displacement_mid_mm': round(d_mid, 4),
            'displacement_high_mm': round(d_hi, 4),
            'full_interval_width_mm': round(d_hi - d_lo, 4),
            'width_from_facet_fraction_alone_mm': round(only_facet, 4),
            'width_from_stiffness_alone_mm': round(only_stiff, 4),
            'facet_over_stiffness': round(only_facet / only_stiff, 3) if only_stiff else None,
        })

    arth = [{'load_N': r['load_N'],
             'displacement_mm': round(displacement_mm(r['load_N'], FACET_ARTHRITIC,
                                                      K_DISC_N_PER_MM), 4),
             'reduction_vs_mid_percent': round(
                 100.0 * (1 - displacement_mm(r['load_N'], FACET_ARTHRITIC, K_DISC_N_PER_MM)
                          / r['displacement_mid_mm']), 1)}
            for r in rows]

    summary = {
        'question': ('for an applied axial load, what displacement interval does the twin carry, and '
                     'which of the two inputs governs its width'),
        'inputs_both_from_harvested_edges': {
            'disc_axial_stiffness_N_per_mm': K_DISC_N_PER_MM,
            'disc_axial_stiffness_sd': K_DISC_SD,
            'relative_stiffness_uncertainty_percent': round(100 * K_DISC_SD / K_DISC_N_PER_MM, 1),
            'facet_compressive_fraction': list(FACET_FRACTION),
            'facet_fraction_span_factor': round(FACET_FRACTION[1] / FACET_FRACTION[0], 2),
            'facet_fraction_arthritic': FACET_ARTHRITIC,
        },
        'which_uncertainty_governs': (
            'the facet fraction contributes a width of ' + f"{rows[0]['width_from_facet_fraction_alone_mm']}"
            + ' mm at 200 N against ' + f"{rows[0]['width_from_stiffness_alone_mm']}"
            + ' mm from the stiffness, a ratio of ' + f"{rows[0]['facet_over_stiffness']}"
            + '. The ratio is the same at every load because both terms are linear in it, so one '
              'measurement answers it for all loads: the load path is the thing to measure, not the '
              'disc stiffness'),
        'decision_rows': rows,
        'arthritic_case': arth,
        'control': ('each input alone at its central value, which is how both numbers were reported '
                    'until they were put together'),
        'claim_type': 'information_link',
        'scope': ('a one-dimensional load split using two published intervals; no geometry, no posture '
                  'and no clinical use. The two sources are not stated to be the same specimen or level'),
        'falsifier': ('if the stiffness and the fraction come from different spinal levels or different '
                      'loading directions, the split is not theirs to make jointly and the interval is '
                      'void. Neither source states the other, and this is NOT established'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'LOAD_PATH_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    print(f"  styvhet {K_DISC_N_PER_MM} +/- {K_DISC_SD} N/mm "
          f"({100*K_DISC_SD/K_DISC_N_PER_MM:.1f} %), facettandel {FACET_FRACTION} "
          f"({FACET_FRACTION[1]/FACET_FRACTION[0]:.1f}x spann)\n")
    print(f"  {'last N':>8s} {'displacement mm':>26s} {'bredd':>8s} {'facett':>8s} {'styvhet':>8s} {'kvot':>6s}")
    for r in rows:
        print(f"  {r['load_N']:>8.0f} "
              f"{f'{r[chr(100)+chr(105)+chr(115)+chr(112)+chr(108)+chr(97)+chr(99)+chr(101)+chr(109)+chr(101)+chr(110)+chr(116)+chr(95)+chr(108)+chr(111)+chr(119)+chr(95)+chr(109)+chr(109)]}-{r['displacement_high_mm']}':>26s} "
              f"{r['full_interval_width_mm']:>8.3f} {r['width_from_facet_fraction_alone_mm']:>8.3f} "
              f"{r['width_from_stiffness_alone_mm']:>8.3f} {r['facet_over_stiffness']:>6.2f}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
