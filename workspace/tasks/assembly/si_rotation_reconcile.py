"""Four independent methods measured the same joint rotation. What interval should a decision carry?

Why this one. Every other reconciliation tonight had two sources at most, and the eye taught the hard
lesson twice: a factor between two methods is either a convention or a real disagreement, and you
cannot tell which from the factor alone. The sacroiliac document carries FOUR independent measurements
of the same rotation by three different methods plus one independent computation, each with a PMID, and
it labels which is which. That is enough to separate method from quantity instead of guessing.

The sources, as the document labels them and with its own cautions kept:

  Egund 1978, PMID 717034       ~2.0 deg max, n=4, first roentgen stereophotogrammetry, and the
                                document notes the figure is at the axis precision limit
  Sturesson 1989, PMID 2922636  2.5 deg mean, range 0.8 to 3.9, n=25, refined RSA in patients
  Jacob & Kissling 1995         1.7 deg, n=24, 3-D photogrammetry with EXTERNAL markers and the
                                document flags explicitly that this is not X-ray
  Sturesson 2000, PMID 10703111 RSA during a loaded functional task
  Zheng 1997, PMID 9136197      quasi-static finite element, and the document flags it as an
                                independent COMPUTATION and not a measurement

The document also pre-registered six gates and reports them all PASS, including a cross-method
over-determination of 1.68x and 1.81x against a 3x bar. So the cross-method factor is already known and
is NOT being discovered here; what is being decided is what a twin should carry given it.

The decision. Three candidate intervals, and the choice between them is the decision:
  all sources pooled            treats a computation and an external-marker method as measurements
  X-ray methods only            keeps the two that image bone directly
  the widest single-study range the one study that reported its own spread, n=25

The control is the pooled reading, which is what using every number in the document would give.
"""
from __future__ import annotations

import json
import statistics
from pathlib import Path

W = Path('')
OUT = W / 'results/ASSEMBLY_SI_ROTATION'
SRC = 'source_documents/SACROILIAC_JOINT.md'

SOURCES = [
    {'study': 'Egund 1978', 'pmid': '717034', 'method': 'RSA', 'images_bone': True,
     'is_measurement': True, 'n': 4, 'rotation_deg': 2.0,
     'caution': 'at the axis precision limit per the document'},
    {'study': 'Sturesson 1989', 'pmid': '2922636', 'method': 'RSA refined', 'images_bone': True,
     'is_measurement': True, 'n': 25, 'rotation_deg': 2.5, 'own_range_deg': [0.8, 3.9],
     'caution': None},
    {'study': 'Jacob & Kissling 1995', 'pmid': '11415579', 'method': '3-D photogrammetry',
     'images_bone': False, 'is_measurement': True, 'n': 24, 'rotation_deg': 1.7,
     'caution': 'external markers, explicitly NOT X-ray per the document'},
    {'study': 'Zheng 1997', 'pmid': '9136197', 'method': 'quasi-static FE', 'images_bone': False,
     'is_measurement': False, 'n': None, 'rotation_deg': None,
     'caution': 'independent computation, not a measurement, per the document'},
]
CROSS_METHOD_FACTOR = [1.68, 1.81]
CROSS_METHOD_BAR = 3.0


def interval(vals: list[float]) -> dict:
    return {'low': min(vals), 'high': max(vals), 'mean': round(statistics.mean(vals), 4),
            'spread_factor': round(max(vals) / min(vals), 3), 'n_sources': len(vals)}


def main() -> int:
    measured = [s for s in SOURCES if s['is_measurement'] and s['rotation_deg'] is not None]
    xray = [s for s in measured if s['images_bone']]

    pooled = interval([s['rotation_deg'] for s in measured])
    xray_only = interval([s['rotation_deg'] for s in xray])
    widest_single = next(s for s in SOURCES if s.get('own_range_deg'))

    candidates = {
        'all_measurements_pooled': {
            **pooled,
            'includes': [s['study'] for s in measured],
            'objection': ('includes a method the document flags as not X-ray, so a disagreement between '
                          'it and the others cannot be attributed to the joint'),
        },
        'x_ray_methods_only': {
            **xray_only,
            'includes': [s['study'] for s in xray],
            'objection': ('one of the two is flagged as sitting at its own axis precision limit, so the '
                          'narrow spread may be two readings of the same limitation'),
        },
        'widest_single_study_range': {
            'low': widest_single['own_range_deg'][0], 'high': widest_single['own_range_deg'][1],
            'mean': widest_single['rotation_deg'], 'n_subjects': widest_single['n'],
            'spread_factor': round(widest_single['own_range_deg'][1]
                                   / widest_single['own_range_deg'][0], 3),
            'includes': [widest_single['study']],
            'objection': ('one study, but the only one that reported its own between-subject spread, and '
                          'n=25 is the largest here'),
        },
    }

    recommendation = ('the widest single-study range, 0.8 to 3.9 degrees. It is the only interval in the '
                      'document that measures BETWEEN-SUBJECT variation rather than between-method '
                      'variation, and a decision about a subject needs the former. The pooled interval '
                      'is narrower at 1.7 to 2.5 and that narrowness is an artefact of comparing study '
                      'means to each other, which averages away exactly the variation a decision faces')

    summary = {
        'question': 'what sacroiliac rotation interval should a decision carry, given four methods',
        'source': {'path': SRC, 'labels_and_cautions_kept': True},
        'sources': SOURCES,
        'cross_method_factor_already_in_the_document': {
            'values': CROSS_METHOD_FACTOR, 'bar': CROSS_METHOD_BAR, 'gate': 'G3 PASS',
            'note': 'not discovered here; it is the premise this decision is made under',
        },
        'candidate_intervals': candidates,
        'recommendation': recommendation,
        'the_trap_this_avoids': ('a pooled mean of study means looks tighter than any single study and '
                                 'is tighter for the wrong reason: averaging across methods removes '
                                 'between-subject spread, which is the only spread a per-subject '
                                 'decision has to survive. Same shape as reading a composed standard '
                                 'deviation as a floor earlier tonight'),
        'control': 'the pooled reading, which is what using every number in the document would give',
        'claim_type': 'information_link',
        'scope': 'reconciliation of published values as labelled; no new measurement and no clinical use',
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'SI_ROTATION_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    for name, c in candidates.items():
        print(f"  {name:28s} [{c['low']}, {c['high']}] degrees, spread factor "
              f"{c['spread_factor']}, {len(c['includes'])} sources")
    print(f"\n  method separation already in the document: {CROSS_METHOD_FACTOR} against limit {CROSS_METHOD_BAR}")
    print(f"  recommendation: {recommendation[:96]}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
