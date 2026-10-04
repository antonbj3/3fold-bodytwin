"""Regulatory failure as a two-sided window: spare the host without losing the ability to discriminate.

This is the seed stated as under- and over-activation of the immune system expressed as regulatory
failure, and it is the first thing built against it. The two sides are already measured in our own net
and neither was ever put against the other:

  host sparing   the self-surface growth rate crosses zero at a recycle fraction of 0.6901150139456149,
                 and below it the host is spared while the activator is lysed throughout
  discrimination activator-over-host discrimination is 486.11 without recycling and 4.969 with it, a
                 factor of 97.83, and published biology supports the recycling end

Read separately, each says recycling is good for one thing. Read together they are a trade-off: the
recycling that spares the host is the same recycling that destroys the ability to tell an activator from
a host surface. Over-activation and under-activation are therefore not two failures of one knob in
opposite directions -- they are one knob with a different failure at each end, and the question is
whether any setting satisfies both.

The decision: given a required discrimination factor, what recycle fractions are admissible, and does
the admissible set intersect the host-sparing range at all? The honest answer may be that it does not,
and that would be the finding: no setting of this parameter is regulatorily sound, which means the
mechanism needs a second control the model does not have.

Everything is computed from the supplied cell, not interpolated between the two endpoint numbers.
The control is each constraint alone, which is how both were reported until now.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

W = Path('')
CELL = W / 'tasks/free48/sources/COMPLEMENT_DISCRIMINATION/cell.py'
OUT = W / 'results/ASSEMBLY_IMMUNE_WINDOW'
HOST_FLIP = 0.6901150139456149

spec = importlib.util.spec_from_file_location('complement', CELL)
c = importlib.util.module_from_spec(spec)
sys.path.insert(0, str(CELL.parent))
spec.loader.exec_module(c)


def growth(surface: str, rho: float) -> float:
    p = c.params(surface=surface)
    p['decay_recycle_fraction'] = rho
    return float(c.growth(p)['lambda_per_min'])


def discrimination(rho: float) -> float:
    """Activator growth over host growth at the same recycle fraction, as a ratio of rates."""
    a, h = growth('activator', rho), growth('host', rho)
    if abs(h) < 1e-12:
        return float('inf')
    return a / h


def main() -> int:
    rows = []
    n = 1001
    for i in range(n):
        rho = i / (n - 1)
        a, h = growth('activator', rho), growth('host', rho)
        rows.append({'rho': rho, 'activator': a, 'host': h,
                     'host_spared': h < 0, 'discrimination': discrimination(rho)})

    spared = [r for r in rows if r['host_spared']]
    # the admissible set for a required discrimination, intersected with host sparing
    requirements = [2.0, 5.0, 10.0, 50.0, 100.0, 486.0]
    verdicts = []
    for req in requirements:
        ok = [r for r in rows if r['host_spared'] and abs(r['discrimination']) >= req]
        verdicts.append({
            'required_discrimination': req,
            'admissible_rho_count_of_1001': len(ok),
            'rho_low': round(min((r['rho'] for r in ok), default=float('nan')), 6) if ok else None,
            'rho_high': round(max((r['rho'] for r in ok), default=float('nan')), 6) if ok else None,
            'window_exists': bool(ok),
        })

    summary = {
        'question': ('is there a recycle fraction that spares the host and still discriminates an '
                     'activator from a host surface'),
        'the_two_sides_each_measured_separately_until_now': {
            'host_sparing_flip': HOST_FLIP,
            'discrimination_without_recycling': 486.11,
            'discrimination_with_recycling': 4.969,
            'ratio': 97.83,
            'reading': ('the recycling that spares the host is the same recycling that destroys '
                        'discrimination, so over- and under-activation are one knob with a different '
                        'failure at each end rather than one failure in two directions'),
        },
        'host_spared_range': {
            'rho_low': round(min(r['rho'] for r in spared), 6) if spared else None,
            'rho_high': round(max(r['rho'] for r in spared), 6) if spared else None,
            'count_of_1001': len(spared),
        },
        'discrimination_at_the_ends': {
            'rho_0': round(discrimination(0.0), 4),
            'rho_at_flip_minus': round(discrimination(HOST_FLIP - 1e-3), 4),
            'rho_1': round(discrimination(1.0), 4),
        },
        'verdicts': verdicts,
        'control': ('each constraint alone, which is how both numbers were reported until now: host '
                    'sparing said recycling below the flip is safe, discrimination said recycling '
                    'costs a factor of 97.83, and neither was put against the other'),
        'claim_type': 'capability',
        'scope': ('algebraic behaviour of one supplied model under a swept parameter; the branch label '
                  'the cell does not carry still conditions every number here'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'IMMUNE_WINDOW_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    hr = summary['host_spared_range']
    print(f"  the host is spared for rho in [{hr['rho_low']}, {hr['rho_high']}], "
          f"{hr['count_of_1001']} of 1001 points")
    d = summary['discrimination_at_the_ends']
    print(f"  discrimination: rho=0 -> {d['rho_0']}, just below the turning point -> "
          f"{d['rho_at_flip_minus']}, rho=1 -> {d['rho_1']}")
    print(f"\n  {'requirement':>8s} {'admissible rho of 1001':>22s} {'interval':>26s}")
    for v in verdicts:
        rng = f"[{v['rho_low']}, {v['rho_high']}]" if v['window_exists'] else 'EMPTY'
        print(f"  {v['required_discrimination']:>8.0f} {v['admissible_rho_count_of_1001']:>22d} {rng:>26s}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
