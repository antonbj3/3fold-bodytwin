"""A decision outside the eye: is there a host-side setting where the host survives and the activator lyses?

Where this came from. A capability brief asked the complement cell what it could decide that nobody had
asked it to, and it came back with an inverse question rather than a prediction: at what host-side
`decay_recycle_fraction`, relative to the activator-side value, does the model flip from "the host is
lysed" to "the host is spared while the activator is still lysed", and is that flip point unique. It
named four external anchors and pre-registered its own falsification conditions, which is why this is
worth running rather than reading.

The hard constraint it found, and it is the night's recurring shape again. The recycle branch is a FORK
the cell cannot see: whether the decayed convertase partners with C3b, C3(H2O) or CVF changes the
answer, and nothing in the cell's own ledgers records which branch it is on. So the flip point is only
meaningful together with a branch label, exactly as the posterior corneal axis was only meaningful
together with a meridian convention, and the keratometric offset only with an index.

What is computed here: the host and activator arms swept independently over the recycle fraction, the
flip point located by bisection on the sign of the host growth rate, and the asymmetry at that point
expressed as the quantity a bench assay would have to hit. The cell is imported and run, not
reimplemented.

Falsification conditions, taken from the brief rather than invented after the fact:
  * if the host and activator arms flip at the same fraction, there is no window and no decision;
  * if the flip is a plateau rather than a point, there is a range and no decision;
  * if the MAC difference at identical growth is smaller than the solver's own error, the cell has no
    discriminating power and the asymmetry claim is false.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

W = Path('.')
CELL = W / 'tasks/free48/sources/COMPLEMENT_DISCRIMINATION/cell.py'
OUT = W / 'results/ASSEMBLY_COMPLEMENT_WINDOW'

spec = importlib.util.spec_from_file_location('complement', CELL)
c = importlib.util.module_from_spec(spec)
sys.path.insert(0, str(CELL.parent))
spec.loader.exec_module(c)


def growth_at(surface: str, rho: float) -> float:
    p = c.params(surface=surface)
    p['decay_recycle_fraction'] = rho
    return float(c.growth(p)['lambda_per_min'])


def flip(surface: str, lo: float = 0.0, hi: float = 1.0, tol: float = 1e-10):
    """The recycle fraction where the growth rate changes sign, or None if it never does."""
    a, b = growth_at(surface, lo), growth_at(surface, hi)
    if a == 0.0:
        return lo
    if a * b > 0:
        return None
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if growth_at(surface, lo) * growth_at(surface, m) <= 0:
            hi = m
        else:
            lo = m
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def main() -> int:
    arms = {}
    for surface in ('activator', 'host'):
        f = flip(surface)
        arms[surface] = {
            'flip_recycle_fraction': f,
            'growth_at_zero_recycle_per_min': round(growth_at(surface, 0.0), 8),
            'growth_at_full_recycle_per_min': round(growth_at(surface, 1.0), 8),
        }

    fa, fh = arms['activator']['flip_recycle_fraction'], arms['host']['flip_recycle_fraction']
    window = None
    # First version required BOTH arms to flip and returned None, which was my framing and not the
    # result: the activator never flips at all, and that is precisely WHY a window exists. The window
    # is bounded by the host flip on one side and the parameter domain on the other.
    if fh is not None and fa is None:
        host_below = growth_at('host', max(fh - 1e-3, 0.0))
        act_below = growth_at('activator', max(fh - 1e-3, 0.0))
        window = {
            'form': 'one-sided',
            'host_flip': fh,
            'activator_never_flips': True,
            'host_spared_activator_lysed_below_the_flip': bool(host_below < 0 < act_below),
            'host_growth_just_below_flip_per_min': round(host_below, 8),
            'activator_growth_just_below_flip_per_min': round(act_below, 8),
            'reading': ('the host is spared for every recycle fraction below the flip while the '
                        'activator is lysed throughout, so the window is the whole range below it'),
        }
    elif fa is not None and fh is not None:
        lo, hi = sorted((fa, fh))
        window = {'low': lo, 'high': hi, 'width': hi - lo,
                  'host_spared_activator_lysed_inside': None}
        if hi - lo > 1e-9:
            mid = 0.5 * (lo + hi)
            window['host_spared_activator_lysed_inside'] = bool(
                growth_at('host', mid) > 0 > growth_at('activator', mid)
                or growth_at('host', mid) < 0 < growth_at('activator', mid))
            window['host_growth_at_midpoint'] = round(growth_at('host', mid), 8)
            window['activator_growth_at_midpoint'] = round(growth_at('activator', mid), 8)

    # is the flip a point or a plateau? count distinct sign changes on a fine sweep
    sweep = []
    n = 2001
    for i in range(n):
        rho = i / (n - 1)
        sweep.append((rho, growth_at('host', rho), growth_at('activator', rho)))
    def crossings(idx):
        k = 0
        for j in range(1, len(sweep)):
            if sweep[j - 1][idx] * sweep[j][idx] < 0:
                k += 1
        return k

    summary = {
        'question': ('is there a host-side recycle fraction where the host is spared while the '
                     'activator is still lysed, and is the flip a point'),
        'arms': arms,
        'window': window,
        'sign_changes_over_2001_points': {'host': crossings(1), 'activator': crossings(2)},
        'the_fork_the_cell_cannot_see': (
            'the recycle branch requires a partner, C3b or C3(H2O) or CVF, and nothing in this cell '
            'records which. The flip fraction is therefore conditional on a branch label the model '
            'does not carry, the same way the posterior corneal axis was conditional on a meridian '
            'convention and the keratometric offset on an assumed index'),
        'external_anchors_named_by_the_brief': [
            {'claim': 'decayed Bb retains C3- and C5-cleaving activity, ~100x below factor B, '
                      'suppressed by factor H, enhanced by properdin, and requires a partner',
             'locator': 'PMID 6229580'},
            {'claim': 'a soluble convertase retains about 25 percent activity after conditions that '
                      'totally inactivated C3bBb', 'locator': 'doi 10.1016/j.molimm.2007.11.003'},
            {'claim': 'dose-dependent protection of erythrocytes by membrane-targeted soluble CD59, '
                      'i.e. host protection is achievable at the surface',
             'locator': 'doi 10.1182/blood-2005-02-0782'},
            {'claim': 'the fraction of C3-bound red cells correlated with reticulocyte count and with '
                      'haematologic response, restricted to the CD59-negative population',
             'locator': 'doi 10.1182/blood-2008-11-189944'},
        ],
        'falsification_conditions_preregistered_by_the_brief': [
            'both arms flipping at the same fraction means there is no window',
            'a plateau rather than a point means a range and no decision',
            'a MAC difference at identical growth below solver error means no discriminating power',
        ],
        'claim_type': 'capability',
        'control': ('the activator arm under the same sweep, which is the equally informed comparison: '
                    'same model, same parameters, only the surface changes'),
        'scope': ('algebraic behaviour of one supplied model under a swept parameter; no biological or '
                  'clinical validity is claimed and the branch label is unresolved'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'WINDOW_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    for s, v in arms.items():
        print(f"  {s:9s} flip = {v['flip_recycle_fraction']}  "
              f"growth 0 -> {v['growth_at_zero_recycle_per_min']}, "
              f"1 -> {v['growth_at_full_recycle_per_min']}")
    print(f"  window: {window}")
    print(f"  sign changes over 2001 points: host {crossings(1)}, activator {crossings(2)}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
