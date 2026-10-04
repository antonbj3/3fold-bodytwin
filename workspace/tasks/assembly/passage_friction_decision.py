#!/usr/bin/env python3
"""A decision outside the eye: what forward friction should a needle plan for on a repeated passage?

Why it can be built now. Until today the net assumed friction is equal across passages, with the
stated reason that no digitised series separates them. That reason fell twice and the replacement
claim -- that a first-pass number was missing -- fell too. The number is printed in the source's own
figure under the wrong name: past the penetration mark no material lies ahead of the tip, so the
cutting term is identically zero on BOTH passages and the measured force is pure shaft-tissue
friction. Kataoka, Washio, Chinzei, Mizuhara, Simone and Okamura, MICCAI 2002, LNCS 2488:216-223:

    first passage   0.349 +/- 0.018 N
    second passage  0.241 +/- 0.004 N

This coordinator recomputed the ratio and propagated the dispersion rather than taking it:
1.4481 with a standard deviation of 0.0785, which is 5.71 sigma from unity. Friction is not equal
across passages, and now that is a measurement rather than an assumption.

What the decision answers. Given a planned pass count and a friction budget -- the largest forward
force the instrument or the operator will accept -- it returns the number of passes that fit inside
the budget, and refuses when the budget cannot even accommodate the first pass.

The model is deliberately thin, because the data supports exactly two points. Pass 1 carries
F1; every later pass in the same hole carries F2. Nothing here interpolates a decay curve, because
there is no third passage in the source to constrain one, and a fitted decay would be the kind of
invention that makes a decision look stronger than its evidence.

The control is the equal-friction assumption the edge used to carry: every pass at F2. That is the
equally informed comparison, since it uses the same source and the same budget, and the difference
between the two is the whole value of the finding.

The falsifier is depth. The source's own second-passage fit is linear in the product of load and
shaft contact area, F_2nd = 10.679/N * (P*S_t) - 0.01333 N over depths 10, 15 and 20 mm, so the
friction at a given pass is depth dependent and the two numbers above are one depth. If the ratio
moves with depth by more than its own 0.0785 dispersion, a single ratio cannot carry the decision
and this file must take depth as an input instead.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

W = Path('.')
OUT = W / 'results/ASSEMBLY_PASSAGE_FRICTION'

F_FIRST_N, F_FIRST_SD_N = 0.349, 0.018
F_SECOND_N, F_SECOND_SD_N = 0.241, 0.004
SOURCE = 'Kataoka et al., MICCAI 2002, LNCS 2488:216-223, Fig. 6 post-penetration segment'


def ratio_with_dispersion() -> tuple[float, float]:
    r = F_FIRST_N / F_SECOND_N
    sd = r * math.hypot(F_FIRST_SD_N / F_FIRST_N, F_SECOND_SD_N / F_SECOND_N)
    return r, sd


def total_force_N(passes: int, equal_friction: bool) -> float:
    """Cumulative forward friction over a pass sequence in one hole."""
    if passes <= 0:
        return 0.0
    if equal_friction:
        return passes * F_SECOND_N
    return F_FIRST_N + (passes - 1) * F_SECOND_N


def decide(budget_N: float, planned_passes: int) -> dict:
    out = {'budget_N': budget_N, 'planned_passes': planned_passes}
    if budget_N < F_FIRST_N - F_FIRST_SD_N:
        out['decision'] = 'NO_PASS_COUNT_RETURNED'
        out['reason'] = (f'budget {budget_N} N is below the first passage even at its lower '
                         f'dispersion bound {F_FIRST_N - F_FIRST_SD_N:.3f} N; no pass count fits')
        return out
    fits = max(n for n in range(0, planned_passes + 1) if total_force_N(n, False) <= budget_N)
    fits_control = max(n for n in range(0, planned_passes + 1)
                       if total_force_N(n, True) <= budget_N)
    out.update({
        'decision': 'PASSES_WITHIN_BUDGET',
        'decision_value_passes': fits,
        'control_equal_friction_passes': fits_control,
        'passes_the_control_would_have_allowed_in_excess': fits_control - fits,
        'total_force_at_decision_N': total_force_N(fits, False),
        'total_force_control_N': total_force_N(fits_control, True),
        'decision_changes_vs_control': fits != fits_control,
    })
    return out


def main() -> int:
    r, sd = ratio_with_dispersion()
    # Six hand-picked budgets gave the same pass count as the control in all six, which looked like
    # the finding buying nothing. It was a sampling artefact: the two readings differ by the constant
    # F_FIRST - F_SECOND = 0.108 N, so the control over-allows exactly when the budget falls in a
    # 0.108 N window below a multiple of 0.241 N. The budget axis is therefore scanned instead of
    # sampled, and the fraction of it where the decision changes is reported.
    rows = [decide(b, 8) for b in (0.30, 0.35, 0.48, 0.60, 0.72, 0.85, 0.96, 1.10, 1.60)]
    scan = [decide(round(0.20 + 0.001 * i, 3), 12) for i in range(1801)]
    scan_live = [x for x in scan if x['decision'] == 'PASSES_WITHIN_BUDGET']
    scan_changed = [x for x in scan_live if x['decision_changes_vs_control']]
    changed = [x for x in rows if x.get('decision_changes_vs_control')]
    summary = {
        'question': 'how many passes in one hole fit inside a forward-friction budget',
        'claim_type': 'capability',
        'source': SOURCE,
        'first_passage_N': [F_FIRST_N, F_FIRST_SD_N],
        'second_passage_N': [F_SECOND_N, F_SECOND_SD_N],
        'passage_friction_ratio': r,
        'passage_friction_ratio_sd': sd,
        'sigma_from_unity': (r - 1.0) / sd,
        'control': 'the equal-friction assumption the edge used to carry: every pass at the second-passage force',
        'rows': rows,
        'budgets_where_the_decision_changes': len(changed),
        'budgets_tested': len(rows),
        'budget_axis_scan': {
            'range_N': [0.201, 2.0], 'step_N': 0.001,
            'budgets_with_a_decision': len(scan_live),
            'budgets_where_the_decision_differs_from_the_control': len(scan_changed),
            'fraction_of_the_budget_axis': len(scan_changed) / len(scan_live) if scan_live else None,
            'why': ('the two readings differ by the constant F_first - F_second = '
                    f'{F_FIRST_N - F_SECOND_N:.3f} N, so the control over-allows one pass whenever '
                    f'the budget falls in that window below a multiple of {F_SECOND_N} N'),
            'closed_form_fraction': (F_FIRST_N - F_SECOND_N) / F_SECOND_N,
            'closed_form_note': ('(F_first - F_second)/F_second is the fraction of each '
                                 'F_second-wide budget period in which the control over-allows, '
                                 'and it confirms the scan independently'),
            'reading': ('the equal-friction assumption lets through one pass too many on this '
                        'fraction of the budget axis; on the rest the two agree, which is why six '
                        'hand-picked budgets showed no difference at all')},
        'falsifier': ('the source fit F_2nd = 10.679/N * (P*S_t) - 0.01333 N is linear over depths '
                      '10, 15 and 20 mm, so friction is depth dependent and these two numbers are '
                      'one depth; if the ratio moves with depth by more than its own 0.0785 '
                      'dispersion, a single ratio cannot carry this decision'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
        'no_claim_of_biological_validation': True,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'DECISION_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f'  passagefriktionskvot {r:.4f} +/- {sd:.4f}, {(r - 1) / sd:.2f} sigma fran ett')
    print(f'{"budget N":>9} {"beslut":>8} {"kontroll":>9} {"overskott":>10} {"kraft N":>9}')
    for x in rows:
        if x['decision'] == 'NO_PASS_COUNT_RETURNED':
            print(f'{x["budget_N"]:9} {"INGET":>8}')
        else:
            print(f'{x["budget_N"]:9} {x["decision_value_passes"]:8} '
                  f'{x["control_equal_friction_passes"]:9} '
                  f'{x["passes_the_control_would_have_allowed_in_excess"]:10} '
                  f'{x["total_force_at_decision_N"]:9.3f}')
    print(f'  beslutet skiljer sig fran kontrollen i {len(changed)} av {len(rows)} valda budgetar')
    print(f'  pa hela budgetaxeln 0,201-2,000 N i steg om 1 mN: '
          f'{len(scan_changed)} av {len(scan_live)} = '
          f'{100 * len(scan_changed) / len(scan_live):.1f} % av budgetarna slapper kontrollen '
          f'igenom en passage for mycket')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
