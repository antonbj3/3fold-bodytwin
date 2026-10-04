#!/usr/bin/env python3
"""Do a published source's own numbers agree with each other, before we score anything against them?

The question arrived by accident and then turned out to be general. The laser chain is scored against
a mass-loss law whose slope is 267 ug/J. The same paper separately reports a heat of ablation of
3740 J/g. Those are the same measurement read two ways, so one must be the reciprocal of the other,
and 1/3740 g/J is 267.38 ug/J -- they agree to 0.14 percent. Nothing in the chain required that check
and no report had made it.

Once asked, the question had an answer in five of our external facits, because a source that reports
two quantities linked by an identity has already tested itself. The identities are elementary and that
is the point: an elementary check that nobody runs is exactly where a bad number survives.

  reciprocal         slope in ug/J against heat of ablation in J/g
  isotropic elastic  shear modulus against Young modulus at a stated Poisson ratio
  dispersion         standard error against standard deviation over root n
  reciprocal         fibre coefficient k against its own denominator g
  product            pressure times area against the stated load

What it does NOT do. It does not reject a source whose numbers disagree -- a disagreement can mean the
two quantities were measured on different specimens, or that a convention differs, and those are
findings rather than errors. It prints the residual and names the identity, and the judgement stays
with the reader. It also makes no claim about sources not listed here: five instances is what has been
verified by hand, not a survey.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

W = Path(__file__).resolve().parents[2]
OUT = W / 'results/ASSEMBLY_FACIT_INTERNAL_CHECK'

CHECKS = [
    {'source': 'Payne BP et al., Lasers Surg Med 1998;23(1):1, PMID 9694144',
     'chain': 'laser ablation',
     'identity': 'slope [ug/J] = 1e6 / heat_of_ablation [J/g]',
     'reported': {'slope_ug_per_J': 267.0, 'heat_of_ablation_J_per_g': 3740.0},
     'derived': lambda r: 1e6 / r['heat_of_ablation_J_per_g'],
     'against': 'slope_ug_per_J'},
    {'source': 'Han & Eriten 2018, doi 10.1098/rsos.172051, with the normal modulus from doi 10.6084/m9.figshare.6231002.v2',
     'chain': 'instrument tissue adhesion',
     'identity': 'G = E / (2 (1 + nu)) at nu = 0.25',
     'reported': {'E_Pa': 32800000.0, 'G_held_Pa': 13280000.0, 'nu': 0.25},
     'derived': lambda r: r['E_Pa'] / (2 * (1 + r['nu'])),
     'against': 'G_held_Pa'},
    {'source': 'the mPTP rabbit preparation, PMID 15642769 and 15698843',
     'chain': 'mitochondrial pore ceiling',
     'identity': 'true SD = reported SEM * sqrt(n), at n = 7',
     'reported': {'reported_pm_pp': 4.0, 'n': 7, 'true_SD_if_SEM_pp': 10.583005244258363},
     'derived': lambda r: r['reported_pm_pp'] * math.sqrt(r['n']),
     'against': 'true_SD_if_SEM_pp'},
    # WITHDRAWN 2026-10-05. This row was never a source self-check and should not have been counted
    # as one. k was DEFINED as 1/g in the same computation, so checking k against 1/g tests the
    # arithmetic and nothing about any source -- compare the laser row, where the slope and the heat
    # of ablation are two separate measurements reported in one paper. Worse, the citation itself is
    # unverified: BT-NET-NUCLEUS-PRESSURE-MULTIPLIER-VS-ROUTE-OF-DERIVATION searched for a published
    # fibre equation yielding 4.877250096040077 and could not attach it to any paper, with Schroeder
    # et al. 2008 (PMID 18327799) and 2006 (PMID 16724211) the nearest candidates and neither
    # reporting that value. The name 'Ngwa' appears nowhere in this repository except in files this
    # coordinator wrote today, which means it entered from a model's output and was propagated as if
    # it were a reference. An untraceable number is an assertion, not a route.
    {'source': 'the Leeds 500 N knee model, doi 10.5518/981',
     'chain': 'meniscus contact',
     'identity': 'sum over condyles of mean pressure times nodal area = applied load',
     'reported': {'A1_mm2': 199.01469891419063, 'p1_MPa': 1.2885913888072562,
                  'A2_mm2': 252.57164890113802, 'p2_MPa': 1.0695118722455845,
                  'applied_load_N': 500.0},
     'derived': lambda r: r['A1_mm2'] * r['p1_MPa'] + r['A2_mm2'] * r['p2_MPa'],
     'against': 'applied_load_N'},
]


def main() -> int:
    rows = []
    for c in CHECKS:
        got = c['derived'](c['reported'])
        ref = c['reported'][c['against']]
        rel = abs(got - ref) / abs(ref) if ref else None
        rows.append({'chain': c['chain'], 'source': c['source'], 'identity': c['identity'],
                     'derived_value': got, 'reported_value': ref,
                     'relative_residual': rel,
                     'reported_quantities': {k: v for k, v in c['reported'].items()}})
    rels = [r['relative_residual'] for r in rows]
    summary = {
        'question': 'do a published source\'s own numbers satisfy the identity that links them',
        'claim_type': 'information_link',
        'why_it_matters': ('a source that reports two quantities linked by an identity has already '
                           'tested itself, and an elementary check nobody runs is where a bad number '
                           'survives'),
        'checks': len(rows),
        'relative_residual': {'min': min(rels), 'max': max(rels),
                              'median': sorted(rels)[len(rels) // 2]},
        'all_below_5_percent': all(x < 0.05 for x in rels),
        'what_a_disagreement_means': ('different specimens, or a different convention; both are '
                                      'findings rather than errors, so this file prints the residual '
                                      'and names the identity and leaves the judgement to the reader'),
        'scope': ('four instances verified by hand, not a survey of all sources in the net. A fifth '
                  'was withdrawn on 2026-10-05: it compared a number against its own definition, '
                  'which tests arithmetic and not a source, and its citation could not be traced to '
                  'any paper'),
        'rows': rows,
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'CHECK_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f'{"kedja":30} {"identitet":46} {"residual":>10}')
    for r in rows:
        print(f"  {r['chain']:28} {r['identity'][:44]:44} {100 * r['relative_residual']:9.3f} %")
    print(f"  {len(rows)} kontroller: min {100 * min(rels):.3f} %, median "
          f"{100 * sorted(rels)[len(rels) // 2]:.3f} %, max {100 * max(rels):.3f} %")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
