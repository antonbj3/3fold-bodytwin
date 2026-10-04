#!/usr/bin/env python3
"""How many extra fields does each chain need before its summary can carry its decision?

Six times on 2026-10-04 and 05 a summary with an identity error of exactly zero hid states that the
decision separates. Each case was found by hand in a different lane, and each was logged as a
surprise. Six surprises in six independent chains is a property of how we summarise, not six
accidents, so it gets a ledger.

PROOF_LANE_MINIMAL_FIELD supplied the arithmetic that makes the ledger more than a list. Let M be the
largest number of distinct decisions inside ONE summary fibre -- one set of states the summary cannot
tell apart. With full-state access the minimum number of extra fields of q symbols each is
ceil(log_q M), and that bound is sharp. Two further results from the same run shape what this file
may and may not do:

  - every f = h(sigma) preserves an existing witness, so no reprocessing of what the summary already
    carries can repair it; the ledger therefore never proposes a derived field as a fix
  - the COARSEST sufficient joint partition is unique, while minimum-SIZED fields need not be, so the
    ledger reports a count and refuses to name the field

What the ledger does not do, deliberately. It does not compute M from first principles -- general
program-level minimisation is undecidable, as the same run established. It reads M off the witness
pairs each lane actually measured, which means M is a LOWER bound on the true M: a chain where two
states were found to differ has M >= 2, and finding a third separated state would raise it. Every row
therefore carries the count of witnesses it rests on, and a row resting on one pair says so.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

W = Path(__file__).resolve().parents[2]
OUT = W / 'results/ASSEMBLY_SUFFICIENCY_LEDGER'

# One row per chain where a zero identity error was measured alongside a decision-relevant gap.
# 'witnesses' is the number of distinct states the lane actually exhibited inside one fibre.
CHAINS = [
    {'chain': 'instrument tissue adhesion', 'lane': 'LANE_TISSUE_ADHESION_CONTACT r4',
     'summary_identity_error': 0.0, 'also_identical': 'normal field, identity error 0.0',
     'gap': 10.0, 'gap_unit': 'ratio of peak force', 'witnesses': 2,
     'decision_scale': 'the peak force the decision consumes',
     'gap_over_scale': 10.0,
     'note': 'the largest ratio measured; summary AND normal field both identical'},
    {'chain': 'complement window', 'lane': 'LANE_IMMUNE_OVERACTIVATION r1',
     'summary_identity_error': 0.0, 'also_identical': 'both identity errors 0.0',
     'gap': 17.60797843701942, 'gap_unit': 'nM MAC', 'witnesses': 2,
     'decision_scale': 'primary numerical interval 0 to 0.5441770553588867 nM',
     'gap_over_scale': 17.60797843701942 / 0.5441770553588867,
     'note': 'gap exceeds the whole interval the decision operates in'},
    {'chain': 'disc electrical alias', 'lane': 'LANE_DISC_LOAD_LABEL r5',
     'summary_identity_error': 0.0, 'also_identical': 'electrical summary exact',
     'gap': 0.08108314647592979, 'gap_unit': 'MPa', 'witnesses': 2,
     'decision_scale': 'facit band width 0.65 - 0.53 = 0.12 MPa',
     'gap_over_scale': 0.08108314647592979 / 0.12,
     'note': 'the only one so far where the alias CHANGES the diagnosis'},
    {'chain': 'delivered electrical power', 'lane': 'LANE_SETTING_TO_DIFFUSIVITY r3',
     'summary_identity_error': 0.0, 'also_identical': 'mean V and mean I both identical',
     'gap': 20.0, 'gap_unit': 'W', 'witnesses': 2,
     'decision_scale': 'product of means, about 40 W',
     'gap_over_scale': 20.0 / 40.0,
     'note': ('witness pair is 50 against 30 W; an earlier coordinator row compared the 20 W gap '
              'against an unrelated 5 W terminal case and reported 4.0x, which was a conflation')},
    {'chain': 'biofilm plateau', 'lane': 'LANE_BIOFILM_RETENTION_EXPONENT r7',
     'summary_identity_error': 0.0, 'also_identical': 'every observed point',
     'gap': 0.1875, 'gap_unit': 'normalised signal', 'witnesses': 2,
     'decision_scale': 'the decision interval itself',
     'gap_over_scale': 1.0,
     'note': 'plateaus 4/5 against 1 give 15/16 against 3/4, gap exactly 3/16'},
    {'chain': 'distribution gate initial field', 'lane': 'LANE_DISTRIBUTION_GATE_43 r9',
     'summary_identity_error': 0.0, 'also_identical': 'initial total field exact',
     'gap': 0.148502787956807, 'gap_unit': 'field units', 'witnesses': 2,
     'decision_scale': None, 'gap_over_scale': None,
     'note': 'scale not stated by the lane, so no ratio is formed here'},
]
# Categorical result from PROOF_LANE_SCALE_BRIDGE, kept separate because it is not a witness pair.
MEAN_VARIANCE_FAILED = {'pairs_tested': 378, 'pairs_failed': 378,
                        'reading': 'mean plus variance is never a sufficient summary in that family'}


def min_fields(M: int, q: int = 2) -> int:
    """PROOF_LANE_MINIMAL_FIELD: with full-state access the minimum is ceil(log_q M), and it is sharp."""
    return 0 if M <= 1 else math.ceil(math.log(M, q))


def main() -> int:
    rows = []
    for c in CHAINS:
        M = c['witnesses']
        rows.append({**c, 'M_lower_bound': M,
                     'min_extra_binary_fields': min_fields(M, 2),
                     'min_extra_fields_is_a_lower_bound_because':
                         'M is read off the witnesses the lane exhibited, not derived; a third '
                         'separated state would raise it'})
    total = sum(r['min_extra_binary_fields'] for r in rows)
    ratios = [r['gap_over_scale'] for r in rows if r['gap_over_scale'] is not None]
    summary = {
        'question': 'per chain, how many extra fields must the summary carry before it can decide',
        'claim_type': 'information_link',
        'bound': 'ceil(log_q M) fields of q symbols, from PROOF_LANE_MINIMAL_FIELD, sharp with full-state access',
        'why_no_field_is_named': ('the coarsest sufficient joint partition is unique but '
                                  'minimum-sized fields need not be, so a named field would be one '
                                  'choice among several presented as the choice'),
        'why_no_derived_field_is_proposed': ('every f = h(sigma) preserves an existing witness, so '
                                            'reprocessing what the summary already carries cannot repair it'),
        'chains': len(rows),
        'minimum_extra_binary_fields_over_all_chains': total,
        'gap_over_decision_scale': {'min': min(ratios), 'max': max(ratios),
                                    'median': sorted(ratios)[len(ratios) // 2]},
        'mean_plus_variance_categorical_failure': MEAN_VARIANCE_FAILED,
        'control': ('carry the raw state, which is always sufficient and never compresses -- what the '
                    'repo does today when it gives up'),
        'falsifier': ('a chain whose summary separates its witnesses after one extra field, when this '
                      'ledger says it needs more, falsifies the M read off that lane; the lane then '
                      'has to exhibit the third state or withdraw the row'),
        'rows': rows,
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
        'no_claim_of_biological_validation': True,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'LEDGER_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f'{"kedja":34} {"vittnen":>8} {"M>=":>5} {"min falt":>9} {"gap/skala":>11}')
    for r in rows:
        g = f"{r['gap_over_scale']:.4f}" if r['gap_over_scale'] is not None else '-'
        print(f"  {r['chain']:32} {r['witnesses']:8} {r['M_lower_bound']:5} "
              f"{r['min_extra_binary_fields']:9} {g:>11}")
    print(f'  sum minimum extra binary fields above {len(rows)} kedjor: {total}')
    print(f"  gap mot beslutets skala: min {min(ratios):.4f}, median "
          f"{sorted(ratios)[len(ratios) // 2]:.4f}, max {max(ratios):.4f}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
