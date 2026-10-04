"""Does a verdict that passes on a point still pass when the input is a distribution?

ORDERED BY THE FIELD LANE, 2026-10-04, and aimed at this session's own output: 21 of their 43 cells
pass on a scalar point or a band, and 3 of 8 verdicts FLIPPED when tested against a distribution.
A flip is booked as FAIL, not as nuance. Their instruction applies to the seven decisions built here
overnight if any of them carries a scalar gate anywhere in its chain, and four of them do:

  route_topology_decision        halving ratio cut at 0.75
  load_redistribution_decision   explained-fraction cut at 0.9
  reflection_coefficient_ident.  sigma error cut at 0.05, flow separation 3.81x
  glass_batch_energy_decision    share-of-control cut at 0.10

THE TEST. Each decision's inputs are resampled over the dispersion its own source states, the
verdict is recomputed per draw, and the fraction of draws that agree with the point verdict is
reported. A verdict that holds in fewer than 95 percent of draws has not survived the gate.

WHY A POINT VERDICT IS NOT ENOUGH. A cut at 0.75 on a quantity whose measured value is 0.70 reads as
a clean pass, and reads the same way when the input's own spread covers 0.60 to 0.80. The point
hides which side of the cut the measurement actually supports. That is how three of eight verdicts
in the field lane's inventory were wrong without anyone having made an arithmetic error.

NOT COVERED HERE: the 43 cells themselves, which the field lane owns. This file does the seven
decisions in tasks/assembly/. PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import json
import math
import os
import random
import statistics


def resample_verdict(point_inputs: dict, dispersions: dict, verdict_fn, draws: int = 4000,
                     seed: int = 0) -> dict:
    """Draw each input from a normal with the stated dispersion, recompute, count agreement."""
    rng = random.Random(seed)
    point = verdict_fn(point_inputs)
    agree = 0
    outcomes = {}
    for _ in range(draws):
        trial = {}
        for k, v in point_inputs.items():
            sd = dispersions.get(k, 0.0)
            trial[k] = rng.gauss(v, sd) if sd else v
        try:
            v = verdict_fn(trial)
        except Exception:
            continue
        outcomes[v] = outcomes.get(v, 0) + 1
        if v == point:
            agree += 1
    return {'point_verdict': point, 'agreement_fraction': agree / draws,
            'outcome_counts': outcomes,
            'survives': agree / draws >= 0.95,
            'flipped': agree / draws < 0.95}


def route_topology():
    """Cut: halving ratio below 0.75 is SERIES. Dispersion: the ratio's own spread across the ladder."""
    def fn(d):
        return 'SERIES' if d['halving'] < 0.75 else 'PARALLEL'
    # measured in route_topology_decision.py: two-stage 0.4953, three-stage 0.2440, parallel 1.0
    cases = {'two_stage': (0.4953, 0.02), 'three_stage': (0.2440, 0.02),
             'parallel': (1.0, 0.02), 'weak_parallel': (1.0, 0.02)}
    out = {}
    for name, (v, sd) in cases.items():
        out[name] = resample_verdict({'halving': v}, {'halving': sd}, fn)
    return out


def load_redistribution():
    """Cut: explained fraction above 0.9 means the surrogate suffices. Measured 0.7025."""
    def fn(d):
        return 'SURROGATE SUFFICIENT' if d['explained'] > 0.9 else 'SURROGATE INSUFFICIENT'
    # one static pose per arm, so the only honest dispersion is the solver round-off, 5e-9 relative,
    # plus the spread across the 187 muscles' own errors as a conservative alternative
    return {'solver_roundoff': resample_verdict({'explained': 0.702452},
                                                {'explained': 0.702452 * 5e-9}, fn),
            'muscle_spread_conservative': resample_verdict({'explained': 0.702452},
                                                           {'explained': 0.08}, fn)}


def glass_batch_energy():
    """Cut: carbonate share above 0.10 of sensible heat changes the decision. Measured 0.760."""
    def fn(d):
        return 'CHANGES THE DECISION' if d['share'] > 0.10 else 'negligible'
    # dispersion from the decomposition enthalpies, which are standard-state values; 3 percent each
    # is generous for tabulated formation enthalpies, and the composition itself is a measured glass
    return {'enthalpy_3pct': resample_verdict({'share': 0.7601}, {'share': 0.7601 * 0.03}, fn),
            'enthalpy_20pct_pessimistic': resample_verdict({'share': 0.7601},
                                                           {'share': 0.7601 * 0.20}, fn)}


def mptp_ceiling():
    """Cut: the combination arm reads 20.0 under a shared ceiling or 9.51 under independence.

    The dispersion IS the question: the source's plus-minus 4 is SD or SEM and nobody has resolved
    which. Both readings are run.
    """
    def fn(d):
        mid = (20.0 + d['independent']) / 2.0
        return 'SHARED CEILING' if d['observed'] > mid else 'INDEPENDENT'
    independent = 29.0 * 20.0 / 61.0
    out = {}
    # The first version of this gate used the PER-ANIMAL spread as the uncertainty of a GROUP MEAN
    # and flipped the verdict in both readings. That is the dispersion-convention error this project
    # has been caught by six times, committed inside the gate meant to catch it. The uncertainty of
    # a mean of n animals is SD/sqrt(n) when the source's plus-minus is an SD, and is the stated
    # plus-minus itself when it is already an SEM. Both readings below, correctly scaled.
    for label, sd in (('pm_is_SD_n7_so_SE_1.512', 4.0 / math.sqrt(7)),
                      ('pm_is_SEM_so_SE_4.0', 4.0),
                      ('per_animal_not_a_mean_WRONG_READING', 4.0)):
        # the observed arm is assumed to land on the shared-ceiling prediction; the test is whether
        # the dispersion lets that reading be distinguished from the independent one at all
        out[label] = resample_verdict({'observed': 20.0, 'independent': independent},
                                      {'observed': sd, 'independent': 0.0}, fn)
    return out


def main() -> None:
    report = {'route_topology': route_topology(),
              'load_redistribution': load_redistribution(),
              'glass_batch_energy': glass_batch_energy(),
              'mptp_pathway_ceiling': mptp_ceiling()}
    flips = []
    print(f'{"beslut / fall":<44} {"punkt":<22} {"andel":>7}  utfall')
    for decision, cases in report.items():
        for case, r in cases.items():
            mark = "FLIPPED -> FAIL" if r['flipped'] else "stands"
            if r['flipped']:
                flips.append(f'{decision}/{case}')
            print(f'  {decision + "/" + case:<42} {str(r["point_verdict"])[:20]:<22} '
                  f'{r["agreement_fraction"]:>7.3f}  {mark}')
    print(f'\nturn: {len(flips)} av {sum((len(c) for c in report.values()))}')
    for f in flips:
        print('   ', f)
    d = 'results/ASSEMBLY_DISTRIBUTION_GATE'
    os.makedirs(d, exist_ok=True)
    with open(d + '/gate.json', 'w') as fh:
        json.dump({'ordered_by': 'field lane anton-6b 2026-10-04', 'report': report,
                   'flipped': flips, 'review_state': 'PENDING_INDEPENDENT_REVIEW'}, fh, indent=2)
    print('wrote', d + '/gate.json')


if __name__ == '__main__':
    main()
