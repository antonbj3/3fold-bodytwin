"""Decide whether a drug and a mechanical protocol share one final step, before spending animals.

THE DECISION. Two interventions reduce infarct in the same preparation: a mechanical one applied at
reperfusion, and a drug that blocks the mitochondrial permeability transition pore. The question that
decides whether a combination arm is worth running is whether they share the final step. Two readings
exist and they are not opinions, because they predict different numbers for the same unrun arm:

    SHARED CEILING   both act through the pore, so the combination cannot beat the better single arm
    INDEPENDENT      they act on separable steps, so the surviving fractions multiply

One scalar comes out of this file: the infarct percent the combination arm must show for each reading,
and the number of standard deviations between them, which is what decides whether the arm is worth
running and at what group size.

WHERE THE NUMBERS COME FROM, all four from one lineage in one preparation, carried into
`results/BT-CONN-MITOSTRESS--Q036/RESULTS.md` section 1b with their locators:

    control, no intervention                     61 +/- 6 % of area at risk
    postconditioning, 4 x (1 min / 1 min)        29 +/- 4 %
    ischemic preconditioning                     18 +/- 4 %
    NIM-811, 5 mg/kg IV, 1 min before reflow     20 +/- 4 %
        Argaud L et al., Circulation 2005;111:116-121.
        PMID 15642769, DOI 10.1161/01.CIR.0000151290.04952.3B
        open-chest NZW rabbit, 30 min occlusion, 4 h reperfusion, infarct by TTC
    mechanism, NIM-811 as a cyclophilin-D/mPTP inhibitor, 10 mg/kg IV:
        Argaud L et al., J Mol Cell Cardiol 2005. PMID 15698843, DOI 10.1016/j.yjmcc.2004.12.001

UNIT CONVENTION, stated because the same word means two things in this project. Every percent here is
a fraction of the ANATOMICAL AREA AT RISK, measured at 4 h by TTC in rabbit. The source report records
a collision on exactly this word: the paired cell's percents are within-cell optical signal fractions
over 300 s, and the two are not commensurable because the denominators are different physical things.
Nothing in this file compares across that boundary. The dose is mg/kg of body weight, not micromolar;
the source report also records that the cell's `nim811_concentration_uM` field has no provenance in
this lineage, so no concentration is used anywhere below.

THE EQUALLY INFORMED CONTROL is the practice that reads these same four numbers and takes the best
single arm as the attainable floor -- the dearer path that already works, since both single arms are
published and need no new animals. It predicts 20 % and stops. This file's content is that the rival
reading is already nearly excluded by data in hand, which is a decision reached without spending an
animal, and that the group size needed to close it is small.

FALSIFIER, written before the arm is run. If the combination arm measures at or below 14.8 % of area
at risk -- the midpoint of the two predictions -- the shared-ceiling reading is dead and this file is
wrong. If it measures above 29 %, i.e. worse than postconditioning alone, both readings are wrong and
an interaction term is required that neither carries. A second falsifier needs no new animals: if the
NIM-811 arm was run against its own control group whose mean is below 29 %, the independent prediction
rises to meet the ceiling prediction and the separation this file reports vanishes; the sweep below
reports exactly where that happens.

review_state: PENDING_INDEPENDENT_REVIEW. No clinical or veterinary recommendation is implied.
"""
from __future__ import annotations

import json
import math
import os

CONTROL_PCT = 61.0          # % of area at risk, PMID 15642769
CONTROL_SD = 6.0
POSTCOND_PCT = 29.0
POSTCOND_SD = 4.0
PRECOND_PCT = 18.0
PRECOND_SD = 4.0
NIM811_PCT = 20.0
NIM811_SD = 4.0
NIM811_DOSE_MG_KG = 5.0     # IV bolus, 1 min before reflow

Z_ALPHA_TWO_SIDED_05 = 1.959964
Z_POWER_80 = 0.841621


def protected_fraction(treated_pct: float, control_pct: float) -> float:
    """Fraction of the control infarct that the intervention removed. Dimensionless."""
    return 1.0 - treated_pct / control_pct


def independent_action_prediction(treated_a_pct: float, treated_b_pct: float,
                                  control_pct: float) -> float:
    """If the two act on separable steps the surviving fractions multiply.

    control * (A/control) * (B/control) = A*B/control, so the prediction is a closed form in the three
    measured means and needs no fitted parameter. Written out because the closed form is what makes the
    cancellation question below answerable: the prediction is driven by the RATIO of the control to the
    single arms, not by any difference that two terms could silently cancel.
    """
    return treated_a_pct * treated_b_pct / control_pct


def shared_ceiling_prediction(treated_a_pct: float, treated_b_pct: float) -> float:
    """If both act through one final step, the combination cannot beat the better single arm."""
    return min(treated_a_pct, treated_b_pct)


def pooled_sd(sd_a: float, sd_b: float) -> float:
    return math.sqrt((sd_a ** 2 + sd_b ** 2) / 2.0)


def group_size_for(delta_pct: float, sd: float,
                   z_alpha: float = Z_ALPHA_TWO_SIDED_05,
                   z_power: float = Z_POWER_80) -> int:
    """Animals per group to separate two means by `delta_pct` at the stated error rates."""
    if delta_pct <= 0.0:
        return 10 ** 9
    n = 2.0 * (z_alpha + z_power) ** 2 * sd ** 2 / delta_pct ** 2
    return max(2, math.ceil(n))


def control_at_which_readings_coincide(treated_b_pct: float) -> float:
    """A*B/C = min(A,B) has the solution C = max(A,B). Below that control the two readings agree and
    the experiment decides nothing, so this is the admissibility gate on the control group, not a
    curiosity."""
    return treated_b_pct


def decide(control_pct: float = CONTROL_PCT,
           arm_a_pct: float = POSTCOND_PCT, arm_a_sd: float = POSTCOND_SD,
           arm_b_pct: float = NIM811_PCT, arm_b_sd: float = NIM811_SD) -> dict:
    indep = independent_action_prediction(arm_a_pct, arm_b_pct, control_pct)
    ceil = shared_ceiling_prediction(arm_a_pct, arm_b_pct)
    sd = pooled_sd(arm_a_sd, arm_b_sd)
    gap = ceil - indep
    return {
        'control_pct_of_aar': control_pct,
        'independent_action_prediction_pct_of_aar': indep,
        'shared_ceiling_prediction_pct_of_aar': ceil,
        'gap_percentage_points': gap,
        'pooled_sd_percentage_points': sd,
        'separation_in_pooled_sd': gap / sd,
        'animals_per_group_80pct_power': group_size_for(gap, sd),
        'falsifier_midpoint_pct_of_aar': (indep + ceil) / 2.0,
        'verdict': ('SHARED CEILING is the admissible reading; the independent prediction lies below '
                    'the best infarct ever measured in this preparation'
                    if indep < PRECOND_PCT else
                    'both readings lie inside the measured range; the combination arm is required'),
    }


def main() -> None:
    d = decide()
    print('DECISION: do a pore-blocking drug and a mechanical protocol share one final step?')
    print(f'  preparation: NZW rabbit, 30 min occlusion, 4 h reperfusion, TTC, % of area at risk')
    print(f'  control                                  : {CONTROL_PCT:.0f} +/- {CONTROL_SD:.0f} %')
    print(f'  postconditioning 4 x (1 min / 1 min)     : {POSTCOND_PCT:.0f} +/- {POSTCOND_SD:.0f} %'
          f'   protected {protected_fraction(POSTCOND_PCT, CONTROL_PCT):.4f}')
    print(f'  NIM-811 {NIM811_DOSE_MG_KG:.0f} mg/kg IV before reflow     : '
          f'{NIM811_PCT:.0f} +/- {NIM811_SD:.0f} %   protected '
          f'{protected_fraction(NIM811_PCT, CONTROL_PCT):.4f}')
    print(f'  ischemic preconditioning (benchmark)     : {PRECOND_PCT:.0f} +/- {PRECOND_SD:.0f} %'
          f'   protected {protected_fraction(PRECOND_PCT, CONTROL_PCT):.4f}')

    print('\nTHE TWO PREDICTIONS FOR THE UNRUN COMBINATION ARM')
    print(f'  shared ceiling, one final step   : {d["shared_ceiling_prediction_pct_of_aar"]:.2f} % of AAR')
    print(f'  independent action, steps multiply: {d["independent_action_prediction_pct_of_aar"]:.2f} % of AAR')
    print(f'  gap {d["gap_percentage_points"]:.2f} percentage points = '
          f'{d["separation_in_pooled_sd"]:.2f} pooled SD '
          f'(pooled SD {d["pooled_sd_percentage_points"]:.2f} pp)')
    print(f'  animals per group at 80 % power, two-sided 0.05, IF the +/- are SDs: '
          f'{d["animals_per_group_80pct_power"]}')
    print('\n  DISPERSION CONVENTION, which changes the answer by a factor and must not be guessed.')
    print('  Cardioprotection papers of this period usually print mean +/- SEM, not SD. The source')
    print('  report carries the +/- without saying which, so both readings are computed and the')
    print('  CONSERVATIVE one is what the action below is sized on. A group size of 3 is the tell that')
    print('  the optimistic reading cannot be right for a whole-animal infarct endpoint.')
    for n_assumed in (6, 7, 8, 10):
        sd_true = d['pooled_sd_percentage_points'] * math.sqrt(n_assumed)
        print(f'    if +/- 4 is SEM on n = {n_assumed:2d} -> true SD {sd_true:5.2f} pp, separation '
              f'{d["gap_percentage_points"] / sd_true:4.2f} SD, n/group '
              f'{group_size_for(d["gap_percentage_points"], sd_true):3d}')
    sd_sem7 = d['pooled_sd_percentage_points'] * math.sqrt(7)
    n_conservative = group_size_for(d['gap_percentage_points'], sd_sem7)
    print(f'  SIZED ON: n = {n_conservative} per group, the SEM reading at the lineage\'s own n of 7.')

    print('\nTHE DECISION THAT COSTS NO ANIMALS')
    margin_sd = (PRECOND_PCT - d['independent_action_prediction_pct_of_aar']) / PRECOND_SD
    margin_sem = (PRECOND_PCT - d['independent_action_prediction_pct_of_aar']) / (PRECOND_SD * math.sqrt(7))
    print(f'  the independent prediction, {d["independent_action_prediction_pct_of_aar"]:.2f} %, lies '
          f'{margin_sd:.2f} SD below the lowest')
    print(f'  (and {margin_sem:.2f} SD below it on the conservative SEM reading, so this leg is '
          'suggestive, not decisive)')
    print(f'  infarct ever measured in this preparation, preconditioning at {PRECOND_PCT:.0f} +/- '
          f'{PRECOND_SD:.0f} %.')
    print('  So the independent reading requires the combination to outperform every published arm in')
    print('  the lineage, including the one that needs advance warning of the occlusion. That is the')
    print('  cheapest available refutation and it is already in hand.')
    print(f'  VERDICT: {d["verdict"]}')
    print('  ACTION: run the combination arm as a CEILING test, sized for the '
          f'{d["gap_percentage_points"]:.1f} pp contrast, not as an additivity test.')

    print('\nCONTROL, equally informed: the best-single-arm reading. It sees the same four means and')
    print(f'  predicts {d["shared_ceiling_prediction_pct_of_aar"]:.0f} % with no interaction model. It '
          'reaches the same verdict here, and that is the')
    print('  point: the content is the SEPARATION and the group size, which it cannot supply, and the')
    print('  control group gate below, which tells it when its own answer stops being decidable.')

    print('\nCANCELLATION SCREEN, and the repair response. The headline gap could be an artefact of one')
    print('  term: the prediction is A*B/C, so a control mean that is too high inflates the gap. Scale')
    print('  the control toward and away from its measured ideal and read whether the headline degrades.')
    sweep = []
    for c in (45.0, 50.0, 55.0, 61.0, 67.0, 73.0, 80.0):
        s = decide(control_pct=c)
        sweep.append({'control_pct': c,
                      'independent_pct': s['independent_action_prediction_pct_of_aar'],
                      'separation_sd': s['separation_in_pooled_sd'],
                      'n_per_group': s['animals_per_group_80pct_power']})
        print(f'    control {c:5.1f} % -> independent {s["independent_action_prediction_pct_of_aar"]:6.2f} %'
              f'  separation {s["separation_in_pooled_sd"]:5.2f} SD   n/group '
              f'{s["animals_per_group_80pct_power"]:3d}')
    c_crit = control_at_which_readings_coincide(POSTCOND_PCT)
    print(f'  the two readings coincide only at a control of {c_crit:.0f} %, which is the '
          'postconditioning arm')
    print('  itself. The measured control stands 32 percentage points above that, so the separation is')
    print('  carried by a measured difference and not by two terms cancelling. The repair does not make')
    print('  the headline worse: over a control range of 45 to 80 % the separation stays above 1.5 SD.')

    print('\nWHAT THIS FILE REFUSES TO DO. It does not convert mg/kg to a free micromolar concentration')
    print('  at the inner membrane; the molecule is heavily protein bound and no partition coefficient')
    print('  is in evidence. It does not compare any percent here with an optical signal percent from')
    print('  another cell: the denominators are different physical things.')

    out = {
        'decision': d,
        'measured_inputs_pct_of_area_at_risk': {
            'control': CONTROL_PCT, 'control_sd': CONTROL_SD,
            'postconditioning': POSTCOND_PCT, 'postconditioning_sd': POSTCOND_SD,
            'preconditioning_benchmark': PRECOND_PCT, 'preconditioning_sd': PRECOND_SD,
            'nim811': NIM811_PCT, 'nim811_sd': NIM811_SD,
            'nim811_dose_mg_per_kg': NIM811_DOSE_MG_KG,
        },
        'protected_fractions': {
            'postconditioning': protected_fraction(POSTCOND_PCT, CONTROL_PCT),
            'nim811': protected_fraction(NIM811_PCT, CONTROL_PCT),
            'preconditioning': protected_fraction(PRECOND_PCT, CONTROL_PCT),
        },
        'independent_prediction_sd_below_best_measured_arm': margin_sd,
        'independent_prediction_sd_below_best_measured_arm_sem_reading': margin_sem,
        'control_at_which_readings_coincide_pct': c_crit,
        'control_sweep_repair_response': sweep,
        'animals_per_group_if_pm_is_sem_n7': n_conservative,
        'true_sd_if_pm_is_sem_n7_pp': sd_sem7,
        'separation_in_sd_if_pm_is_sem_n7': d['gap_percentage_points'] / sd_sem7,
        'dispersion_convention_unresolved': True,
        'source': {'pmid': ['15642769', '15698843'],
                   'doi': ['10.1161/01.CIR.0000151290.04952.3B', '10.1016/j.yjmcc.2004.12.001'],
                   'locator': 'results/BT-CONN-MITOSTRESS--Q036/RESULTS.md section 1b',
                   'preparation': 'open-chest NZW rabbit, 30 min occlusion, 4 h reperfusion, TTC',
                   'denominator': 'anatomical area at risk'},
        'falsifier': ('combination arm at or below the midpoint prediction kills the shared-ceiling '
                      'reading; a combination above the postconditioning arm kills both readings; a '
                      'separate NIM-811 control group below 29 % removes the separation'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    dd = 'results/ASSEMBLY_MPTP_PATHWAY_CEILING'
    os.makedirs(dd, exist_ok=True)
    with open(dd + '/decision.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('\nwrote', dd + '/decision.json')


if __name__ == '__main__':
    main()
