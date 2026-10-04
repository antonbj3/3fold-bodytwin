#!/usr/bin/env python3
"""A decision outside the eye: which face of a tubular cell is at complement risk, on measured data?

Why this replaces an axis rather than adding one. The complement window this repo carried rested on a
C3b recycling fraction, and LANE_IMMUNE_OVERACTIVATION r1 found rho_status =
UNMEASURED_MODEL_CLOSURE -- the axis is a closure of the model, not a measured quantity, and the lysis
thresholds behind it are synthetic. Two edges were downgraded from TIGHT the same day. A decision
needs an axis someone has measured, and the swarm supplied one.

The measured axis. Complement C3 from human proximal tubular epithelial cells in a Transwell is
apically dominated at rest, basolateral-over-apical ratio 0.45 +/- 0.16, and transferrin in the lumen
moves it to 0.93 +/- 0.24 (Tang et al., Am J Kidney Dis 2001;37(1):94-103, doi
10.1053/ajkd.2001.20593, PMID 11136173). So the asymmetry is not a constant -- it is a quantity that
shifts under a named, published intervention.

The mechanism behind it is measured too, and it is asymmetric by construction rather than by
assumption. The only membrane-bound complement regulator in tubule is Crry, polarised to the
basolateral face: immunogold labelling significantly higher basolaterally than in the brush border
(P < 0.001), the 53 kDa band present in the basolateral fraction and absent from the brush border
fraction, controlled against an apical marker; DAF and CD59 are not expressed in tubule at all
(Thurman et al., J Clin Invest 2006;116(2):357-68, doi 10.1172/JCI24521, PMID 16444293). The apical
face is instead regulated by fluid-phase factor H, and blocking factor H raises activation on BOTH
faces while Crry protects only basolaterally (Renner et al., J Immunol 2010;185(5):3086-94, doi
10.4049/jimmunol.1000111, PMID 20675597).

What this coordinator recomputed before building on it. The shift is 0.48 with a combined dispersion
of 0.2884, so the two states are 1.6641 sigma apart -- real but not overwhelming, and the decision
must not claim more. At rest the ratio is below 1 by 3.438 sigma, which is the firm part. With
transferrin it is 0.292 sigma from 1, which means the loaded state is indistinguishable from
symmetric. So the decision can say which face dominates at rest and must refuse to say it under
luminal protein load.

The control is a single unpolarised regulator, which is what a twin without this measurement would
assume: then the ratio is 1 by construction and no face is ever at differential risk. The falsifier
is the 1.6641 sigma separation itself -- if an independent measurement of the same ratio under the
same intervention lands within one combined dispersion of the resting value, the shift is not there
and this decision has no axis.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

W = Path('.')
OUT = W / 'results/ASSEMBLY_POLARISED_COMPLEMENT'

REST_BA, REST_SD = 0.45, 0.16
LOADED_BA, LOADED_SD = 0.93, 0.24
# Below this many dispersions from 1 the faces are not separable and no face is named.
MIN_SIGMA_FROM_SYMMETRY = 2.0


def sigma_from_symmetry(ratio: float, sd: float) -> float:
    return abs(1.0 - ratio) / sd if sd else float('inf')


def decide(ratio: float, sd: float, label: str) -> dict:
    s = sigma_from_symmetry(ratio, sd)
    out = {'state': label, 'basolateral_over_apical': ratio, 'dispersion': sd,
           'sigma_from_symmetry': s}
    if s < MIN_SIGMA_FROM_SYMMETRY:
        out['decision'] = 'NO_FACE_NAMED'
        out['reason'] = (f'the ratio sits {s:.4f} dispersions from 1, below {MIN_SIGMA_FROM_SYMMETRY}; '
                         f'the two faces are not separable in this state')
        return out
    face = 'APICAL' if ratio < 1 else 'BASOLATERAL'
    out.update({
        'decision': 'FACE_AT_RISK',
        'decision_value_face': face,
        'regulator_on_that_face': ('fluid-phase factor H, not membrane bound'
                                   if face == 'APICAL' else 'Crry, membrane bound and polarised here'),
        'reading': (f'C3 deposition is {1 / ratio:.4f} times higher on the {face.lower()} face'
                    if ratio < 1 else f'{ratio:.4f} times higher basolaterally'),
    })
    return out


def main() -> int:
    rows = [decide(REST_BA, REST_SD, 'rest'), decide(LOADED_BA, LOADED_SD, 'luminal transferrin')]
    shift = LOADED_BA - REST_BA
    comb = math.hypot(REST_SD, LOADED_SD)
    summary = {
        'question': 'which face of a tubular cell carries the complement risk, and in which state',
        'claim_type': 'information_link',
        'measured_axis': {
            'quantity': 'basolateral over apical C3 deposition ratio',
            'rest': [REST_BA, REST_SD], 'luminal_transferrin': [LOADED_BA, LOADED_SD],
            'source': 'Tang et al., Am J Kidney Dis 2001;37(1):94-103, doi 10.1053/ajkd.2001.20593, PMID 11136173',
            'why_this_axis': ('it replaces the C3b recycling fraction, which LANE_IMMUNE_OVERACTIVATION '
                              'r1 found to be an unmeasured model closure with synthetic thresholds'),
        },
        'shift': shift, 'combined_dispersion': comb, 'separation_sigma': shift / comb,
        'asymmetric_mechanism': {
            'basolateral': 'Crry, membrane bound, immunogold labelling higher than brush border at P < 0.001, '
                           '53 kDa band absent from the brush border fraction; doi 10.1172/JCI24521, PMID 16444293',
            'apical': 'fluid-phase factor H; blocking it raises activation on both faces while Crry protects '
                      'only basolaterally; doi 10.4049/jimmunol.1000111, PMID 20675597',
            'not_expressed_in_tubule': ['DAF', 'CD59'],
        },
        'control': ('a single unpolarised regulator, which makes the ratio 1 by construction and no face '
                    'ever at differential risk; that is what a twin without this measurement assumes'),
        'falsifier': (f'the separation is {shift / comb:.4f} sigma. If an independent measurement of the '
                      f'same ratio under luminal protein load lands within one combined dispersion '
                      f'({comb:.4f}) of the resting {REST_BA}, the shift is not there and this decision '
                      f'has no axis'),
        'refusal': {'min_sigma_from_symmetry': MIN_SIGMA_FROM_SYMMETRY,
                    'states_refused': [r['state'] for r in rows if r['decision'] == 'NO_FACE_NAMED']},
        'rows': rows,
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
        'no_claim_of_biological_validation': True,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'DECISION_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f'  measured axle: B/A {REST_BA} +/- {REST_SD} i vila, {LOADED_BA} +/- {LOADED_SD} med transferrin; skifte {shift:.4f}, kombinerad spridning {comb:.4f}, separation {shift / comb:.4f} sigma')
    for r in rows:
        if r['decision'] == 'NO_FACE_NAMED':
            print(f"  {r['state']:22} {r['sigma_from_symmetry']:7.4f} sigma from 1  INGEN SIDA NAMNGES")
        else:
            print(f"  {r['state']:22} {r['sigma_from_symmetry']:7.4f} sigma from 1  {r['decision_value_face']} i risk, {r['reading']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
