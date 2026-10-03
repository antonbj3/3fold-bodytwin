"""The second independent measurement per eye turns out to exist: two instruments measured the same eyes.

What this resolves. The decorrelation cell ended on a sharp negative: an additive stage decomposition of
a total observed only once is not identifiable, so the toric margin could not be certified, and the
instrument component could not be stacked because it had been measured on 300 eyes of a different
cohort with zero overlap. That conclusion was right about the 300-eye cohort and wrong about the data --
the surgical workbook carries TWO instruments on the SAME eyes. IOLMaster and CASIA both report central
corneal thickness, anterior chamber depth, lens thickness, white-to-white, posterior corneal cylinder
and total corneal cylinder, on all 89 eyes and all 69 toric eyes.

Why that is the key and not just more data. The disagreement between two devices measuring the same eye
is an instrument error component that is MEASURED per unit, carries no shared term with the outcome, and
sits on the same eyes as the outcomes. It is therefore admissible in the stage matrix in a way the
300-eye cohort never was, and it is not constructed by subtraction from the total, which is what made
the earlier toric correlation an artefact.

What is measured here, and nothing is assumed about it: the per-eye disagreement in dioptres of
cylinder for the two corneal quantities the toric decision actually consumes, plus the device
disagreement on the three biometric lengths. The posterior cylinder is converted from the IOLMaster's
two radii using the SAME constants and the same stromal index the toric decision uses, imported rather
than restated, because inventing a second conversion chain is what produced a 5.2 D residual once
already.

A systematic offset is expected and is reported rather than removed: the two devices differ by 9.90 um
in corneal thickness with a spread of 6.26, which is a bias and not noise. In the graph engine's terms
that is a partially shared systematic rather than an independent component, so the mean and the spread
are reported separately and the composition is left to the verifier.
"""
from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path

import numpy as np
import openpyxl

W = Path('.')
R15 = W / 'results/LANE_EYE_OPTICAL_TWIN/r15'
ENGINE = Path('the public staging tree/3fold-graph-engine/src')
OUT = W / 'results/ASSEMBLY_INSTRUMENT_COMPONENT'
SPEC_D = 0.5

sys.path.insert(0, str(R15))
sys.path.insert(0, str(R15.parent))
sys.path.insert(0, str(ENGINE))
from hydration_r4 import material, H0                                        # noqa: E402
from graph_engine.stage_decorrelation_verifier import (                      # noqa: E402
    covariance_aware_margin, heavy_tail_index, stage_decorrelation_verifier)

N_AQUEOUS = 1.336


def cyl_from_radii(r1_mm: float, r2_mm: float, n_from: float, n_to: float) -> float:
    """Cylinder in dioptres from two principal radii, with the same convention as the toric cell."""
    k1 = (n_to - n_from) / (r1_mm / 1000.0)
    k2 = (n_to - n_from) / (r2_mm / 1000.0)
    return k1 - k2


def jvec(cyl_D: float, axis_deg: float) -> tuple[float, float]:
    """The two cylinder components of the power vector, so axis and magnitude cannot be confused.

    A signed cylinder cannot be differenced across instruments without its axis: these two devices
    report ORTHOGONAL posterior meridians, 86 degrees apart in median, so a raw difference of the
    signed cylinders measures the convention rather than the disagreement. That was the first version
    of this cell and it reported 0.4263 D of posterior-cylinder disagreement, which is withdrawn.
    """
    a = math.radians(2.0 * axis_deg)
    return (-(cyl_D / 2.0) * math.cos(a), -(cyl_D / 2.0) * math.sin(a))


def main() -> int:
    nc = material(H0)['n_stroma']
    rows = list(openpyxl.load_workbook(next(R15.glob('*.xlsx')), data_only=True).active.values)
    h = [str(x) for x in rows[1]]
    recs = [dict(zip(h, r)) for r in rows[2:]]
    num = lambda v: isinstance(v, (int, float))                              # noqa: E731

    lengths = []
    for name, a, b, unit in (('CCT', 'IOLM_CCT', 'CASIA_PR_CCT', 'um'),
                             ('ACD', 'IOLM_ACD', 'CASIA_PR_ACD', 'mm'),
                             ('LT', 'IOLM_LT', 'CASIA_PR_LT', 'mm'),
                             ('WTW', 'IOLM_WTW', 'CASIA_PR_WTW', 'mm')):
        d = [r[a] - r[b] for r in recs if num(r.get(a)) and num(r.get(b))]
        lengths.append({'quantity': name, 'unit': unit, 'eyes_with_both_devices': len(d),
                        'mean_difference': round(statistics.mean(d), 5),
                        'sd_of_difference': round(statistics.stdev(d), 5),
                        'max_abs_difference': round(max(abs(x) for x in d), 5),
                        'reading': 'a nonzero mean is a device bias, i.e. a shared systematic, not noise'})

    toric = [r for r in recs if num(r.get('IOLT')) and r['IOLT']]
    stages, kept = [], []
    for r in toric:
        if not all(num(r.get(k)) for k in ('IOLM_PR1', 'IOLM_PR2', 'IOLM_PRAXIS',
                                           'CASIA_PR_PostCyl', 'CASIA_PR_PostAXIS')):
            continue
        post_iolm = cyl_from_radii(r['IOLM_PR1'], r['IOLM_PR2'], nc, N_AQUEOUS)
        ji = jvec(post_iolm, r['IOLM_PRAXIS'])
        jc = jvec(r['CASIA_PR_PostCyl'], (r['CASIA_PR_PostAXIS'] + 90.0) % 180.0)
        stages.append([ji[0] - jc[0], ji[1] - jc[1]])
        kept.append(r['IOLM_Eye'] if 'IOLM_Eye' in r else None)

    X = np.asarray(stages)
    v = stage_decorrelation_verifier(X, rho_tol=0.15, min_units=20, n_boot=400, seed=0)
    cam = covariance_aware_margin(X, SPEC_D, ci_pctl=99.0, n_boot=300, seed=0, min_units=20)
    hti = heavy_tail_index(X, frac=0.15, hill_thresh=0.35)

    summary = {
        'finding': ('the second independent measurement per eye exists in the surgical workbook: two '
                    'instruments measured the same eyes on six shared quantities'),
        'corrects': ('the decorrelation cell declared the instrument component unstackable because it had '
                     'been measured on 300 eyes of a different cohort; that was true of that cohort and '
                     'false of the data'),
        'why_admissible': ('a between-device disagreement is measured per unit, shares no term with the '
                           'outcome, and sits on the same eyes as the outcomes, so it is not constructed '
                           'by subtraction from the total the way the withdrawn toric stages were'),
        'length_quantities_both_devices': lengths,
        'cylinder_stage_matrix': {
            'units_of_each_column': 'D of cylinder',
            'columns': ['posterior_cylinder_disagreement_J0', 'posterior_cylinder_disagreement_J45'],
            'convention': ('power-vector components with CASIA rotated to the same meridian; the two '
                           'devices print orthogonal posterior axes, 86 degrees apart in median'),
            'withdrawn_first_version': ('a raw difference of the signed cylinders gave 0.4263 D and '
                                        'measured the convention, not the disagreement'),
            'eyes': int(X.shape[0]),
            'per_stage_mean_D': [round(float(m), 5) for m in X.mean(axis=0)],
            'per_stage_sd_D': [round(float(s), 5) for s in X.std(axis=0, ddof=1)],
            'conversion': ('IOLMaster radii converted with the stromal index and aqueous index the toric '
                           'decision itself imports, not a second chain'),
        },
        'decorrelation_verifier': v,
        'covariance_aware_margin': cam,
        'heavy_tail': hti,
        'spec_D': SPEC_D,
        'claim_type': 'capability',
        'control': ('the module own endpoints, spec/sum(sd) at rho=1 and spec/sqrt(sum sd^2) at rho=0; '
                    'and the realised residual cylinder of the toric decision as the external check'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'INSTRUMENT_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=str))

    for L in lengths:
        print(f"  {L['quantity']:5s} n={L['eyes_with_both_devices']:3d} "
              f"bias {L['mean_difference']:+9.4f} sd {L['sd_of_difference']:8.4f} {L['unit']}")
    print(f"  cylinder stages: {X.shape[0]} eyes, sd "
          f"{[round(float(s), 4) for s in X.std(axis=0, ddof=1)]} D")
    print(f"  rho_hat {v.get('rho_hat')}  {v.get('verdict')} -> {v.get('composition')}  "
          f"verified={v.get('decorrelation_verified')}")
    print(f"  margin {cam.get('margin')}  bracket {cam.get('rss_margin')} .. {cam.get('sum_margin')}  "
          f"({cam.get('verdict')})   heavy tail {hti.get('heavy_tail')}")
    print(f"  written: {OUT / 'INSTRUMENT_V1.json'}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
