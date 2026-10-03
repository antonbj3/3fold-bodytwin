"""A third independent source for the key input: the total cornea measured in one step.

The chain builds the combined cornea by adding an anterior surface computed from two radii to a
posterior surface computed from two more. CASIA reports `CASIA_PR_Real Cyl` -- the total corneal
cylinder measured directly, posterior included. So the sum our model assembles can be checked against a
single measurement of the same thing, which is an external check of the corneal model itself rather
than of a decision made from it.

Then the same substitution that worked for the posterior surface: run the toric decision with the
measured total cornea in place of the assembled one and see whether the residual holds. That makes
three independent sources for the input the decision depends on.

The convention is tested, not assumed. The two devices turned out to report ORTHOGONAL posterior
meridians -- 86 degrees apart in median -- and taking the printed axis cost 0.47 D of apparent accuracy
and nearly produced a false refutation of our own result. So this cell measures the axis relationship
between our assembled cornea and CASIA's total before using it, reports it, and evaluates both the
printed axis and the rotated one rather than choosing silently.

Here the conventions AGREE -- median axis gap 3.55 degrees, within 15 degrees in 67 of 69 eyes and not
one in the 75-105 band -- which is the opposite of the posterior case. Testing instead of assuming was
therefore worth it twice over: the answer is different for the two quantities, and a rule carried over
from one would have been wrong.

One limitation of this cell, stated because it would otherwise look like a second confirmation: the
printed-axis and rotated-axis residuals come out at 0.3112 and 0.3109 D, essentially identical, and
that is invariant BY CONSTRUCTION rather than evidence. The toric lens is aligned to the combined
corneal axis, so rotating the cornea rotates the lens with it and the residual cannot change. The axis
agreement is established by the direct comparison above, not by that pair of numbers.
"""
from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path

import openpyxl

W = Path('.')
R15 = W / 'results/LANE_EYE_OPTICAL_TWIN/r15'
OUT = W / 'results/ASSEMBLY_TOTAL_CORNEA_CHECK'

sys.path.insert(0, str(W / 'tasks/assembly'))
sys.path.insert(0, str(R15))
sys.path.insert(0, str(R15.parent))
from toric_decision import (N_AQUEOUS, from_vector, meridional_refraction,         # noqa: E402
                            surface_power, to_vector)
from hydration_r4 import H0, material                                              # noqa: E402


def axis_gap(a: float, b: float) -> float:
    d = abs(a - b) % 180.0
    return min(d, 180.0 - d)


def main() -> int:
    nc = material(H0)['n_stroma']
    rows = list(openpyxl.load_workbook(next(R15.glob('*.xlsx')), data_only=True).active.values)
    h = [str(x) for x in rows[1]]
    need = ['IOLM_CCT', 'IOLM_AL', 'IOLM_R1', 'IOLM_R2', 'IOLM_AXIS', 'IOLM_PR1', 'IOLM_PR2',
            'IOLM_PRAXIS', 'IOLM_rpmean', 'IOLP', 'IOLT', 'CASIA_POR_AQD post', 'sph', 'cyl',
            'CASIA_PR_Real Cyl', 'CASIA_PR_RealAXIS']
    recs = []
    for n, raw in enumerate(rows[2:], 3):
        d = dict(zip(h, raw))
        if not all(isinstance(d.get(k), (int, float)) for k in need):
            continue
        if not d['IOLT']:
            continue
        d['sheet_row'] = n
        recs.append(d)

    gaps, model_cyl, measured_cyl = [], [], []
    errors = {'assembled_cornea': [], 'measured_total_as_printed': [], 'measured_total_rotated': []}
    out = []
    for d in recs:
        ant = surface_power(d['IOLM_R1'], d['IOLM_R2'], d['IOLM_AXIS'], 1.0, nc)
        post = surface_power(d['IOLM_PR1'], d['IOLM_PR2'], d['IOLM_PRAXIS'], nc, N_AQUEOUS)
        asm_s, asm_c, asm_a = from_vector(*[a + b for a, b in zip(ant, post)])

        gaps.append(axis_gap(asm_a, d['CASIA_PR_RealAXIS']))
        model_cyl.append(abs(asm_c))
        measured_cyl.append(abs(d['CASIA_PR_Real Cyl']))

        variants = {
            'assembled_cornea': (asm_s, asm_c, asm_a),
            'measured_total_as_printed': (asm_s, d['CASIA_PR_Real Cyl'], d['CASIA_PR_RealAXIS']),
            'measured_total_rotated': (asm_s, d['CASIA_PR_Real Cyl'],
                                       (d['CASIA_PR_RealAXIS'] + 90.0) % 180.0),
        }
        row = {'sheet_row': d['sheet_row'], 'measured_cyl_D': d['cyl'],
               'assembled_corneal_cyl_D': round(asm_c, 4),
               'assembled_corneal_axis_deg': round(asm_a, 2),
               'casia_total_cyl_D': d['CASIA_PR_Real Cyl'],
               'casia_total_axis_deg': d['CASIA_PR_RealAXIS'],
               'axis_gap_deg': round(gaps[-1], 2)}
        for name, (cs, cc, ca) in variants.items():
            _, pc, _ = meridional_refraction(d, cs, cc, ca, d['IOLP'], abs(d['IOLT']), ca)
            row[name] = round(pc, 5)
            errors[name].append(abs(pc - d['cyl']))
        out.append(row)

    n = len(recs)
    mae = {k: round(statistics.mean(v), 4) for k, v in errors.items()}
    summary = {
        'question': ('does the cornea our model assembles from two surfaces agree with the total cornea '
                     'measured in one step, and does the decision hold on the measured total'),
        'eyes': n,
        'axis_relationship_tested_not_assumed': {
            'median_axis_gap_deg': round(statistics.median(gaps), 2),
            'mean_axis_gap_deg': round(statistics.mean(gaps), 2),
            'eyes_within_15_deg': sum(1 for g in gaps if g <= 15),
            'eyes_between_75_and_105_deg': sum(1 for g in gaps if 75 <= g <= 105),
            'of': n,
            'why': ('the two devices report orthogonal POSTERIOR meridians, and taking the printed axis '
                    'there cost 0.47 D and nearly produced a false refutation of our own result'),
        },
        'cylinder_magnitudes_D': {
            'assembled_mean': round(statistics.mean(model_cyl), 4),
            'measured_total_mean': round(statistics.mean(measured_cyl), 4),
            'mean_abs_difference': round(statistics.mean(
                [abs(a - b) for a, b in zip(model_cyl, measured_cyl)]), 4),
        },
        'residual_cylinder_mae_D': mae,
        'claim_type': 'information_link',
        'control': ('the cornea assembled by our own model from two surfaces; CASIA total is a third '
                    'independent source for the same input, neither being a facit for the other'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
        'rows': out,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'TOTAL_CORNEA_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    ar = summary['axis_relationship_tested_not_assumed']
    print(f"  eyes {n}")
    print(f"  axis gap: median {ar['median_axis_gap_deg']} deg, within 15 deg in "
          f"{ar['eyes_within_15_deg']} of {n}, between 75 and 105 in "
          f"{ar['eyes_between_75_and_105_deg']}")
    cm = summary['cylinder_magnitudes_D']
    print(f"  cylinder magnitude: assembled {cm['assembled_mean']} D, measured total "
          f"{cm['measured_total_mean']} D, mean difference {cm['mean_abs_difference']} D")
    for k, v in mae.items():
        print(f"    {k:28s} residual MAE {v} D")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
