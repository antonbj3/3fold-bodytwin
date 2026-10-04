"""Does the toric gain survive a change of instrument? The decisive control, and it was not run before.

Why it has to be run. The toric decision's reported result is a residual cylinder of 0.282 D with the
measured posterior cornea against 0.626 D with a population estimate, on 69 eyes. The instrument cell
then measured something that puts that straight into question: the two devices in the same workbook
disagree about posterior corneal cylinder with a standard deviation of 0.426 D, which is LARGER than the
0.344 D improvement the decision claims. If the gain came from this quantity, and the quantity is only
known to half a dioptre, the gain may belong to the IOLMaster rather than to the measurement.

The answer is that it survives, and the first run of this cell said the opposite because of a
convention. Taking CASIA's posterior axis as printed gives 0.758 D, worse than the population
estimate, which looked like a clean refutation. The two devices report ORTHOGONAL meridians -- median
axis difference 86 degrees, 65 of 69 eyes between 75 and 105, not one within 15 -- so the printed axis
applies the cylinder 90 degrees out. Aligned, CASIA gives 0.293 D against the IOLMaster's 0.282 D and
the population estimate's 0.626 D. The gain is a property of measuring the posterior cornea, not of one
device. The earlier reading is withdrawn; it was the same axis-convention trap that inflated this
cell's first residual to 5.2 D, caught this time before it was reported.

So this is the equally informed control for an information_link claim: the same chain, the same eyes,
the same outcomes, with the posterior cylinder taken from the OTHER device. Everything else is held
fixed, including the mean posterior power, which is kept from the IOLMaster radii so that only the
disagreed-upon quantity changes. CASIA reports the posterior cylinder in the same sign convention the
cell already uses (its values run -0.89 to -0.06 D against the population constant's -0.30), so no
conversion is introduced.

The reviewed cell is imported, not edited: its helpers and its chain are used exactly as they stand.
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

import openpyxl

W = Path('')
R15 = W / 'results/LANE_EYE_OPTICAL_TWIN/r15'
OUT = W / 'results/ASSEMBLY_TORIC_CROSS_DEVICE'

sys.path.insert(0, str(W / 'tasks/assembly'))
sys.path.insert(0, str(R15))
sys.path.insert(0, str(R15.parent))
from toric_decision import (N_AQUEOUS, POP_POSTERIOR_AXIS_DEG, POP_POSTERIOR_CYL_D,   # noqa: E402
                            from_vector, meridional_refraction, surface_power, to_vector)
from hydration_r4 import H0, material                                                 # noqa: E402


def main() -> int:
    nc = material(H0)['n_stroma']
    rows = list(openpyxl.load_workbook(next(R15.glob('*.xlsx')), data_only=True).active.values)
    h = [str(x) for x in rows[1]]
    need = ['IOLM_CCT', 'IOLM_AL', 'IOLM_R1', 'IOLM_R2', 'IOLM_AXIS', 'IOLM_PR1', 'IOLM_PR2',
            'IOLM_PRAXIS', 'IOLM_rpmean', 'IOLP', 'IOLT', 'CASIA_POR_AQD post', 'sph', 'cyl',
            'CASIA_PR_PostCyl', 'CASIA_PR_PostAXIS']
    recs = []
    for n, raw in enumerate(rows[2:], 3):
        d = dict(zip(h, raw))
        if not all(isinstance(d.get(k), (int, float)) for k in need):
            continue
        if not d['IOLT']:
            continue
        d['sheet_row'] = n
        recs.append(d)

    errors = {'iolmaster_posterior': [], 'casia_posterior': [], 'casia_posterior_as_printed': [],
              'population_estimate': []}
    rows_out = []
    for d in recs:
        ant = surface_power(d['IOLM_R1'], d['IOLM_R2'], d['IOLM_AXIS'], 1.0, nc)
        mean_post_power = (N_AQUEOUS - nc) / (((d['IOLM_PR1'] + d['IOLM_PR2']) / 2) / 1000.0)
        variants = {
            'iolmaster_posterior': surface_power(d['IOLM_PR1'], d['IOLM_PR2'], d['IOLM_PRAXIS'],
                                                 nc, N_AQUEOUS),
            # only the disagreed-upon quantity changes: the mean posterior power is the IOLMaster's
            # CASIA reports the ORTHOGONAL meridian: the two devices' posterior axes differ by 86
            # degrees in median, with 65 of 69 eyes in the 75-105 band and NONE within 15. Taking the
            # axis as printed applies the cylinder 90 degrees off, which is what made the first run of
            # this cell report the gain as instrument-specific. That reading is withdrawn.
            'casia_posterior_as_printed': to_vector(mean_post_power, d['CASIA_PR_PostCyl'],
                                                    d['CASIA_PR_PostAXIS']),
            'casia_posterior': to_vector(mean_post_power, d['CASIA_PR_PostCyl'],
                                         (d['CASIA_PR_PostAXIS'] + 90.0) % 180.0),
            'population_estimate': to_vector(mean_post_power, -POP_POSTERIOR_CYL_D,
                                             POP_POSTERIOR_AXIS_DEG),
        }
        row = {'sheet_row': d['sheet_row'], 'measured_cyl_D': d['cyl']}
        for name, post in variants.items():
            cs, cc, ca = from_vector(*[a + b for a, b in zip(ant, post)])
            _, pc, _ = meridional_refraction(d, cs, cc, ca, d['IOLP'], abs(d['IOLT']), ca)
            row[name] = round(pc, 5)
            errors[name].append(abs(pc - d['cyl']))
        rows_out.append(row)

    mae = {k: round(statistics.mean(v), 4) for k, v in errors.items()}
    n = len(recs)
    better = {k: sum(1 for a, b in zip(errors[k], errors['population_estimate']) if a < b)
              for k in ('iolmaster_posterior', 'casia_posterior', 'casia_posterior_as_printed')}
    device_disagreement = [abs(r['iolmaster_posterior'] - r['casia_posterior']) for r in rows_out]

    summary = {
        'question': 'does the toric gain survive taking the posterior cylinder from the other device',
        'eyes': n,
        'residual_cylinder_mae_D': mae,
        'gain_over_population_estimate_D': {
            'iolmaster_posterior': round(mae['population_estimate'] - mae['iolmaster_posterior'], 4),
            'casia_posterior': round(mae['population_estimate'] - mae['casia_posterior'], 4),
            'casia_posterior_as_printed': round(
                mae['population_estimate'] - mae['casia_posterior_as_printed'], 4)},
        'axis_convention': ('the two devices report orthogonal posterior meridians: median axis '
                            'difference 86 degrees, 65 of 69 eyes between 75 and 105, none within 15. '
                            'The printed-axis variant is retained only to show what that costs.'),
        'eyes_closer_than_population_estimate': {**better, 'of': n},
        'predicted_cylinder_disagreement_between_devices_D': {
            'mean': round(statistics.mean(device_disagreement), 4),
            'max': round(max(device_disagreement), 4)},
        'what_was_held_fixed': ('the chain, the eyes, the outcomes, the anterior cornea, the implanted '
                                'lens and the mean posterior power; only the posterior cylinder '
                                'magnitude and axis change'),
        'claim_type': 'information_link',
        'control': ('the population estimate, as in the original cell, plus the second device as an '
                    'equally informed alternative source for the same quantity'),
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
        'rows': rows_out,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'CROSS_DEVICE_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    print(f"  eyes {n}")
    for k, v in mae.items():
        print(f"    {k:22s} residual MAE {v} D")
    print(f"  gain over population estimate: IOLMaster "
          f"{summary['gain_over_population_estimate_D']['iolmaster_posterior']} D, "
          f"CASIA {summary['gain_over_population_estimate_D']['casia_posterior']} D")
    print(f"  closer than the estimate: {better} of {n}")
    print(f"  the two devices' predictions differ by "
          f"{summary['predicted_cylinder_disagreement_between_devices_D']['mean']} D on average, "
          f"max {summary['predicted_cylinder_disagreement_between_devices_D']['max']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
