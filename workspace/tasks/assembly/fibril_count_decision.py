"""How many collagen fibrils must act in parallel to carry a stated load, and does the spread matter?

THE DECISION, and it sits in the surgical and laser track rather than in the eye. A cut, a suture
pull-out or a laser-driven impulse is resisted by fibrils in parallel. The twin is asked how much
load a given cross-section carries. The cheap answer divides the load by a mean per-fibril capacity.
The expensive answer says a bundle fails at its weakest members, so the mean overstates what the
bundle holds. If the two answers agree, the mean is enough and nothing more needs measuring. If they
disagree, every bundle-level force the twin states is wrong in a known direction.

THE EXTERNAL DATA, measured, not modelled: 38 single collagen fibril tensile curves from
doi 10.6084/m9.figshare.c.4126559 (CC0, paper PMID 30351303), acquired by LANE_EXTERNAL_SOLVER in
round 31 and already on disk. Each curve is strain against dry-area-referenced stress. Per-specimen
geometry is in the same release: dry area 0.008 to 0.077 um^2 and length 36.6 to 72.7 um.

WHAT IS COMPUTED. Per fibril, the peak of its own curve gives a failure stress; the peak times that
fibril's dry area gives its failure force. The count needed for a stated load follows three ways:

  mean      N = L / mean(F_i)                      -- the cheap reading
  weakest   N such that the N weakest fibrils drawn from this sample carry L
  order     N = L / F_(k) at the k-th order statistic, k chosen as the 5th percentile

THE CONTROL is the mean reading, which sees the same 38 curves and the same geometry and uses only
their average. It is scored, not dismissed: if it lands within the sample spread of the weakest-link
reading, it is the right tool and the extra machinery is waste.

THE FALSIFIER, stated before the run: if the mean and 5th-percentile counts differ by less than the
count's own uncertainty from resampling the 38, the distribution does not matter at this sample size
and this decision collapses to the control. A difference larger than that is the finding.

WHAT THIS IS NOT. A fibril is one scale below tissue. Nothing here claims a tissue modulus, a
boundary load or a surgical force; it claims a fibril count for a load the CALLER states. Combining
fibrils into tissue needs a packing assumption that this file does not supply and does not hide.
PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import json
import math
import os
import random
import statistics
import zipfile

W = '.'
ARCHIVE = f'{W}/results/LANE_EXTERNAL_SOLVER/FIGSHARE_SINGLE_FIBRIL_FILE_11881202_R31.bin'
INSPECTION = f'{W}/results/LANE_EXTERNAL_SOLVER/PUBLIC_DATASET_INSPECTION_R31.json'
LOCATOR = 'doi 10.6084/m9.figshare.c.4126559 (CC0), paper PMID 30351303'


def load_curves() -> dict[str, list[tuple[float, float]]]:
    z = zipfile.ZipFile(ARCHIVE)
    out = {}
    for name in z.namelist():
        if not name.endswith('.txt') or name.startswith('__MACOSX'):
            continue
        rows = []
        for line in z.read(name).decode('utf8', 'replace').splitlines():
            parts = line.replace(',', '.').split()
            if len(parts) != 2:
                continue
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
        if rows:
            out[name.split('/')[-1][:-4]] = rows
    return out


def specimen_areas() -> dict[str, float]:
    """Dry area per specimen, in um^2, from the release's own record list."""
    d = json.load(open(INSPECTION))
    areas = {}
    for rec in d.get('records', []):
        sid = rec.get('specimen_id')
        for key in ('dry_area_um2', 'dry_cross_sectional_area_um2', 'area_um2'):
            if isinstance(rec.get(key), (int, float)):
                areas[sid] = float(rec[key])
                break
    return areas


def failure_forces(curves, areas) -> tuple[dict[str, float], list[str]]:
    """Peak stress times that fibril's own area. Stress is MPa on a um^2 area, so the product is uN."""
    forces, missing = {}, []
    for sid, rows in curves.items():
        peak = max(s for _, s in rows)
        a = areas.get(sid)
        if a is None:
            missing.append(sid)
            continue
        forces[sid] = peak * a          # MPa * um^2 = uN
    return forces, missing


def counts_for_load(load_uN: float, forces: list[float]) -> dict:
    mean_f = statistics.mean(forces)
    srt = sorted(forces)
    # 5th percentile by the sample's own order statistics, no distribution assumed
    k = max(0, int(math.floor(0.05 * (len(srt) - 1))))
    p05 = srt[k]
    # weakest-link: fill the bundle from the weakest end until the load is carried
    total, n_weak = 0.0, 0
    for f in srt:
        if total >= load_uN:
            break
        total += f
        n_weak += 1
    return {'n_mean': load_uN / mean_f, 'n_p05': load_uN / p05,
            'n_weakest_link_fill': n_weak, 'mean_force_uN': mean_f,
            'p05_force_uN': p05, 'min_force_uN': srt[0], 'max_force_uN': srt[-1]}


def resample_spread(load_uN: float, forces: list[float], draws: int = 2000) -> dict:
    rng = random.Random(0)
    ns = []
    for _ in range(draws):
        sample = [rng.choice(forces) for _ in forces]
        ns.append(load_uN / statistics.mean(sample))
    ns.sort()
    lo, hi = ns[int(0.025 * draws)], ns[int(0.975 * draws)]
    return {'n_mean_ci95': [lo, hi], 'n_mean_ci_width': hi - lo}


def main() -> None:
    curves = load_curves()
    areas = specimen_areas()
    forces, missing = failure_forces(curves, areas)
    print(f'curves read: {len(curves)}   with geometry: {len(forces)}   without: {len(missing)}')
    if not forces:
        print('no specimen areas in the release record; the per-fibril force cannot be formed here.')
        print('missing ids:', missing[:6])
        return
    vals = list(forces.values())
    print(f'failure force per fibril, uN: min {min(vals):.4f}  median '
          f'{statistics.median(vals):.4f}  max {max(vals):.4f}  spread {max(vals)/min(vals):.2f}x')
    out = {'source': LOCATOR, 'n_specimens': len(vals), 'review_state': 'PENDING_INDEPENDENT_REVIEW',
           'loads': {}}
    for load in (1.0, 10.0, 100.0):
        c = counts_for_load(load, vals)
        r = resample_spread(load, vals)
        gap = abs(c['n_p05'] - c['n_mean'])
        verdict = ('DISTRIBUTION MATTERS' if gap > r['n_mean_ci_width'] else
                   'mean is sufficient at this sample size')
        print(f'\nload {load:g} uN:')
        print(f'  count from the mean      {c["n_mean"]:.2f}   95 % interval from resampling '
              f'[{r["n_mean_ci95"][0]:.2f}, {r["n_mean_ci95"][1]:.2f}]')
        print(f'  count from 5th pct       {c["n_p05"]:.2f}')
        print(f'  gap {gap:.2f} against interval width {r["n_mean_ci_width"]:.2f}  ->  {verdict}')
        out['loads'][f'{load:g}_uN'] = dict(c, **r, gap=gap, verdict=verdict)
    d = f'{W}/results/ASSEMBLY_FIBRIL_COUNT'
    os.makedirs(d, exist_ok=True)
    with open(d + '/decision.json', 'w') as fh:
        json.dump(out, fh, indent=2)
    print('\nwrote', d + '/decision.json')


if __name__ == '__main__':
    main()
