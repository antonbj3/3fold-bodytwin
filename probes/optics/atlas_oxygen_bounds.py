"""Illustrative atlas oxygen observable with explicit capped-path tail bounds.

This is a static two-wavelength diagnostic, not clinical calibration. Capture
files and coordinates remain caller-owned; only aggregate results are written.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

EXTINCTION = np.array([[586., 1548.52], [1058., 691.32]])
FRACTION = np.array([.02, .01, 0, .04, .04, 0])
BACKGROUND = np.array([.002, .004, .0004, .004, .004, 0])


def observe(paths, selected, unfinished, saturation):
    scale = np.log(10.) * 15 / 64500.
    blood = scale * (EXTINCTION[:, 0] * saturation + EXTINCTION[:, 1] * (1 - saturation))
    mu = BACKGROUND[None, :] + blood[:, None] * FRACTION[None, :]
    dmu = scale * (EXTINCTION[:, 0] - EXTINCTION[:, 1])[:, None] * FRACTION[None, :]
    weights = np.exp(-(paths @ mu.T))
    lower = weights[selected].sum(axis=0) / len(paths)
    tail = weights[unfinished].sum(axis=0) / len(paths)
    upper = lower + tail
    mean_path = (weights[selected].T @ paths[selected]) / weights[selected].sum(axis=0)[:, None]
    dlog = -(mean_path * dmu).sum(axis=1)
    return dict(saturation=saturation, intensity_lower=lower.tolist(), intensity_upper=upper.tolist(),
                unfinished_upper=tail.tolist(), conditional_ratio=float(lower[0] / lower[1]),
                ratio_lower=float(lower[0] / upper[1]), ratio_upper=float(upper[0] / lower[1]),
                conditional_log_slope=float(dlog[0] - dlog[1])), weights


def run(captures):
    legs = []
    hashes = []
    provenance = []
    counts = []
    for filename in captures:
        with np.load(filename, allow_pickle=False) as a:
            paths = a['paths']; terminal = a['terminal']; exits = a['exits']; source = a['source']
        if paths.shape != (len(terminal), 6) or not np.isfinite(paths).all() or (paths < 0).any():
            raise ValueError('Six finite nonnegative region paths required')
        if (terminal[:, 0] != 0).any() or (terminal[:, 3] != 0).any():
            raise ValueError('Zero absorption and zero leak capture required')
        if not np.all(terminal[:, 1] + terminal[:, 2] == 2**30):
            raise ValueError('Complete escaped/residual proposal ledger required')
        escaped = terminal[:, 1] > 0
        unfinished = terminal[:, 2] > 0
        radius = np.linalg.norm(exits[:, :3] - source, axis=1)
        selected = escaped & (radius >= 5) & (radius < 40)
        if not selected.any(): raise ValueError('Empty preregistered collection shell')
        rows = []; weights = []
        for sat in (.7, .8, .9, 1.):
            row, w = observe(paths, selected, unfinished, sat)
            if sat < 1:
                h = 1e-5
                plus = observe(paths, selected, unfinished, sat + h)[0]['conditional_ratio']
                minus = observe(paths, selected, unfinished, sat - h)[0]['conditional_ratio']
                fd = (np.log(plus) - np.log(minus)) / (2 * h)
                row['slope_relative_error'] = float(abs(fd - row['conditional_log_slope']) / abs(row['conditional_log_slope']))
            rows.append(row); weights.append(w)
        hashes.append(hashlib.sha256(np.array(weights).tobytes()).hexdigest())
        provenance.append({k: hashlib.sha256(v.tobytes()).hexdigest() for k, v in dict(paths=paths, terminal=terminal, exits=exits).items()})
        counts.append(dict(selected=int(selected.sum()), unfinished=int(unfinished.sum()), launched=len(paths)))
        legs.append(rows)
    gates = dict(full_repeat=legs[0] == legs[1] and hashes[0] == hashes[1] and provenance[0] == provenance[1],
                 positive_intervals=all(0 < min(r['intensity_lower']) <= max(r['intensity_upper']) <= 1 for leg in legs for r in leg),
                 separated_monotonic=all(a['ratio_upper'] < b['ratio_lower'] for leg in legs for a, b in zip(leg, leg[1:])),
                 derivative=all(r['slope_relative_error'] <= 1e-7 for leg in legs for r in leg if 'slope_relative_error' in r),
                 matching_published_calibration=False)
    return dict(legs=legs, weight_hashes=hashes, input_hashes=provenance, counts=counts, gates=gates,
                blood_fraction=FRACTION.tolist(), background=BACKGROUND.tolist(), hb_g_dl=15,
                extinction_source='https://omlc.org/spectra/hemoglobin/summary.html',
                matching_calibration_anchors=0,
                scope='Illustrative static oxygen model on public-atlas captures; capped continuation bound only, excludes sampling, geometry and physiological uncertainty.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('captures', nargs=2)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = run(args.captures)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
    raise SystemExit(0 if all(result['gates'].values()) else 2)
