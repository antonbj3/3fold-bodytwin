"""Measure fixed detector support before placement optimisation.

Surface centers depend on training exits only. No raw coordinates are written.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def candidates(exits, terminal, source):
    n = len(exits)
    training = (np.arange(n) % 2 == 0) & (terminal[:, 1] > 0)
    escaped = terminal[:, 1] > 0
    points = exits[training, :3]
    if len(points) == 0: raise ValueError('No training exits')
    angles = np.arange(12) * 2 * np.pi / 12
    targets = np.asarray(source) + np.column_stack((20 * np.cos(angles), 20 * np.sin(angles), np.zeros(12)))
    centers = np.array([points[np.argmin(np.sum((points - t)**2, axis=1))] for t in targets])
    masks = (np.linalg.norm(exits[:, None, :3] - centers[None, :, :], axis=2) <= 5) & escaped[:, None]
    return centers, masks


def run(captures):
    rows = []; hashes = []
    for filename in captures:
        with np.load(filename, allow_pickle=False) as a:
            centers, masks = candidates(a['exits'], a['terminal'], a['source'])
        overlap = masks.astype(np.int64).T @ masks.astype(np.int64)
        np.fill_diagonal(overlap, 0)
        rows.append(dict(training_counts=masks[::2].sum(axis=0).tolist(), holdout_counts=masks[1::2].sum(axis=0).tolist(),
                         max_pair_overlap=int(overlap.max()), overlapping_pairs=int(np.count_nonzero(np.triu(overlap, 1))),
                         unique_centers=len(np.unique(centers, axis=0)),
                         minimum_center_separation=float(min(np.linalg.norm(centers[i] - centers[j]) for i in range(12) for j in range(i)))))
        hashes.append(dict(centers=hashlib.sha256(centers.tobytes()).hexdigest(), masks=hashlib.sha256(masks.tobytes()).hexdigest()))
    gates = dict(full_repeat=rows[0] == rows[1] and hashes[0] == hashes[1],
                 supported=all(min(r['training_counts'] + r['holdout_counts']) >= 20 for r in rows),
                 disjoint=all(r['max_pair_overlap'] == 0 for r in rows))
    return dict(rows=rows, hashes=hashes, gates=gates, target_azimuth_degrees=list(range(0, 360, 30)),
                scope='Training-selected surface exits, fixed Euclidean acceptance; support diagnostic only, no anatomical or information-gain claim.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('captures', nargs=2)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); result = run(args.captures)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
    raise SystemExit(0 if all(result['gates'].values()) else 2)
