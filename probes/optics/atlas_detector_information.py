"""Constrained detector information under explicit optical/noise assumptions."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from atlas_detector_support import candidates
from atlas_oxygen_bounds import EXTINCTION, FRACTION, BACKGROUND


def matrices(paths, masks):
    scale = np.log(10.) * 15 / 64500.
    mu = BACKGROUND[None, :] + (scale * (EXTINCTION[:, 0] * .8 + EXTINCTION[:, 1] * .2))[:, None] * FRACTION[None, :]
    dmu = (scale * (EXTINCTION[:, 0] - EXTINCTION[:, 1]))[:, None] * FRACTION[None, :]
    weights = np.exp(-(paths @ mu.T))
    derivative = -(paths @ dmu.T) * weights
    out = []
    for mask in masks.T:
        signal = weights[mask].sum(axis=0)
        if np.any(signal <= 0): raise ValueError('Positive detector support required')
        dlog = derivative[mask].sum(axis=0) / signal
        out.append(np.sqrt(signal)[:, None] * np.column_stack((dlog, np.ones(2))))
    return np.array(out)


def score(matrix, triple):
    return float(np.linalg.svd(matrix[list(triple)].reshape(6, 2), compute_uv=False)[-1])


def run(captures):
    rows = []; hashes = []
    baseline = (0, 4, 8)
    for filename in captures:
        with np.load(filename, allow_pickle=False) as a:
            centers, masks = candidates(a['exits'], a['terminal'], a['source'])
            paths = a['paths']
        counts = [masks[::2].sum(axis=0), masks[1::2].sum(axis=0)]
        distance = np.linalg.norm(centers[:, None] - centers[None, :], axis=2)
        train = matrices(paths[::2], masks[::2]); holdout = matrices(paths[1::2], masks[1::2])
        eligible = np.flatnonzero(counts[0] >= 20).tolist()
        triples = [t for t in itertools.combinations(eligible, 3) if all(distance[i, j] > 10 for i, j in itertools.combinations(t, 2))]
        if not triples: raise ValueError('No feasible triple')
        training_scores = [score(train, t) for t in triples]
        chosen = min(zip(triples, training_scores), key=lambda x: (-x[1], x[0]))[0]
        b = [score(train, baseline), score(holdout, baseline)]
        c = [score(train, chosen), score(holdout, chosen)]
        rows.append(dict(selected_ids=list(chosen), selected_target_angles=[i * 30 for i in chosen],
                         reference_ids=list(baseline), feasible_triples=len(triples), training_candidates=eligible,
                         selected_sigma_min=c, reference_sigma_min=b, gain=[c[i] / b[i] for i in (0, 1)],
                         selected_counts=[v[list(chosen)].tolist() for v in counts],
                         reference_counts=[v[list(baseline)].tolist() for v in counts],
                         minimum_separation=min(float(distance[i, j]) for t in (baseline, chosen) for i, j in itertools.combinations(t, 2)),
                         all_training_scores=[dict(ids=list(t), sigma_min=s) for t, s in zip(triples, training_scores)]))
        hashes.append(dict(jacobians=hashlib.sha256(np.array([train, holdout]).tobytes()).hexdigest(),
                           training_scores=hashlib.sha256(np.array(training_scores).tobytes()).hexdigest()))
    gates = dict(full_repeat=rows[0] == rows[1] and hashes[0] == hashes[1],
                 disjoint=all(r['minimum_separation'] > 10 for r in rows),
                 support=all(min(sum(r['selected_counts'] + r['reference_counts'], [])) >= 20 for r in rows),
                 positive_rank=all(np.isfinite(r['selected_sigma_min'] + r['reference_sigma_min']).all() and min(r['selected_sigma_min'] + r['reference_sigma_min']) > 0 for r in rows),
                 heldout_gain=all(r['gain'][1] > 1 for r in rows))
    return dict(rows=rows, hashes=hashes, gates=gates,
                scope='Conditional Poisson information with fixed illustrative optics and parameter scales; held-out packet halves, no clinical, anatomical or sampling-uncertainty certificate.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('captures', nargs=2)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); result = run(args.captures)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(gates=result['gates'], row={k:v for k,v in result['rows'][0].items() if k != 'all_training_scores'})))
    raise SystemExit(0 if all(result['gates'].values()) else 2)
