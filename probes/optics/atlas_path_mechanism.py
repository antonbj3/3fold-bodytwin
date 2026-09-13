"""Aggregate complete transport captures before oxygen or detector design.

No anatomy coordinates or raw arrays are emitted by this observer. The caller
owns capture provenance, the input source point and any explicit report writes.
"""
import hashlib
import numpy as np


def summarize(out, source):
    terminal = out['terminal']
    paths = out['paths']
    exits = out['exits']
    escaped = terminal[:, 1] > 0
    radius = np.linalg.norm(exits[escaped, :3] - np.asarray(source), axis=1)
    edges = [0, 5, 10, 20, 40, 80, np.inf]
    counts = np.histogram(radius, edges)[0]
    sums = terminal[:, :3].sum(axis=0, dtype=np.int64)
    return dict(
        hashes={k: hashlib.sha256(v.tobytes()).hexdigest() for k, v in out.items()},
        photons=len(terminal), absorbed=int(sums[0]), escaped=int(sums[1]),
        residual=int(sums[2]), leaks=int(terminal[:, 3].sum()),
        caps=int(out['counters'][2]), events=int(out['counters'][0]),
        hits=int(out['counters'][1]), escaped_packets=int(escaped.sum()),
        zero_escaped_packets=int((~escaped).sum()),
        shell_edges=[0, 5, 10, 20, 40, 80, 'infinity'],
        shell_counts=counts.tolist(),
        region_visited=np.count_nonzero(paths > 0, axis=0).tolist(),
        region_mean_length=paths.mean(axis=0).tolist(),
        region_max_length=paths.max(axis=0).tolist(),
        finite_nonnegative=bool(np.isfinite(paths).all() and (paths >= 0).all()),
    )


def evaluate(rows, packet):
    return dict(
        full_repeat=all(rows[i] == rows[i + 1] for i in (0, 2)),
        exact_energy=all(r['absorbed'] + r['escaped'] + r['residual'] == r['photons'] * packet for r in rows),
        no_failures=all(r['leaks'] == r['caps'] == r['residual'] == 0 for r in rows),
        finite_nonnegative=all(r['finite_nonnegative'] for r in rows),
        uncensored_proposal=all(r['zero_escaped_packets'] == r['absorbed'] == 0 and r['escaped'] == r['photons'] * packet for r in rows[2:]),
    )


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('captures', nargs=4, help='Original twice, zero absorption twice; NPZ with source and all transport arrays')
    parser.add_argument('--packet', type=int, default=2**30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for filename in args.captures:
        with np.load(filename, allow_pickle=False) as capture:
            rows.append(summarize({k: capture[k] for k in capture.files if k != 'source'}, capture['source']))
    gates = evaluate(rows, args.packet)
    args.output.write_text(json.dumps(dict(rows=rows, gates=gates), indent=2) + '\n')
    raise SystemExit(0 if all(gates.values()) else 2)
