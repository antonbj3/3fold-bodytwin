import subprocess
import numpy as np
from freeze import ROOT

def simplify(curve, k):
    body = f'{len(curve)} {k} 10\n' + '\n'.join((' '.join((format(float(v), '.17g') for v in p)) for p in curve)) + '\n'
    result = subprocess.run(['java', '-Xmx128m', '-cp', str(ROOT / 'comparison_java'), 'Compare'], input=body, capture_output=True, text=True, timeout=20)
    if result.returncode:
        raise ValueError('UPSTREAM_JAVA_FAILED: ' + result.stderr)
    lines = result.stdout.splitlines()
    first = next((i for (i, v) in enumerate(lines) if v.startswith('NODES ')), None)
    if first is None:
        raise ValueError('UPSTREAM_NO_NODES: ' + result.stdout)
    n = int(lines[first].split()[1])
    nodes = np.array([[float(x) for x in line.split()] for line in lines[first + 1:first + 1 + n]])
    return (nodes, dict(stdout=result.stdout, stderr=result.stderr, count=n - 1, length_mm=float(np.linalg.norm(np.diff(nodes, axis=0), axis=1).sum()), shortest_segment_mm=float(np.linalg.norm(np.diff(nodes, axis=0), axis=1).min()), count_constraint_met=n - 1 <= k, minimum_length_constraint_met=bool(np.linalg.norm(np.diff(nodes, axis=0), axis=1).min() >= 10), scope='Published OsteoOpt geometric simplifier only; not its full union/chewing optimizer'))
