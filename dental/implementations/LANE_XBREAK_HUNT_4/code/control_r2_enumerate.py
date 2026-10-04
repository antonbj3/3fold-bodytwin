"""Additional source-frozen validation: direct finite integer table enumeration, no LP."""
import json, itertools, hashlib, time
from pathlib import Path
start = time.perf_counter()
s = json.loads(Path('raw/R1_SOURCE_COUNTS.json').read_text())
r = json.loads(Path('raw/R2_RESULTS.json').read_text())

def compositions(n, cols):
    for a in range(n + 1):
        for b in range(n - a + 1):
            z = [a, b, n - a - b]
            if all((z[i] <= cols[i] for i in range(3))):
                yield z
worlds = []
for d11 in range(265):
    fixedD = [307 - d11, d11, 0]
    fixedV = [116 - fixedD[0], 264 - d11, 0]
    options = []
    for (d, row) in enumerate([fixedV, fixedD]):
        if min(row) < 0 or sum(row) != s['cold_margins'][d][1]:
            break
        remainder = [s['percussion_margins'][d][i] - row[i] for i in range(3)]
        if min(remainder) < 0:
            break
        options.append([[[remainder[i] - missing[i] for i in range(3)], row, missing] for missing in compositions(s['cold_margins'][d][2], remainder)])
    if len(options) == 2:
        worlds.extend(itertools.product(*options))
checks = []
for case in r['cases']:
    (c, p) = (case['C'], case['P'])
    values = [w[1][c][p] / (w[0][c][p] + w[1][c][p]) for w in worlds if w[0][c][p] + w[1][c][p] > 0]
    endpoints = [min(values), max(values)]
    error = max((abs(endpoints[i] - case['posterior_interval'][i]) for i in range(2)))
    checks.append({'C': c, 'P': p, 'direct_enumeration_endpoints': endpoints, 'LP_endpoints': case['posterior_interval'], 'max_error': error})
out = {'scope': 'Additional independent R2 exact-count validation, frozen original metrics unchanged', 'integer_tables': len(worlds), 'checks': checks, 'gate': 'PASS' if max((v['max_error'] for v in checks)) <= 1e-08 else 'FAIL', 'seconds': time.perf_counter() - start}
Path('raw/R2_INDEPENDENT_ENUMERATION.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
