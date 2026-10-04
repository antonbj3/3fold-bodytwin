"""Accept independently measured signed losses; never invent missing measurements."""
import argparse, csv, json, math
from pathlib import Path
REQUIRED = ['independent_block_id', 'guide', 'planned_gap_mm', 'measured_min_tool_gap_mm', 'metrology_loss_upper_mm', 'source_locator', 'measurement_kind', 'population_transfer_verified']

def required_blocks(x=0.01, alpha=0.05, candidates=1):
    n = math.ceil(math.log(alpha / candidates) / math.log1p(-x))
    while 1 - (alpha / candidates) ** (1 / n) >= x:
        n += 1
    return n

def tail_upper(n, alpha=0.05, candidates=1):
    return 1 - (alpha / candidates) ** (1 / n)

def calibrate(rows, x=0.01, alpha=0.05, candidates=4, anatomy_budget=None):
    seen = set()
    groups = {}
    if not rows:
        return {'status': 'UNKNOWN', 'reason': 'No independent measured signed losses', 'required_blocks_per_guide': required_blocks(x, alpha, candidates)}
    for row in rows:
        missing = [k for k in REQUIRED if k not in row or row[k] == '']
        if missing:
            raise ValueError('Missing measured fields: ' + ','.join(missing))
        if row['measurement_kind'] != 'INDEPENDENT_PHYSICAL':
            raise ValueError('Synthetic/source-derived observations cannot calibrate an empirical tail')
        if str(row['population_transfer_verified']).lower() != 'true':
            raise ValueError('Population transfer not established')
        key = (row['guide'], row['independent_block_id'])
        if key in seen:
            raise ValueError('Repeated block; use its maximum loss once')
        seen.add(key)
        err = float(row['metrology_loss_upper_mm'])
        if err < 0:
            raise ValueError('Negative error bound')
        loss = float(row['planned_gap_mm']) - float(row['measured_min_tool_gap_mm']) + err
        if not math.isfinite(loss):
            raise ValueError('Nonfinite measurement')
        groups.setdefault(row['guide'], []).append(loss)
    result = {}
    for (guide, loss) in groups.items():
        n = len(loss)
        q = max(loss)
        up = tail_upper(n, alpha, candidates)
        result[guide] = {'n_independent_blocks': n, 'observed_loss_upper_mm': q, 'tail_upper_at_max': up, 'tail_gate_pass': up < x, 'required_nominal_gap_mm': None if anatomy_budget is None or up >= x else 1 + q + anatomy_budget, 'anatomical_budget_supplied': anatomy_budget is not None, 'anatomical_budget_validated': False, 'anatomical_budget_contract': 'Requires external evidence; this CSV operator cannot validate the supplied B', 'clinical_injury_risk': 'UNKNOWN'}
    return {'status': 'CONDITIONAL_GEOMETRY_ONLY', 'guides': result, 'clinical_injury_risk': 'UNKNOWN'}
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', type=Path)
    ap.add_argument('--x', type=float, default=0.01)
    ap.add_argument('--confidence', type=float, default=0.95)
    ap.add_argument('--candidates', type=int, default=4)
    ap.add_argument('--anatomy-budget-mm', type=float)
    ap.add_argument('--output', type=Path)
    a = ap.parse_args()
    if not (0 < a.x < 1 and 0 < a.confidence < 1 and (a.candidates >= 1)):
        raise SystemExit('Invalid probability/confidence/candidate count')
    if a.anatomy_budget_mm is not None and a.anatomy_budget_mm < 0:
        raise SystemExit('Negative anatomical budget')
    rows = list(csv.DictReader(a.csv.open())) if a.csv else []
    out = calibrate(rows, a.x, 1 - a.confidence, a.candidates, a.anatomy_budget_mm)
    s = json.dumps(out, indent=2) + '\n'
    if a.output:
        a.output.write_text(s)
    print(s, end='')
