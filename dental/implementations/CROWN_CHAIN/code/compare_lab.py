import argparse, json, hashlib
from pathlib import Path
from contracts import compare_observation
R = Path(__file__).resolve().parents[1]

def main():
    p = argparse.ArgumentParser()
    p.add_argument('observations')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    f = R / 'FROZEN_PREDICTIONS.json'
    if hashlib.sha256(f.read_bytes()).hexdigest() != f.with_suffix('.sha256').read_text().split()[0]:
        raise ValueError('Prediction freeze drift')
    outpath = Path(a.output)
    if outpath.exists():
        raise ValueError('Use a new output; preserve earlier comparison')
    fr = json.loads(f.read_text())
    obs = json.loads(Path(a.observations).read_text())
    index = {(x['key'], x['quantity']): x for x in fr['laboratory_predictions']}
    rows = []
    for x in obs:
        f = index.get((x.get('key'), x.get('quantity')))
        rows.append(dict(key=x.get('key'), quantity=x.get('quantity'), comparison=compare_observation(f, x) if f else dict(status='REJECTED_UNKNOWN_SPECIMEN_OR_QUANTITY')))
    outpath.write_text(json.dumps(dict(rows=rows, status='UNKNOWN_NO_OBSERVATIONS' if not rows else 'SCORED_BEFORE_REFIT', scientific_admission=False), indent=2) + '\n')
if __name__ == '__main__':
    main()
