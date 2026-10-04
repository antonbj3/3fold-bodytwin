import argparse
import json
import numpy as np
from common import DATA, write
from round1 import score_rows

def main():
    p = argparse.ArgumentParser(description='Score any candidate-ROI detector against private voxel-derived depth truth')
    p.add_argument('--predictions', required=True)
    p.add_argument('--truth', default=str(DATA / 'private/truth_index.json'))
    p.add_argument('--threshold', type=float, required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    rows = json.load(open(a.truth))
    pred = json.load(open(a.predictions))
    if not np.isfinite(a.threshold):
        raise ValueError('Non-finite threshold')
    ids = [r['sample_id'] for r in pred]
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate prediction id')
    mapping = {r['sample_id']: float(r['probability']) for r in pred}
    truth_ids = {r['sample_id'] for r in rows}
    if set(mapping) != truth_ids:
        raise ValueError('Missing or unknown sample ids; provide all split predictions')
    if any((not np.isfinite(v) or not 0 <= v <= 1 for v in mapping.values())):
        raise ValueError('Invalid probability')
    write(a.output, score_rows(rows, np.array([mapping[r['sample_id']] for r in rows]), a.threshold))
if __name__ == '__main__':
    main()
