from dental_release.paths import expand as _release_expand
from common import *
from data import load_geometry
from experiment import matrix
import argparse, pickle

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default=_release_expand('@DENTAL_CASE_ID@'))
    ap.add_argument('--out', default='raw/DEMO_CASE.json')
    args = ap.parse_args()
    g = load_geometry()
    assert args.case in g, 'Unknown case ID'
    rows = [dict(case_id=args.case, fdi=f) for f in sorted(g[args.case])]
    r1 = read('FROZEN_PREDICTIONS_R1.json')
    assert sha(r1['model_path']) == r1['model_sha256']
    model = pickle.loads(Path(r1['model_path']).read_bytes())
    scores = model['models']['shape'].predict_proba(matrix(rows, g, model['columns']))[:, 1]
    r3 = read('FROZEN_PREDICTIONS_R3.json')
    assert sha(r3['model_path']) == r3['model_sha256']
    absence = pickle.loads(Path(r3['model_path']).read_bytes())
    a_scores = absence['models']['shape'].predict_proba(matrix(rows, g, absence['columns']))[:, 1]
    ports = {(r['case_id'], r['fdi']): r for r in read('raw/MATERIAL_APPLICABILITY_PORT_R5.json')}
    output = []
    for (r, s, a) in zip(rows, scores, a_scores):
        output.append(dict(**r, shape_restoration_score=float(s), shape_restoration_flag=bool(s >= model['thresholds']['shape']), shape_absence_score=float(a), shape_absence_flag=bool(a >= absence['thresholds']['shape']), observed_material_state=ports.get((r['case_id'], r['fdi']))))
    write(args.out, dict(case_id=args.case, rows=output, geometry_source=next((r for r in read('raw/GEOMETRY_INPUT_MANIFEST.json') if r['case_id'] == args.case)), interpretation='Research scores, not calibrated diagnosis probabilities. Material law and restoration boundary/volume UNKNOWN. Report-linked state is the exported information.', resolution='PER_TOOTH', timescale='SIMULTANEOUS'))
    print('Case demo', args.case, 'saved to', args.out, flush=True)
if __name__ == '__main__':
    main()
