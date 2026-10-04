"""Compare later independent point measurements with already-frozen candidates."""
import json, argparse
from pathlib import Path
import numpy as np
from contact import P, sha, write

def prepare():
    m = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    op = json.load(open(P / 'raw/OPPOSITION_FEATURES_R2.json'))
    records = []
    for c in m['numerical_panel']:
        from contact import LABELS
        meta = json.load(open(LABELS / c / 'labels+landmarks.json'))
        for jaw in ['upper', 'lower']:
            for (k, t) in meta['arches'][jaw]['landmarks']['teeth'].items():
                if int(k) % 10 not in [1, 6]:
                    continue
                roles = {'incisal': t['incisal_candidate_xyz_mm']} if int(k) % 10 == 1 else {'mesiobuccal_cusp': t['mesiobuccal_cusp_candidate_xyz_mm'], 'buccal_groove': t['buccal_groove_candidate_xyz_mm']}
                for (role, xyz) in roles.items():
                    if xyz is not None:
                        records.append(dict(case_id=c, jaw=jaw, fdi=int(k), role=role, xyz_mm=xyz, source_sha256=meta['arches'][jaw]['input_sha256']))
    write(P / 'raw/MEASUREMENT_PREDICTIONS.json', dict(records=records))
    write(P / 'FROZEN_MEASUREMENT_PANEL.json', dict(frozen_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(), before_any_new_point_measurement=True, predictions_sha256=sha(P / 'raw/MEASUREMENT_PREDICTIONS.json'), cases=m['numerical_panel'], maximum_point_error_mm=0.5, current_external_measurement_status='NOT_RUN', code_sha256=sha(Path(__file__))))

def compare(path, out):
    freeze = json.load(open(P / 'FROZEN_MEASUREMENT_PANEL.json'))
    assert sha(P / 'raw/MEASUREMENT_PREDICTIONS.json') == freeze['predictions_sha256']
    pred = json.load(open(P / 'raw/MEASUREMENT_PREDICTIONS.json'))
    obs = json.load(open(path))
    key = lambda r: (r['case_id'], r['jaw'], int(r['fdi']), r['role'])
    lut = {key(r): r for r in pred['records']}
    rr = []
    for r in obs['records']:
        k = key(r)
        if k not in lut:
            raise ValueError('UNBOUND_POINT ' + str(k))
        p = lut[k]
        if r['source_sha256'] != p['source_sha256']:
            raise ValueError('WRONG_SOURCE ' + str(k))
        e = float(np.linalg.norm(np.array(r['xyz_mm']) - np.array(p['xyz_mm'])))
        rr.append(dict(key=list(k), error_mm=e, gate='PASS' if e <= 0.5 else 'FAIL'))
    write(out, dict(records=rr, observation_sha256=sha(path), frozen_prediction_sha256=freeze['predictions_sha256'], claim_type='capability', gate='UNKNOWN_NO_MEASUREMENTS' if not rr else 'PASS' if all((x['gate'] == 'PASS' for x in rr)) else 'FAIL', unmeasured_frozen_points=len(lut) - len({tuple(x['key']) for x in rr}), refitting=False, measurement_independence='Must be documented by observation provider; software cannot establish independence.'))
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('measurement', nargs='?', type=Path)
    ap.add_argument('--prepare', action='store_true')
    ap.add_argument('--out', type=Path, default=P / 'raw/MEASUREMENT_COMPARISON.json')
    a = ap.parse_args()
    prepare() if a.prepare else compare(a.measurement, a.out)
