"""Read one completed case, or compute a conditional map from your labelled pair.

For an external pair, labels JSON must be X11 labels+landmarks.json with source
SHA256, unique-vertex indexing and labels.path. Scanner/label accuracy stays UNKNOWN.
"""
import argparse, json, datetime, time
from pathlib import Path
import numpy as np
from contact import P, D, sha, write, prepare, boxes, buckets, extrema, make_map, summaries

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case', help='Existing anonymous Bite2Text case ID')
    ap.add_argument('--upper', type=Path)
    ap.add_argument('--lower', type=Path)
    ap.add_argument('--segmentation', type=Path)
    ap.add_argument('--out', type=Path, default=P / 'raw/DEMO_CASE.json')
    a = ap.parse_args()
    if a.case:
        row = json.load(open(P / 'raw/cases' / (a.case + '.json')))
        tight = P / 'raw/refined_certificates' / (a.case + '.json')
        cert = json.load(open(tight if tight.exists() else P / 'raw/certificates' / (a.case + '.json')))
        opp = json.load(open(P / 'raw/OPPOSITION_FEATURES_R2.json'))[a.case]
        write(a.out, dict(case_id=a.case, geometry=row, minimum_gap_certificates=cert, opposition_candidate=opp, physical_contact='UNKNOWN', clinical_point_accuracy='UNKNOWN'))
        print(a.out)
        return
    if not all([a.upper, a.lower, a.segmentation]):
        ap.error('Provide --case, or --upper --lower --segmentation together')
    start = time.perf_counter()
    meta = json.load(open(a.segmentation))
    fr = meta['arches']['lower']['frame']
    R = np.column_stack([fr['right_unit'], fr['anterior_unit'], fr['superior_unit']])
    center = np.array(fr['center_mm'])
    if not np.allclose(R.T @ R, np.eye(3), atol=1e-10):
        raise ValueError('Nonorthonormal frame')
    arches = [prepare(path.read_bytes(), meta['arches'][jaw], jaw, R, center) for (path, jaw) in [(a.upper, 'upper'), (a.lower, 'lower')]]
    (U, ul, uc, uf, us) = arches[0]
    (L, ll, lc, lf, ls) = arches[1]
    ub = boxes(U)
    lb = boxes(L)
    origin = np.minimum(ub[:, :2].min(0), lb[:, :2].min(0)) - 0.8
    hi = np.maximum(ub[:, 3:5].max(0), lb[:, 3:5].max(0))
    shape = np.ceil((hi - origin) / 0.8).astype(np.int64) + 2
    (ptr, ids) = buckets(lb, origin, 0.8, shape)
    (best, wi, xy, n, ne) = extrema(U, L, ul, ll, uc, lc, ub, lb, origin, 0.8, shape, ptr, ids, True)
    m = make_map(U, L, ul, ll, uc, lc)
    (pairs, features) = summaries(m, best, wi, xy, U, L, uc, lc, uf, lf)
    root = D / 'demos' / datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%f')
    root.mkdir(parents=True)
    mp = root / 'map.npz'
    np.savez_compressed(mp, **m, minimum_gap=best)
    write(a.out, dict(units='mm', model='Rounded labelled fixed-pose projected geometry; floating minima, no new exact certificate in this single-pair path', pairs=pairs, features=features, map_path=str(mp), map_sha256=sha(mp), frame=dict(center=center, R=R), stats=dict(upper=us, lower=ls), wall_s=time.perf_counter() - start, physical_contact='UNKNOWN', clinical_diagnosis='NOT_PREDICTED_FOR_EXTERNAL_PAIR', segmentation_accuracy='UNKNOWN'))
    print(a.out)
if __name__ == '__main__':
    main()
