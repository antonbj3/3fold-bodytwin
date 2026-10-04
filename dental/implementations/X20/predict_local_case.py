"""Freeze the same geometric query on a NEW anonymous laboratory specimen.

NPZ keys: labels(z,y,x), spacing(z,y,x)mm. Label2=maxillary bone,
5/6=sinus,11..28=teeth,8/9/10=restorations. Do not reuse TF2 predictions
for another specimen. Clinically valid outcome and physical error stay UNKNOWN.
"""
import argparse
import datetime
import json
from pathlib import Path
import numpy as np
from measure import Local, trace, unit, dump, sha, scalar_decide
from capacity_map import candidate

def vec(s):
    return np.asarray([float(v) for v in s.split(',')])
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--labels', type=Path, required=True)
    p.add_argument('--case', required=True)
    p.add_argument('--fdi', type=int, required=True)
    p.add_argument('--crest', type=vec, required=True, help='Query frame mm in array z,y,x order')
    p.add_argument('--axis', type=vec, required=True, help='Inward depth direction')
    p.add_argument('--normal', type=vec, required=True, help='Approximate buccopalatal direction')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    assert not a.out.exists(), 'Never overwrite frozen predictions; use a new version path'
    data = np.load(a.labels)
    labels = data['labels']
    spacing = data['spacing']
    assert labels.ndim == 3 and len(spacing) == 3 and np.all(spacing > 0)
    axis = unit(a.axis)
    normal = unit(a.normal - np.dot(a.normal, axis) * axis)
    tangent = unit(np.cross(normal, axis))
    local = Local(labels, spacing, a.crest, axis)
    r = trace(local, a.crest, axis, normal)
    r['scalar_screen'] = scalar_decide(r)
    r['physical_or_clinical_decision'] = 'UNKNOWN'
    pred = []
    if r['valid']:
        for k in ('bone', 'sinus'):
            assert local.occ[k].any()
        local.signed('bone')
        local.distance('sinus')
        if local.distance('neighbor') is None:
            local.sdf['neighbor_dist'] = np.full(local.lab.shape, 99, dtype=np.float32)
        c = candidate(local, np.asarray(r['crest_mm']), axis, normal, tangent, 4.0, 6.0)
        r['main_candidate'] = c
        h = r['bone_height_mm']
        w = r['widths_mm']['1']
        pred = [{'case': a.case, 'fdi': a.fdi, 'crest_coordinate_mm_zyx': r['crest_mm'], 'axis_zyx': axis.tolist(), 'height_mm': h, 'height_interval_mm': [max(0, h - 0.6), h + 0.6], 'width1_mm': w, 'width1_interval_mm': [max(0, w - 0.6), w + 0.6], 'coordinate_validity': 'Independently selected experimental query frame required'}]
    out = {'claim_type': 'capability', 'frozen_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'physical_data_observed': False, 'input_labels_sha256': sha(a.labels), 'raw_geometry': r, 'predictions': pred, 'scope': 'Prototype failed validation gates on TF2; prospective lab test, not a calibrated physical prediction', 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://doi.org/10.1111/j.1708-8208.2008.00083.x', 'compared_quantity': 'CBCT vs independently registered physical height/width on the SAME specimen', 'refutes_us': True}}
    dump(a.out, out)
    a.out.with_suffix('.sha256').write_text(sha(a.out) + '\n')
    print('Frozen:', a.out, 'physical validity UNKNOWN')
