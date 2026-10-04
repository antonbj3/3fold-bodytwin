"""Research-only bite assessment. Local scans, explicit categorical sets and missing-landmark status."""
import argparse, json, pickle, time, zipfile, sys
from pathlib import Path
import numpy as np
from geometry import P, geometry, load_stl_bytes, proximity, clean
from benchmark import MF
STATUS = {'overbite': 'incisal-extreme proxy; numeric landmark accuracy UNKNOWN', 'overjet': 'incisal-extreme proxy; numeric landmark accuracy UNKNOWN', 'crossbite': 'radial outline proxy; cusp-by-tooth relationship UNKNOWN', 'midlines': 'notch proxy; interincisal landmark identity UNKNOWN', 'spee': 'profile/chord depth proxy; exact tooth endpoints UNKNOWN', 'molar_right': 'statistical text-category estimate; MB cusp / M1 buccal-groove identity UNKNOWN', 'molar_left': 'statistical text-category estimate; MB cusp / M1 buccal-groove identity UNKNOWN', 'canine_right': 'statistical text-category estimate; upper cusp / lower C-P1 embrasure identity UNKNOWN', 'canine_left': 'statistical text-category estimate; upper cusp / lower C-P1 embrasure identity UNKNOWN'}

def assess(u, l, case_id='local'):
    start = time.perf_counter()
    (features, g) = geometry(u, l, case_id)
    bundle = pickle.loads((P / 'raw/models.pkl').read_bytes())
    ths = json.load(open(P / 'raw/THRESHOLDS_R2.json'))
    X = np.array([[features.get(k, np.nan) for k in bundle['columns']]])
    out = {}
    for (f, md) in bundle['models'].items():
        if md.get('gb') is None:
            out[f] = dict(status='UNKNOWN_INSUFFICIENT_TRAIN_CLASSES')
            continue
        pr = md['gb'].predict_proba(X)[0]
        classes = md['classes']
        q = ths[f]['patient_max_q']
        S = [y for (y, p) in zip(classes, pr) if 1 - p <= q + 1e-12]
        out[f] = dict(category_set=S, decision='single_category' if len(S) == 1 else 'review_required', probabilities=dict(zip(classes, pr.tolist())), certainty='Marginal field calibration on local report judgments; joint whole-bite correctness not established (R4 FAIL); numeric landmark and scanner uncertainty UNKNOWN', measurement_kind=STATUS[f], supporting_measures={k: clean(features.get(k)) for k in MF[f]}, single_forced_baseline=classes[int(np.argmax(pr))])
    return dict(case_id=case_id, assessment_type='research categorical estimate; not clinical recommendation', fields=out, measures=clean(features), frame=clean(g.frame.__dict__), frame_limit='Axis-selected published frame: arbitrary rotations may invalidate values; see raw/RESULTS_R3.json. Laterality follows a convention without independent FDI verification in this archive. Vector frame diagnostics do not calibrate this model.', point_proximity=clean(proximity(g.upper_points, g.lower_points)), wall_seconds=time.perf_counter() - start)

def landmark_measures(spec):
    """Independent physical points open the identity bottleneck; categories need calibration."""

    def v(p):
        return np.asarray(p, dtype=float)

    def unit(p):
        a = v(p)
        return a / np.linalg.norm(a)
    s = unit(spec['superior_unit'])
    a = unit(spec['anterior_unit'])
    r = np.cross(a, s)
    if abs(np.dot(s, a)) > 1e-06:
        raise ValueError('Frame anterior/superior must be perpendicular')
    out = {'source': 'provided independent landmark measurements', 'categorical_Angle': 'UNKNOWN until external cutpoint/landmark agreement calibration', 'mm_error': 'UNKNOWN unless independently measured point error is provided'}
    if 'upper_inc_edge' in spec and 'lower_inc_edge' in spec:
        d = v(spec['upper_inc_edge']) - v(spec['lower_inc_edge'])
        out.update(overbite_mm=float(-d @ s), overjet_mm=float(d @ a))
    if 'upper_midline' in spec and 'lower_midline' in spec:
        out['midline_mm'] = float((v(spec['upper_midline']) - v(spec['lower_midline'])) @ r)
    for f in ('molar_right', 'molar_left', 'canine_right', 'canine_left'):
        if f in spec:
            x = spec[f]
            out[f + '_mesial_relation_mm'] = float((v(x['upper_cusp']) - v(x['lower_reference'])) @ unit(x['mesial_unit']))
    if 'spee' in spec:
        out['spee_depth_mm'] = {}
        for (side, x) in spec['spee'].items():
            A = v(x['canine_cusp'])
            B = v(x['distal_molar_cusp'])
            pts = v(x['posterior_cusps'])
            d = B - A
            dh = d - s * np.dot(d, s)
            t = (pts - A) @ dh / np.dot(dh, dh)
            chord = A + np.clip(t, 0, 1)[:, None] * d
            out['spee_depth_mm'][side] = float(np.max((chord - pts) @ s))
    if 'crossbite' in spec:
        out['crossbite_relations_mm'] = [dict(tooth=x['tooth'], upper_minus_lower_buccal_mm=float((v(x['upper_buccal_cusp']) - v(x['lower_buccal_cusp'])) @ unit(x['buccal_unit']))) for x in spec['crossbite']]
    return out
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case')
    ap.add_argument('--upper', type=Path)
    ap.add_argument('--lower', type=Path)
    ap.add_argument('--landmarks', type=Path)
    ap.add_argument('--out', type=Path)
    ar = ap.parse_args()
    if ar.landmarks:
        res = landmark_measures(json.loads(ar.landmarks.read_text()))
    elif ar.case:
        mf = json.load(open(P / 'raw/DATA_MANIFEST.json'))
        with zipfile.ZipFile(mf['zip']) as z:
            (u, _) = load_stl_bytes(z.read(ar.case + '/ios/ios_upper.stl'))
            (l, _) = load_stl_bytes(z.read(ar.case + '/ios/ios_lower.stl'))
        res = assess(u, l, ar.case)
    elif ar.upper and ar.lower:
        (u, _) = load_stl_bytes(ar.upper.read_bytes())
        (l, _) = load_stl_bytes(ar.lower.read_bytes())
        res = assess(u, l)
    else:
        ap.error('Provide --case OR --upper + --lower OR --landmarks')
    txt = json.dumps(clean(res), indent=2, allow_nan=False) + '\n'
    if ar.out:
        ar.out.write_text(txt)
    else:
        print(txt)
