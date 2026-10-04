"""Observable disagreement diagnostics. Do not infer causes from category errors."""
from dental_release.paths import expand as _release_expand
import json
import numpy as np
from benchmark import P, write
from parser import FIELDS
from diagnose import STATUS
from benchmark import MF

def run():
    ys = json.load(open(P / 'raw/TEST_LABELS.json'))
    pred = {r['case_id']: r for r in json.load(open(P / 'raw/PREDICTIONS_R1.json')) if r['split'] == 'test'}
    g = {r['case_id']: r for r in map(json.loads, (P / 'raw/geometry.jsonl').read_text().splitlines())}
    rows = []
    summ = {}
    groups = {}
    for f in FIELDS:
        rr = []
        gg = {}
        for (c, reps) in ys.items():
            if not reps or f not in pred[c]['fields']:
                continue
            first = reps[0]['labels'][f]
            point = pred[c]['fields'][f]['point']
            allref = [r['labels'][f] for r in reps if r['labels'][f] is not None]
            if first is None:
                continue
            feats = {k: g[c]['features'].get(k) for k in MF[f]}
            gg.setdefault(first, []).append(feats)
            if point != first:
                secondary = point in allref[1:]
                disc = len(set(allref)) > 1
                missing = [k for (k, v) in feats.items() if v is None]
                rationale = 'Prediction agrees with another report of the same measured case; source judgments differ.' if secondary else 'Prediction differs from every readable report; numeric landmark error vs categorical model error vs parser error cannot be separated without independent tooth landmarks.'
                if missing:
                    rationale += ' Field-specific proxy(s) missing; model used other correlated features.'
                r = dict(case_id=c, field=f, first_reference=first, model_category=point, all_readable_reports=allref, patient_report_discordance=disc, matches_other_report=secondary, missing_primitive_measures=missing, supporting_geometry=feats, observable_explanation=rationale, landmark_status=STATUS[f], causal_explanation='UNKNOWN')
                rows.append(r)
                rr.append(r)
        summ[f] = dict(errors=len(rr), matches_another_report=sum((r['matches_other_report'] for r in rr)), discordant_reference_patient=sum((r['patient_report_discordance'] for r in rr)), missing_primitive_measure=sum((bool(r['missing_primitive_measures']) for r in rr)))
        groups[f] = {}
        for (cat, vals) in gg.items():
            groups[f][cat] = {k: dict(n=sum((v[k] is not None for v in vals)), median=float(np.median([v[k] for v in vals if v[k] is not None])) if any((v[k] is not None for v in vals)) else None) for k in MF[f]}
    write('raw/ERROR_ANALYSIS.json', dict(summary=summ, rows=rows, categorical_numeric_medians=groups, scope='Observed first-reference test errors; explanations distinguish measured source discordance from unresolved causes. No clinician was adjudicated incorrect.'))
    write('OCCLUSION_PORT.json', dict(status='PENDING_INDEPENDENT_REVIEW', source='Bite2Text registered upper/lower STL', case_records='raw/geometry.jsonl', coordinate_frame='per-case canonical rotation matrix stored in record.frame.rotation; coordinates transform with rotation; no relative jaw repositioning', units='mm, degrees', point_proximity='record.contact: paired canonical xyz and sampled nearest vertex Euclidean distance', surface_distance_validity='Each sampled point distance upper-bounds global full-surface minimum; sampling resolution error UNKNOWN. Far points do not exclude near triangle contact.', FDI_identity='UNKNOWN', load_N=None, force_distribution=None, consumers=[_release_expand('X2'), 'LANE_NEXT_P_OCCLUSION_VALIDATION R6'], next_operation='Use original full triangles plus independent FDI cusp/groove landmarks; no forced-load or support values inferred.'))
    print('Error analysis', summ)
if __name__ == '__main__':
    run()
