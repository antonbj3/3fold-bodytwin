import json, pathlib, datetime, hashlib, time
import numpy as np
from sweep import *
R = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def freeze():
    keys = [x['key'] for x in json.load(open(R / 'rounds/R1/COHORT.json'))['rows'] if x['full_path_status'] == 'COLLISION']
    axes = sorted([[float(x), float(y), 40.0] for x in [-4, -2, 0, 2, 4] for y in [-4, -2, 0, 2, 4]], key=lambda d: (d[0] ** 2 + d[1] ** 2, d[0], d[1]))
    p = R / 'PREREG_R4.json'
    p.write_text(json.dumps(dict(round='R4', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), claim_type='capability', capability='Find the least-tilted surface-clear straight axis in a declared finite class while preserving complete crown and margin', obstacle='R2inside-only repair is mathematically unable to change fixed exterior collision', changed_operation='Change translation direction, keep rigid orientation and all original mesh coordinates; sweep every triangle pair for each axis', keys=keys, axes_mm=axes, path='final+sD,sin[0,1],no rotation', metric='Whole-surface clearance rigorously certified or collision witness or UNKNOWN for each tested direction;25declared candidate axes', selection='Ordered by squared transverse component then lexicographic; stop at first fully certified neighbour surface-clear direction. Minimality in this discrete class only if all earlier directions have exact COLLISION decisions.', budget_per_surface_s=20, gate={'exact_mesh_and_margin_identity_error': 0, 'unsupported_scene_abstain_fraction': 1, 'control_boolean_agreement': 'exact'}, strongest_equally_informed_control='Same barycentric LP and rational witnesses; conservative separating enclosure required for absence', falsifier='Selected clear-axis record mutated to nominal axis must be rejected by a known source-bound collision witness', full_cost={'fit': 0, 'preparation': 'R1frozen meshes', 'discovery': 'up to25axes per5blocked sites, nominal reused exactly', 'validation': 'every swept facet pair or exact collision', 'query': 'source-neighbour surfaces only', 'fallback': 'UNKNOWN if any uncovered candidate', 'physical': 'UNKNOWN actual preparation, registration and solid scene'}, resolution='PER_POINT', timescale='SIMULTANEOUS'), indent=2) + '\n')
    p.with_suffix('.sha256').write_text(sha(p) + '\n')
    print('R4FROZEN25axes5sites')

def run(outdir):
    tick = time.perf_counter()
    pre = json.load(open(R / 'PREREG_R4.json'))
    co = json.load(open(R / 'FROZEN_COHORT.json'))
    r1 = json.load(open(R / 'rounds/R1/COHORT.json'))
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for key in pre['keys']:
        rec = next((q for q in co['rows'] if q['key'] == key))
        original = next((q for q in r1['rows'] if q['key'] == key))
        trials = []
        selected = None
        with np.load(rec['mesh_path']) as m, np.load(rec['public_path']) as p:
            tri = m['vertices'][m['faces']]
            roles = m['roles']
            obs = {side: Obstacle(p[side]) for side in ['mesial', 'distal']}
            for D0 in pre['axes_mm']:
                D = np.array(D0)
                results = {}
                if D0 == [0.0, 0.0, 40.0]:
                    status = 'COLLISION'
                    results = original['neighbours']
                    reused = True
                else:
                    reused = False
                    status = 'SURFACE_CLEAR_CERTIFIED'
                    for side in ['mesial', 'distal']:
                        ow = original['neighbours'][side]['path'].get('witness')
                        quick = None
                        if ow:
                            i = ow['source_face_id']
                            j = ow['obstacle_face_id']
                            (w, st) = lp_witness(tri[i], p[side][j], D)
                            if w:
                                w.update(source_face_id=i, obstacle_face_id=j, source_role=int(roles[i]), source_triangle_mm=tri[i].tolist(), obstacle_triangle_mm=p[side][j].tolist(), translation_mm=D0, independent_control=independent_lp(tri[i], p[side][j], D))
                                quick = dict(status='COLLISION', witness=w, fast_reused_face=True)
                        q = quick or check_surface(tri, obs[side], D, roles=roles, budget_s=pre['budget_per_surface_s'])
                        results[side] = q
                        if q['status'] == 'COLLISION':
                            status = 'COLLISION'
                            break
                        if q['status'] != 'SURFACE_CLEAR_CERTIFIED':
                            status = 'UNKNOWN'
                trial = dict(translation_mm=D0, status=status, results=results, reused_nominal=reused)
                trials.append(trial)
                if status == 'SURFACE_CLEAR_CERTIFIED':
                    selected = trial
                    break
            before = trials[:-1] if selected else trials
            row = dict(key=key, status='NEIGHBOUR_SURFACE_AXIS_FOUND' if selected else 'NO_CERTIFIED_AXIS_IN_TESTED_CLASS', selected_translation_mm=selected['translation_mm'] if selected else None, selected_axis_tilt_degrees=float(np.degrees(np.arctan(np.linalg.norm(selected['translation_mm'][:2]) / 40))) if selected else None, finite_class_minimum_certified=bool(selected and all((x['status'] == 'COLLISION' for x in before))), tested_axes=len(trials), trials=trials, crown_identity_error=0.0, margin_identity_error_mm=0.0, full_prep_scene_path='UNKNOWN; surface-only straight-axis change may collide with virtual/actual preparation', physical_insertion='UNKNOWN', resolution='PER_POINT')
            rows.append(row)
            (outdir / 'AXES.json').write_text(json.dumps(dict(rows=rows, seconds=time.perf_counter() - tick), indent=2) + '\n')
            print(key, row['status'], row['selected_translation_mm'], len(trials), round(time.perf_counter() - tick, 2), flush=True)
    return rows
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--freeze', action='store_true')
    ap.add_argument('--out', default=str(R / 'rounds/R4'))
    a = ap.parse_args()
    freeze() if a.freeze else run(pathlib.Path(a.out))
