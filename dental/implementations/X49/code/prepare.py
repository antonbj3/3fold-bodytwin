"""Prepare bounded inputs and freeze expected outputs BEFORE gate evaluation."""
from dental_release.paths import expand as _release_expand
import datetime as dt
import json
import sys
from pathlib import Path
import numpy as np
import trimesh
from parents import HERE, RESULTS, sha, transport as T
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X49_design_gate'))
OLD = Path(_release_expand('@DENTAL_WORK_ROOT@/X34_stl_design_gate'))

def save_case(name, tri, base, changes=None):
    out = DATA / name
    out.mkdir(parents=True, exist_ok=True)
    meta = dict(base)
    meta.update(changes or {})
    for (k, t) in tri.items():
        T.write_stl(out / (k + '.stl'), t)
    meta['input_sha256'] = {k: sha(out / (k + '.stl')) for k in tri}
    T.dump(out / 'contract.json', meta)
    return out

def cup():
    ring = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)
    vertices = np.vstack([np.c_[2 * ring, np.zeros(4)], np.c_[2 * ring, np.full(4, 3.0)], np.c_[ring, np.zeros(4)], np.c_[0.75 * ring, np.full(4, 2.5)]])
    faces = []
    labels = []

    def add(f, label):
        faces.append(f)
        labels.append(label)
    for i in range(4):
        j = (i + 1) % 4
        for f in ([i, j, 4 + j], [i, 4 + j, 4 + i]):
            add(f, 'exterior')
        for f in ([8 + i, 12 + j, 8 + j], [8 + i, 12 + i, 12 + j]):
            add(f, 'intaglio')
        for f in ([i, 8 + j, j], [i, 8 + i, 8 + j]):
            add(f, 'margin')
    for f in ([4, 5, 6], [4, 6, 7]):
        add(f, 'occlusal')
    for f in ([12, 14, 13], [12, 15, 14]):
        add(f, 'intaglio')
    tri = vertices[np.array(faces)]
    prep = trimesh.creation.box(extents=[1.0, 1.0, 1.5])
    prep.apply_translation([0, 0, 0.75])
    a = np.array([[[-3, -3, 3.025], [3, -3, 3.025], [3, 3, 3.025]], [[-3, -3, 3.025], [3, 3, 3.025], [-3, 3, 3.025]]])
    regions = {k: [i for (i, l) in enumerate(labels) if l == k] for k in ('intaglio', 'exterior', 'occlusal', 'margin')}
    regions['exterior'] += regions['occlusal']
    return (dict(crown=tri, prep=prep.triangles, antagonist=a), regions)

def main():
    if (HERE / 'FROZEN_PREDICTIONS.json').exists():
        raise ValueError('Already frozen; reuse inputs or preregister a new round')
    DATA.mkdir(parents=True, exist_ok=True)
    cases = {}
    for name in ('valid', 'thin', 'film_low', 'film_high', 'penetration'):
        old = OLD / name
        tri = {k: T.load_mesh(old / (k + '.stl'), 'mm')[0] for k in ('prep', 'crown', 'antagonist')}
        meta = json.loads((old / 'contract_R2.json').read_text())
        meta['surface_error_mm'] = {k: 0.001 for k in tri}
        meta['error_source'] = 'Synthetic supplied envelope; NOT measured scanner calibration'
        cases[name] = str(save_case(name, tri, meta))
    tri = {k: T.load_mesh(DATA / 'valid' / (k + '.stl'), 'mm')[0] for k in ('prep', 'crown', 'antagonist')}
    meta = json.loads((DATA / 'valid/contract.json').read_text())
    cases['near_limit'] = str(save_case('near_limit', tri, meta, dict(surface_error_mm={k: 0.1 for k in tri})))
    cases['wrong_scale'] = str(save_case('wrong_scale', {k: 1000 * t for (k, t) in tri.items()}, meta))
    cases['explicit_um'] = str(save_case('explicit_um', {k: 1000 * t for (k, t) in tri.items()}, meta, dict(units='um')))
    cases['open_crown'] = str(save_case('open_crown', dict(tri, crown=tri['crown'][:-1]), meta))
    cases['reverse_winding'] = str(save_case('reverse_winding', dict(tri, crown=tri['crown'][:, ::-1]), meta))
    out = save_case('reverse_stored_normal', tri, meta)
    b = bytearray((out / 'crown.stl').read_bytes())
    dtype = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
    rows = np.frombuffer(b, dtype=dtype, count=len(tri['crown']), offset=84)
    rows['n'][0] *= -1
    (out / 'crown.stl').write_bytes(b)
    meta['input_sha256'] = {k: sha(out / (k + '.stl')) for k in tri}
    T.dump(out / 'contract.json', meta)
    cases['reverse_stored_normal'] = str(out)
    old = OLD / 'stepped_cavity_closed'
    tri = {k: T.load_mesh(old / (k + '.stl'), 'mm')[0] for k in ('prep', 'crown', 'antagonist')}
    meta = json.loads((old / 'contract_R2.json').read_text())
    cases['undercut'] = str(save_case('undercut', tri, meta))
    (tri, regions) = cup()
    meta = dict(units='mm', common_frame_confirmed=True, input_sha256={}, regions=regions, contact_axis=[0, 0, 1], convex_cavity_complete=True, extraction_axis=[0, 0, 1], extraction_travel_mm=4.0, crown_solid_integrity_asserted=True, prep_solid_integrity_asserted=True, bur_radii_final_mm=[0.3], milling_allowance_mm=0.3, ifu_profile='katana-ht', indication='posterior', occlusal_gap_limits_mm=[-0.03, 0.1], region_semantics='Analytical square tapered cup, complete cavity; synthetic')
    cases['convex_cup'] = str(save_case('convex_cup', tri, meta))
    cases['oversized_bur'] = str(save_case('oversized_bur', tri, meta, dict(bur_radii_final_mm=[2.0])))
    a = DATA / 'X1b_missing_antagonist_surrogate.stl'
    T.write_stl(a, tri['antagonist'] + np.array([0, 0, 100.0]))
    x1b = []
    for name in ('D1', 'D2', 'D3', 'M1', 'M2'):
        p = RESULTS / 'LANE_X1B_CROWN_LOOP/exports' / name
        x1b.append(dict(case=name, paths={k: str(v) for (k, v) in dict(prep=p / 'preparation.stl', crown=p / 'crown.stl', antagonist=a).items()}, registered_antagonist='MISSING; surrogate only exercises ingestion; no occlusal validation'))
    expected = {'valid': {'mesh_health': 'PASS', 'scale': 'PASS', 'material_wall': 'PASS', 'film_min': 'PASS', 'film_max': 'PASS'}, 'thin': {'material_wall': 'FAIL'}, 'film_low': {'film_min': 'FAIL'}, 'film_high': {'film_max': 'FAIL'}, 'penetration': {'occlusal_contact': 'UNKNOWN'}, 'near_limit': {'material_wall': 'PASS', 'film_min': 'UNKNOWN', 'film_max': 'UNKNOWN'}, 'wrong_scale': {'scale': 'FAIL'}, 'explicit_um': {'scale': 'PASS', 'material_wall': 'PASS', 'film_min': 'PASS'}, 'open_crown': {'mesh_health': 'FAIL'}, 'reverse_winding': {'mesh_health': 'FAIL'}, 'reverse_stored_normal': {'mesh_health': 'FAIL'}, 'undercut': {'insertion': 'FAIL', 'ball_milling': 'FAIL'}, 'convex_cup': {'mesh_health': 'PASS', 'insertion': 'PASS', 'ball_milling': 'PASS'}, 'oversized_bur': {'ball_milling': 'FAIL'}}
    for name in ('projected_valid', 'projected_penetration'):
        src = DATA / ('valid' if name == 'projected_valid' else 'penetration')
        tri = {k: T.load_mesh(src / (k + '.stl'), 'mm')[0] for k in ('prep', 'crown', 'antagonist')}
        meta = json.loads((src / 'contract.json').read_text())
        meta['surface_error_mm'] = {k: 0.0 for k in tri}
        cases[name] = str(save_case(name, tri, meta))
        expected[name] = {'occlusal_contact': 'PASS' if name == 'projected_valid' else 'FAIL'}
    freeze = dict(round='R1', claim_type='capability', frozen_utc=dt.datetime.now(dt.timezone.utc).isoformat(), prereg_sha256=sha(HERE / 'PREREG_R1.json'), expected=expected, cases=cases, x1b=x1b, externals='DEPENDENCIES.json and PREREG_R1.json', sufficiency=dict(a=[0.625, 0.875], b=[0.375, 1.125], summary='equal-area mean thickness', expected_summary_error=0.0, expected_min_difference_mm=0.25), prospective_lab_measurements=False, physical_measurement_status='NONE', code={str(p.relative_to(HERE)): sha(p) for p in sorted((HERE / 'code').glob('*.py'))}, files=[dict(path=str(p), sha256=sha(p)) for p in sorted(DATA.rglob('*')) if p.is_file()] + [dict(path=p, sha256=sha(p)) for row in x1b for (k, p) in row['paths'].items() if k != 'antagonist'])
    T.dump(HERE / 'FROZEN_PREDICTIONS.json', freeze)
    (HERE / 'FROZEN_PREDICTIONS.sha256').write_text(sha(HERE / 'FROZEN_PREDICTIONS.json') + '\n')
    T.dump(HERE / 'CURRENT_WORK_STATE.json', dict(lane='X49-design-gate', status='R1_PREDICTIONS_FROZEN', latest_gate='No gate query run yet', next_operation='Run fixed fault matrix, independent controls and X1b ingestion', updated_utc=freeze['frozen_utc']))
    print('Frozen', len(cases), 'cases;', len(freeze['files']), 'inputs')
if __name__ == '__main__':
    main()
