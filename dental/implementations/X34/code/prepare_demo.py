"""Prepare tiny inputs and freeze expected decisions before running the gate."""
from dental_release.paths import expand as _release_expand
import datetime as dt
import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from designgate.vendor.geometry import export_stl, shell, parse_stl
from initialize import dump, sha, H
D = Path(_release_expand('@DENTAL_WORK_ROOT@/X34_stl_design_gate'))

def contract(case, paths, regions, **kw):
    x = dict(units='mm', common_frame_confirmed=True, contact_axis=[0, 0, 1], input_sha256={k: sha(p) for (k, p) in paths.items()}, ifu_profile='katana-ht', indication='posterior', film_limits_mm=[0.02, 0.12], film_limits_source='Declared demonstration lab protocol; not universal clinical limits', occlusal_gap_limits_mm=[-0.03, 0.1], geometry_scope='matched_height_graph_patches', regions=regions, physical_geometry_status='closed-form fixture', **kw)
    dump(D / case / 'contract.json', x)

def main():
    if (H / 'FROZEN_PREDICTIONS.json').exists():
        raise RuntimeError('Refusing to refreeze')
    D.mkdir(parents=True, exist_ok=True)
    xy = np.array([[-2, -2], [2, -2], [-2, 2], [2, 2]], float)
    f = np.array([[0, 1, 2], [1, 3, 2]])
    for (name, wall, film, gap) in [('valid', 0.8, 0.05, 0.025), ('thin', 0.2, 0.05, 0.025), ('film_low', 0.8, 0.005, 0.025), ('film_high', 0.8, 0.2, 0.025), ('penetration', 0.8, 0.05, -0.2)]:
        out = D / name
        out.mkdir(exist_ok=True)
        (v, fs) = shell(xy, np.full(4, film + wall), np.full(4, film), f)
        paths = {k: out / (k + '.stl') for k in ('prep', 'crown', 'antagonist')}
        export_stl(paths['crown'], v, fs)
        export_stl(paths['prep'], np.c_[xy, np.zeros(4)], f)
        export_stl(paths['antagonist'], np.c_[xy, np.full(4, film + wall + gap)], f)
        contract(name, paths, dict(intaglio=[2, 3], exterior=[0, 1], occlusal=[0, 1], preparation=[0, 1], antagonist=[0, 1]))
    out = D / 'real_X18_110_34'
    out.mkdir(exist_ok=True)
    src = Path(_release_expand('@DENTAL_WORK_ROOT@/X18_crown_antagonist/110_34_R2_informed.stl'))
    shutil.copyfile(src, out / 'crown.stl')
    (v, fs) = parse_stl(src.read_bytes())
    tri = v[fs]
    n = 800
    plane_z = float(tri[n:2 * n, :, 2].min() - 0.05)
    prep = tri[:n].copy()
    prep[:, :, 2] = plane_z
    export_stl(out / 'prep.stl', prep.reshape(-1, 3), np.arange(3 * n).reshape(-1, 3))
    zpath = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/Bits2Bites/Bits2Bites_v01.zip'))
    with zipfile.ZipFile(zpath) as z:
        member = next((x for x in z.namelist() if x.endswith('/110/upper.stl')))
        blob = z.read(member)
    (uv, uf) = parse_stl(blob)
    ut = uv[uf]
    lo = tri[:n, :, :2].min((0, 1))
    hi = tri[:n, :, :2].max((0, 1))
    keep = np.all(ut[:, :, :2].max(1) >= lo, axis=1) & np.all(ut[:, :, :2].min(1) <= hi, axis=1)
    det = np.cross(ut[:, 1, :2] - ut[:, 0, :2], ut[:, 2, :2] - ut[:, 0, :2])
    keep &= np.abs(det) > 1e-12
    selected = ut[keep]
    export_stl(out / 'antagonist.stl', selected.reshape(-1, 3), np.arange(3 * len(selected)).reshape(-1, 3))
    paths = {k: out / (k + '.stl') for k in ('prep', 'crown', 'antagonist')}
    contract(out.name, paths, dict(intaglio=list(range(n, 2 * n)), exterior=list(range(n)), occlusal=list(range(n)), preparation=list(range(n)), antagonist=list(range(len(selected)))))
    meta = json.loads((out / 'contract.json').read_text())
    meta['physical_geometry_status'] = 'X18 conditional generated roof shell; preparation is a NEW VIRTUAL HORIZONTAL PLANE, not a patient preparation; antagonist facets from published registered scan'
    dump(out / 'contract.json', meta)
    manifest = dict(dataset='Bits2Bites', license='CC BY-NC-SA 4.0; use local research copy, attribute dataset authors', locator='https://ditto.ing.unimore.it/bits2bites/', source_zip=str(zpath), source_member=member, source_member_sha256=hashlib.sha256(blob).hexdigest(), source_crown=str(src), source_crown_sha256=sha(src), compared_quantity='Signed projected original antagonist facet to X18 generated crown exterior separation', selection='XY bounding-box overlap (conservative) with crown; nondegenerate projections', selected_source_face_ids=np.flatnonzero(keep).tolist(), total_source_faces=len(ut), retained=len(selected), rejected=len(ut) - len(selected), rejection_reason='Outside crown projected query domain or zero projected area; no outcome-dependent selection', files={k: dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size) for (k, p) in paths.items()})
    dump(H / 'raw/REAL_DATA_MANIFEST.json', manifest)
    files = []
    for p in sorted(D.glob('*/*')):
        if p.is_file():
            files.append(dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size))
    predictions = dict(round='R1', claim_type='capability', frozen_utc=dt.datetime.now(dt.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R1.json'), expected=dict(valid={'material_wall': 'PASS', 'film_min': 'PASS', 'film_max': 'PASS', 'occlusal_contact': 'PASS', 'insertion': 'PASS', 'ball_milling': 'PASS'}, thin={'material_wall': 'FAIL'}, film_low={'film_min': 'FAIL'}, film_high={'film_max': 'FAIL'}, penetration={'occlusal_contact': 'FAIL'}), analytical=dict(valid=dict(minimum_wall_mm=0.8, film_min_mm=0.05, film_max_mm=0.05, minimum_occlusal_gap_mm=0.025), crossing_triangles_minimum_gap_mm=-0.1), real_prediction=dict(occlusal_gap_mm=0.025, tolerance_mm=1e-05, provenance='Inherited X18 R2 repair target already known; retrospective verification, not blind empirical prediction', other_rules='No frozen positive claim; full clinical restoration UNKNOWN'), files=files, code={str(p.relative_to(H)): sha(p) for p in (H / 'code').rglob('*.py')}, physical_measurement_performed=False)
    dump(H / 'FROZEN_PREDICTIONS.json', predictions)
    (H / 'FROZEN_PREDICTIONS.sha256').write_text(sha(H / 'FROZEN_PREDICTIONS.json') + '\n')
    dump(H / 'CURRENT_WORK_STATE.json', dict(lane='X34-stl-design-gate', status='R1_PREDICTIONS_FROZEN', latest_gate='Closed-form fixtures and retrospective real gap frozen; no validation run yet', next_operation='Run injected geometry, metadata, coordinate transforms and inherited full controls', updated_utc=predictions['frozen_utc']))
    print('Frozen inputs:', len(files), 'local data bytes:', sum((x['bytes'] for x in files)))
if __name__ == '__main__':
    main()
