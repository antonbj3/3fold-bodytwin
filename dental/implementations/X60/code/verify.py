from common import *
import sys, importlib.util, zipfile, xml.etree.ElementTree as ET
sys.path.insert(0, str(V4 / 'payload/participant_code'))
from task_io import attach
from legacy.checks import all_checks
from legacy.generators import submit
import trimesh

def roundtrip(path):
    with zipfile.ZipFile(path) as z:
        r = ET.fromstring(z.read('3D/3dmodel.model'))
    if r.attrib['unit'] != 'millimeter':
        raise ValueError('wrong 3MF unit')
    ns = {'c': 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'}
    v = np.array([[float(a.attrib[k]) for k in ['x', 'y', 'z']] for a in r.findall('.//c:vertex', ns)])
    f = np.array([[int(a.attrib[k]) for k in ['v1', 'v2', 'v3']] for a in r.findall('.//c:triangle', ns)])
    return trimesh.Trimesh(v, f, process=False)

def run():
    rows = read(ROOT / 'raw/SCORED_ROWS.json')
    ctrl = []
    done = set()
    for row in rows:
        if row['participant'] != 'constraint_optimizer' or row['L1'] != 'PASS':
            continue
        fam = row['family']
        if fam in done:
            continue
        done.add(fam)
        key = row['case_key']
        ts = read(V4 / 'payload/public/tasks' / f'{key}.json')
        i = next((j for (j, t) in enumerate(ts) if t['task_id'] == row['task_id']))
        t = ts[i]
        with np.load(V4 / 'payload/public' / t['geometry_file'], allow_pickle=False) as a:
            t = attach(t, dict(a))
        with np.load(V4 / 'payload/predictions/constraint_optimizer' / f'{key}.npz') as a:
            z = a['outer_' + str(i)]
            inner = a['inner_' + str(i)]
        rr = t['requirements']
        faults = [('wall', submit(t, inner + rr['wall_mm'] / 2, inner), 'FAIL'), ('film_min', submit(t, z, t['preparation_z'] - 1), 'FAIL'), ('film_max', submit(t, z, t['preparation_z'] + 1), 'FAIL'), ('nesting', submit(t, inner - 1, inner), 'FAIL')]
        if t['A'].shape[0]:
            shift = max(1.0, float(np.max(np.asarray(t['obstacle_b']) - t['A'] @ z)) + 1.0)
            bad = submit(t, z + shift, inner)
            faults.append(('antagonist', bad, 'FAIL'))
        if fam == 'bridge3':
            faults.append(('connector_area', submit(t, inner + 0.001, inner), 'FAIL'))
        x = np.asarray(t['xy'])[:, 0]
        tilt = 0.01 * (x - x.mean()) / max(abs(x - x.mean()).max(), 1.0)
        faults.append(('milling', submit(t, np.maximum(z, inner + rr['wall_mm'] + 1), t['preparation_z'] + 0.08 + tilt), 'UNKNOWN'))
        good = submit(t, z, inner)
        for (keyfield, badval) in [('units', 'cm'), ('frame', 'wrong_frame')]:
            bad = dict(good)
            bad[keyfield] = badval
            faults.append((keyfield, bad, 'INVALID'))
        bad = dict(good)
        vv = np.asarray(good['inner_vertices'])
        vv[:, 0] += 0.01
        bad['inner_vertices'] = vv.tolist()
        faults.append(('insertion_footprint', bad, 'INVALID'))
        for (field, bad, expected) in faults:
            result = all_checks(t, bad)
            actual = result['checks'].get(field, {}).get('status') if expected != 'INVALID' else result['validity']
            ctrl.append(dict(task_id=t['task_id'], field=field, expected=expected, actual=actual, overall=result['validity'], rejected=actual == expected, source_checker_sha256=sha(V4 / 'payload/participant_code/legacy/checks.py'), external_referent_kind='our_own_fixture'))
    exportchecks = []
    for path in sorted((ROOT / 'exports').glob('whole*/*.stl')):
        m = trimesh.load_mesh(path, process=True)
        other = roundtrip(path.with_suffix('.3mf'))
        dv = abs(m.volume / other.volume - 1)
        faultdir = DATA / 'verification_faults'
        faultdir.mkdir(exist_ok=True)
        badpath = faultdir / (path.parent.name + '_' + path.stem + '_wrong_unit.3mf')
        with zipfile.ZipFile(path.with_suffix('.3mf')) as zin, zipfile.ZipFile(badpath, 'w') as zout:
            for name in zin.namelist():
                blob = zin.read(name)
                if name == '3D/3dmodel.model':
                    blob = blob.replace(b'unit="millimeter"', b'unit="centimeter"')
                zout.writestr(name, blob)
        try:
            roundtrip(badpath)
            unit_reject = False
        except ValueError:
            unit_reject = True
        broken = trimesh.Trimesh(m.vertices, m.faces[1:], process=False)
        scaled = other.copy()
        scaled.apply_scale(1.6)
        badvolume = abs(m.volume / scaled.volume - 1)
        exportchecks.append(dict(path=str(path.relative_to(ROOT)), watertight=bool(m.is_watertight and other.is_watertight), winding_consistent=bool(m.is_winding_consistent and other.is_winding_consistent), volume_relative_error=dv, pass_gate=bool(m.is_watertight and other.is_watertight and (m.volume > 0) and (other.volume > 0) and (dv <= 1e-05)), wrong_units_rejected=unit_reject, missing_face_rejected=not broken.is_watertight, wrong_volume_rejected=badvolume > 1e-05))
    out = dict(valid_coordinate_faults=ctrl, all_faults_rejected=all((r['rejected'] for r in ctrl)), exports=exportchecks, all_exports_pass=all((r['pass_gate'] for r in exportchecks)), all_export_faults_rejected=all((r['wrong_units_rejected'] and r['missing_face_rejected'] and r['wrong_volume_rejected'] for r in exportchecks)), original_R2_control_result_preserved=read(ROOT / 'raw/R2.json')['all_faults_rejected'], full_milling_wall_antagonist_certificates='UNKNOWN', review_state='PENDING_INDEPENDENT_REVIEW')
    dump(ROOT / 'raw/VERIFICATION.json', out)
    print('valid-coordinate controls', len(ctrl), out['all_faults_rejected'], 'exports', out['all_exports_pass'])
    return out
if __name__ == '__main__':
    run()
