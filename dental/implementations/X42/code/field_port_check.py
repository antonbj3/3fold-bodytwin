from common import *
from model_fit import loadmodel
from field_adapter import generate_mesh, cpu_field
import resource
FIELD = ROOT.parent.parent / 'references/field_engine/src/field_engine'
sys.path.insert(0, str(FIELD))

def run():
    tic = time.perf_counter()
    key = split()['dev_round1'][0]
    cfg = json.load(open(ROOT / 'FINAL_CONFIG.json'))
    b = cfg['branch_methods']['x42_contact_branch']
    base = DATA / 'FINAL_MODELS'
    rows = []
    ports = {}
    for fam in ['molar_crown', 'lattice_onlay']:
        orig = next((t for t in tasks(key) if t['family'] == fam and t['level'] == 'normal'))
        t = scene(orig)
        ports[fam] = generate_mesh(t, loadmodel(base / b['shape'] / (fam + '.npz')), loadmodel(base / b['contact'] / (fam + '.npz')), b['beta'])
    p = ports['molar_crown']
    v = p['outer_vertices']
    f = next((scene(task) for task in tasks(key) if task['family'] == 'molar_crown' and task['level'] == 'normal'))['faces']
    xy = v[:, :2]
    area = np.abs(np.cross(xy[f[:, 1]] - xy[f[:, 0]], xy[f[:, 2]] - xy[f[:, 0]])) / 2
    height = v[:, 2] - p['inner_vertices'][:, 2]
    exact = float(np.sum(area * np.mean(height[f], axis=1)))
    lo = p['vertices'].min(0) - 0.4
    frozen = {'frozen_utc': now(), 'analytic_volume_mm3': exact, 'candidate_generator_sha256': sha(ROOT / 'FROZEN_GENERATOR.json'), 'port_vertices_sha256': hashlib.sha256(np.ascontiguousarray(p['vertices']).tobytes()).hexdigest(), 'port_faces_sha256': hashlib.sha256(np.ascontiguousarray(p['faces']).tobytes()).hexdigest(), 'engineering_fixture_not_physical_measurement': True, 'predicted_winding_parity_disagreement_voxels': 0, 'tolerance_relative_volume_at0p1mm': 0.05}
    dump(ROOT / 'FROZEN_CPU_FIELD_PREDICTIONS.json', frozen)
    (ROOT / 'FROZEN_CPU_FIELD_PREDICTIONS.sha256').write_text(sha(ROOT / 'FROZEN_CPU_FIELD_PREDICTIONS.json') + '  FROZEN_CPU_FIELD_PREDICTIONS.json\n')
    for pitch in [0.2, 0.1]:
        ((gmin, shape, surface, solid, sd), diag) = cpu_field(p, pitch, lo)
        ((g2, s2, sf2, par, sd2), diag2) = cpu_field(p, pitch, lo, 'raypar_paritet')
        volume = float(solid.sum() * pitch ** 3)
        delta = int(np.sum(solid != par))
        fn = DATA / ('CPU_FIELD_' + str(pitch) + '.npz')
        np.savez_compressed(fn, gmin=gmin, origin=lo + np.array(gmin) * pitch, pitch=pitch, solid=solid, sd=sd)
        rows.append(dict(pitch_mm=pitch, shape=list(shape), analytic_volume_mm3=exact, voxel_volume_mm3=volume, relative_volume_error=abs(volume - exact) / exact, winding_parity_disagreement_voxels=delta, sha256=sha(fn), path=str(fn), bytes=fn.stat().st_size, diagnostic=diag, resolution='PER_TOOTH'))
    guards = []
    for (name, port) in [('actual_open_lattice', ports['lattice_onlay']), ('injected_missing_facet', dict(p, faces=p['faces'][:-1]))]:
        try:
            cpu_field(port, 0.1, lo)
            reject = False
        except ValueError:
            reject = True
        guards.append(dict(name=name, rejected=reject))
    good = rows[-1]['relative_volume_error'] <= 0.05 and all((r['winding_parity_disagreement_voxels'] == 0 for r in rows)) and all((r['rejected'] for r in guards))
    out = dict(PASS=good, engineering_fixture=True, external_referent=json.load(open(ROOT / 'PREREG_FIELD_PORT.json'))['external_referent'], rows=rows, guards=guards, seconds=time.perf_counter() - tic, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, field_source_sha256=sha(FIELD / 'faltkarna_v1_mesh_to_sdf.py'), native_GPU_execution='NOT_RUN', physical_cement_fit='UNKNOWN; pitches coarser than40um nominal film', prereg_sha256=sha(ROOT / 'PREREG_FIELD_PORT.json'))
    dump(ROOT / 'raw/CPU_FIELD_PORT.json', out)
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    run()
