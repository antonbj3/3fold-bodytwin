"""R6: preserved author-default run + source-context ownership adapter.
Generation uses public inputs only. Scoring begins strictly after output freeze.
"""
from common import *
from scipy.spatial import cKDTree
from scipy import ndimage
from skimage.measure import marching_cubes
import trimesh
from contracts import observation_gate

def context_measure(vertices, faces, neighbor_triangles):
    from functional import sample, distance, make_mesh
    points = sample(neighbor_triangles, 512)
    dist = distance(make_mesh(vertices[faces]), points)
    error = float(np.quantile(dist, 0.95)) if len(dist) else None
    return dict(p95_mm=error, pass_gate=error is not None and observation_gate(0, error, 0.25), probes=len(points))

def run():
    from functional import score_mesh, distance, make_mesh, sample
    selected = read(ROOT / 'EXTERNAL_INPUTS.json')['rows']
    dest = DATA / 'external_support'
    dest.mkdir(exist_ok=True)
    generated = []
    for r in selected:
        key = r['key']
        p = npz(DATA / 'external_inputs' / f'{key}.npz')
        pred = npz(DATA / 'external_predictions_100' / f'{key}.npz')['tsdf']
        scale = float(p['scale'])
        axis = -1 + np.arange(64) / 32
        world = np.stack(np.meshgrid(axis, axis, axis, indexing='ij'), -1) / scale + p['center']
        flat = world.reshape(-1, 3)
        dt = cKDTree(p['preparation_vertices']).query(flat, workers=1)[0].reshape((64,) * 3)
        dn = cKDTree(p['neighbor_triangles'].reshape(-1, 3)).query(flat, workers=1)[0].reshape((64,) * 3)
        owner = dt <= dn
        preserved = np.where(owner, pred, p['tsdf'])
        identity = float(np.max(abs(preserved[~owner] - p['tsdf'][~owner])))
        for tag in ['author100', 'context_owned100']:
            folder = dest / tag / key
            folder.mkdir(parents=True, exist_ok=True)
            crop = np.maximum(np.max(p['roi_low'] - world, axis=-1), np.max(world - p['roi_high'], axis=-1))
            exterior = np.maximum(pred / scale, crop)
            if tag == 'context_owned100':
                exterior = np.maximum(exterior, dt - dn)
            field = np.maximum(exterior, -p['cavity'])
            (v, f, _, _) = marching_cubes(field, 0, spacing=(2 / 64 / scale,) * 3, allow_degenerate=False)
            v += -1 / scale + p['center']
            pts = v[f].mean(1)
            ix = ((pts - p['center']) * scale + 1) * 32
            vals = np.stack([ndimage.map_coordinates(exterior, ix.T, order=1), ndimage.map_coordinates(-p['cavity'], ix.T, order=1)])
            roles = np.argmax(vals, axis=0).astype(np.int8)
            if trimesh.Trimesh(v, f, process=False).volume < 0:
                f = f[:, ::-1]
            np.savez_compressed(folder / 'mesh.npz', vertices=v, faces=f, face_roles=roles)
            trimesh.Trimesh(v @ p['source_R'].T + p['source_base'], f, process=False).export(folder / 'crown.stl')
            generated.append(dict(key=key, tag=tag, path=str(folder / 'mesh.npz'), sha256=sha(folder / 'mesh.npz'), outside_context_identity_error=identity if tag == 'context_owned100' else None, normalization_scale=scale, owner_scope='Nearest sampled public surface vertices, discrete ownership only; not a rigorous anatomical margin'))
        np.savez_compressed(dest / (key + '_context.npz'), preserved_context=preserved, owner=owner)
    frozen = ROOT / 'FROZEN_PREDICTIONS_R6.json'
    if not frozen.exists():
        freeze(frozen, dict(rows=generated, input_lock_sha256=sha(ROOT / 'EXTERNAL_INPUTS.json'), prereg_sha256=sha(ROOT / 'PREREG_R6_EXTERNAL_SUPPORT.json'), old100_prediction_sha256=sha(ROOT / 'history/external_100_steps/FROZEN_PREDICTIONS.json'), source_target_queries_before_freeze=0))
    elif read(frozen)['rows'] != clean(generated):
        raise ValueError('R6 frozen geometry drift')
    scored = []
    context_tests = []
    md = {r['key']: r for r in selected}
    for r in generated:
        m = npz(r['path'])
        s = score_mesh(dict(md[r['key']], participant='ToothCraft_' + r['tag']), m['vertices'], m['faces'], m['face_roles'])
        key = r['key']
        p = npz(DATA / 'external_inputs' / f'{key}.npz')
        scale = float(p['scale'])
        whole = npz(DATA / 'external_predictions_100' / f'{key}.npz')['tsdf'] if r['tag'] == 'author100' else npz(dest / (key + '_context.npz'))['preserved_context']
        (cv, cf, _, _) = marching_cubes(whole, 0, spacing=(2 / 64 / scale,) * 3, allow_degenerate=False)
        cv += -1 / scale + p['center']
        cm = context_measure(cv, cf, p['neighbor_triangles'])
        error = cm['p95_mm']
        if r['tag'] == 'context_owned100':
            bad = context_measure(cv, cf, p['neighbor_triangles'] + np.array([20.0, 20.0, 20.0]))
            context_tests.append(dict(key=key, normal_accepted=cm['pass_gate'], injected_translated_observation_rejected=not bad['pass_gate'], normal_p95_mm=error, bad_p95_mm=bad['p95_mm']))
        scored.append(dict(track='R6', **s, known_context_identity_error=r['outside_context_identity_error'], known_neighbour_surface_p95_mm=error, known_neighbour_surface_gate=cm['pass_gate'], known_neighbour_probes=cm['probes'], known_neighbour_enclosure='MISSING continuous sampling bound'))
    primary = {r['key']: r for r in read(ROOT / 'raw/FUNCTIONAL_ROWS.json') if r['track'] == 'ToothCraft'}
    contrasts = []
    for key in md:
        a = next((r for r in scored if r['key'] == key and r['participant'] == 'ToothCraft_author100'))
        b = next((r for r in scored if r['key'] == key and r['participant'] == 'ToothCraft_context_owned100'))
        ten = primary[key]
        contrasts.append(dict(key=key, ten_step_p95_mm=ten['reconstruction_p95_mm'], hundred_step_p95_mm=a['reconstruction_p95_mm'], owned_hundred_p95_mm=b['reconstruction_p95_mm'], budget_improvement_mm=ten['reconstruction_p95_mm'] - a['reconstruction_p95_mm'], ownership_improvement_mm=a['reconstruction_p95_mm'] - b['reconstruction_p95_mm'], outside_context_identity_error=b['known_context_identity_error'], resolution='PER_TOOTH'))
    out = dict(claim_type='capability', rows=scored, contrasts=contrasts, budget_explains_failure=all((r['budget_improvement_mm'] >= 0.25 for r in contrasts)), known_context_discrete_preserved=all((r['outside_context_identity_error'] == 0 for r in contrasts)), geometric_context_surface_gate=all((r['known_neighbour_surface_gate'] for r in scored if r['participant'] == 'ToothCraft_context_owned100')), context_fault_tests=context_tests, physical_force_and_film='UNKNOWN', next_operation='Independently annotate target margin/coordinate alignment and verify source SDF preprocessing before any external-method ranking')
    dump(ROOT / 'raw/R6_EXTERNAL_SUPPORT.json', out)
    return out
if __name__ == '__main__':
    run()
