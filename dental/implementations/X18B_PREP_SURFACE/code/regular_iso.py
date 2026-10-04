"""Conservative regular-level preparation export after preserved R4 singularity."""
from common import *
from implicit_prep import mesh_field
from scipy import ndimage
from scipy.spatial import cKDTree
import trimesh

def run():
    old = json.loads((H / 'rounds/R4.json').read_text())
    assert not old['gates']['manifold']
    path = H / 'PREREG_R5.json'
    if not path.exists():
        dump(path, dict(round='R5', claim_type='capability', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), capability='Export a closed preparation specimen at a strictly interior generic distance level after R4 exact-lattice tangency', obstacle='R4 depth1.2mm intersects a lattice critical level, leaving one edge with4 incident faces', changed_operation='Choose1.205mm conservative erosion level before extraction; preparation shrinks relative to1.2mm, retaining the same measured exterior and explicit virtual base. Do not repair/drop triangles.', erosion_depth_mm=1.205, nominal_depth_mm=1.2, metrics=dict(manifold_required=True, distance_tolerance_mm=0.2, exterior_max_mm=1e-06), strongest_equally_informed_control='R4 original level retained, independent edge incidence count and nearest-background distances', falsifiers=['Any edge with incidence other than2 admitted', 'Injected single face removal passes manifold gate', 'Generic-level perturbation called new measured anatomy'], external_referent=dict(kind='published_code', locator='https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.marching_cubes', compared_quantity='Isosurface mesh of existing frozen scalar field; geometry validity checked by independent edge incidence', refutes_us=True), physical_preparation_facit='UNKNOWN; source scan only supports exterior, not virtual cervical plane', resolution='PER_TOOTH', full_cost=dict(preparation='Reuse frozen R4 field with hash check', fit='none', discovery='R4 failure preserved', validation='Independent topology/distance and removed-face injection', questions='One iso extraction', fallback='FAIL for any topological defect'), limits=dict(threads=1, ram_gb=3, intermediate_bytes=3000000000)))
    pr = json.loads(path.read_text())
    verify_freeze('R4')
    a = np.load(old['prediction']['field_arrays']['path'])
    step = float(a['grid_step_mm'])
    depth = pr['erosion_depth_mm']
    (v, f) = mesh_field(a['phi'], a['origin'], step, -depth)
    p = D / '127_36_regular_implicit_preparation.stl'
    write_stl(p, v, f)
    raw = D / '127_36_regular_implicit_mesh.npz'
    np.savez_compressed(raw, vertices=v, faces=f)
    files = [artifact(p), artifact(raw), old['prediction']['field_arrays']]
    dump(H / 'FROZEN_PREDICTIONS_R5.json', dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(path), files=files, code_sha256={str(H / 'code/regular_iso.py'): sha(H / 'code/regular_iso.py')}, physical_measurement_status='NOT_RUN'))

    def nbad(ff):
        e = np.sort(np.concatenate([ff[:, [0, 1]], ff[:, [1, 2]], ff[:, [2, 0]]]), axis=1)
        ct = np.unique(e, axis=0, return_counts=True)[1]
        return int(np.sum(ct != 2))
    nm = nbad(f)
    bad = nbad(f[:-1])
    occ = a['occupied']
    background = np.argwhere(ndimage.binary_dilation(occ, structure=np.ones((3, 3, 3))) & ~occ) * step + a['origin']
    distance = cKDTree(background).query(v, workers=1)[0]
    err = float(np.max(np.abs(distance - depth)))
    m = trimesh.Trimesh(vertices=v, faces=f, process=False)
    gates = dict(manifold=nm == 0 and bool(m.is_watertight), distance=err <= pr['metrics']['distance_tolerance_mm'], removed_face_injection=bad > 0)
    out = dict(round='R5', claim_type='capability', external_referent=pr['external_referent'], decision='REGULAR_SPECIMEN_PASS' if all(gates.values()) else 'REGULAR_SPECIMEN_FAIL', gates=gates, depth_mm=depth, preparation_stl=artifact(p), nonmanifold_edges=nm, injected_missing_face_bad_edges=bad, nearest_background_depth_error_mm=err, vertices=len(v), faces=len(f), physical_preparation_validation='UNKNOWN; virtual basal plane and insertion envelope remain closures', scope='Numerical export correction after failed R4; no physical or algorithm superiority claim')
    dump(H / 'rounds/R5.json', out)
    (H / 'HANDOFF_R5.md').write_text(f"R5 {out['decision']}: conservative generic1.205mm erosion level, nonmanifold edges={nm}, injected face removal gives{bad} bad edges. R4 critical-level failure preserved. Next construction requires independent cervical/tooth annotations and actual prep+registered pressure measurement; no fabricated acquisition.\n")
    state('R5_DECIDED', out['decision'], 'Package source carrier, force intervals, calibration ports and manufactured-specimen limits')
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    run()
