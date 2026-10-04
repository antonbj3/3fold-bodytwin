from pose_join import *
from measure_geometry import closest
import trimesh

def sample(t, n=8192):
    area = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2
    idx = np.searchsorted(np.cumsum(area), (np.arange(n) + 0.5) * area.sum() / n)
    return t[idx].mean(1)

def run():
    pr = json.loads((ROOT / 'PREREG_VALIDATION.json').read_text())
    metrics = pr['metrics']
    exports = json.loads((ROOT / 'inputs/R4_FROZEN_EXPORTS.json').read_text())['exports']
    inputs = {r['key']: r for r in json.loads((ROOT / 'inputs/R4_FROZEN_INPUTS_B.json').read_text())['records']}
    r4 = json.loads((ROOT / 'raw/R4_RESULTS.json').read_text())
    rows = []
    for (e, row) in zip(exports, r4['rows']):
        source = inputs[e['key']]
        assert sha(source['private_path']) == source['private_sha256']
        target = np.load(source['private_path'], allow_pickle=False)['target']
        r = dict(np.load(e['mesh_path'], allow_pickle=False))
        outer = r['vertices'][r['faces'][r['roles'] == 0]]
        st = time.perf_counter()
        d1 = closest(trimesh.Trimesh(target.reshape(-1, 3), np.arange(len(target) * 3).reshape(-1, 3), process=False), sample(outer))[1]
        d2 = closest(trimesh.Trimesh(outer.reshape(-1, 3), np.arange(len(outer) * 3).reshape(-1, 3), process=False), sample(target))[1]
        p95 = float(max(np.quantile(d1, 0.95), np.quantile(d2, 0.95)))
        expected = e['predictions']['in_situ']['p95_mm']
        error = abs(p95 - expected)
        assert error <= metrics['original_surface_p95_vs_R4_frozen_mm_max']
        a = dict(np.load(ROOT / row['array_path'], allow_pickle=False))
        m = trimesh.Trimesh(a['vertices'], a['faces'], process=False)
        loaded = trimesh.load(ROOT / row['stl_path'], force='mesh', process=True)
        rt = closest(loaded, m.vertices)[1]
        rterr = float(rt.max())
        volerr = abs(float(loaded.volume - m.volume))
        assert rterr <= metrics['stl_roundtrip_max_distance_mm'] and volerr <= metrics['stl_volume_difference_mm3_max'] and loaded.is_watertight
        broken = trimesh.Trimesh(m.vertices, m.faces[:-1], process=False)
        assert not broken.is_watertight
        rows.append({'family': e['family'], 'original_surface_p95_mm': p95, 'inherited_reference_p95_mm': expected, 'difference_mm': error, 'original_reference_path': source['private_path'], 'original_reference_sha256': source['private_sha256'], 'resolution': 'PER_TOOTH summary of PER_POINT', 'stl_roundtrip_max_mm': rterr, 'stl_volume_difference_mm3': volerr, 'stl_loaded_watertight': bool(loaded.is_watertight), 'injected_open_face_rejected': True, 'source_comparison_seconds': time.perf_counter() - st, 'scope': 'natural anatomical source frame, no clinical tolerance or real preparation'})
        print('reference', e['family'], p95, 'exporterror', rterr, flush=True)
    out = {'status': 'PASS_SCOPED_REFERENCE_AND_EXPORT_CONTROL', 'rows': rows, 'physical_measurement': 'NOT_PERFORMED'}
    dump(ROOT / 'raw/FINAL_VALIDATION.json', out)
    return out
if __name__ == '__main__':
    run()
