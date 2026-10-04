"""Query the candidate acquired pair. Physical/independent registration remains ABSTAIN."""
from dental_release.paths import expand as _release_expand
import argparse
from registration import *

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--jaw', choices=['upper', 'lower'], default='lower')
    ap.add_argument('--points', help='JSON N x3 IOS native points, assumed mm')
    ap.add_argument('--demo', action='store_true')
    args = ap.parse_args()
    w = np.load(DATA / 'r2_witnesses.npz')
    jaw = args.jaw
    T = w[jaw + '_T']
    p = w[jaw + '_source_points']
    fdi = w[jaw + '_fdi']
    ids = []
    if args.demo:
        ids = [int(np.flatnonzero(fdi == i)[0]) for i in sorted(set(fdi))]
        query = p[ids]
    elif args.points:
        query = np.array(json.loads(pathlib.Path(args.points).read_text()), float)
    else:
        ap.error('--points or --demo required')
    if query.ndim != 2 or query.shape[1] != 3 or (not np.isfinite(query).all()):
        raise ValueError('FINITE_XYZ_REQUIRED')
    declared = json.loads((ROOT / 'raw/R2_RESULTS.json').read_text())
    if sha(DATA / 'r2_witnesses.npz') != declared['witness_sha256']:
        raise ValueError('STALE_REGISTRATION_ARRAY')
    source = json.loads((ROOT / 'raw/HAO_ACQUISITION_MANIFEST.json').read_text())['files']
    cbct = next((x for x in source if x['member'].endswith('out_smoothed.stl')))
    if sha(cbct['local_path']) != cbct['sha256']:
        raise ValueError('STALE_CBCT_SOURCE')
    tri = stl(DATA / 'hao_demo/Demo/Demo_1/CBCT/reconstruction/out_smoothed.stl')
    cent = tri.mean(1)
    tri = tri[cent[:, 2] < 39] if jaw == 'lower' else tri[cent[:, 2] >= 39]
    mapped = transform(query, T)
    (d, closest, cell) = nearest_surface(mapped, tri.reshape(-1, 3), np.arange(tri.size // 3).reshape(-1, 3))
    result = {'patient_id': _release_expand('@DENTAL_CASE_ID@'), 'jaw': jaw, 'unit': 'mm_ASSUMED_SOURCE_STL_SCALE', 'resolution': 'PER_POINT', 'time_scale': 'SIMULTANEOUS', 'physical_acquisition_timing': 'UNKNOWN', 'source_frame': 'IOS_NATIVE_' + jaw, 'target_frame': 'CBCT_DEMO1', 'source_to_target': T.tolist(), 'mapped_CBCT_xyz_mm': mapped.tolist(), 'closest_CBCT_surface_xyz_mm': closest.tolist(), 'CBCT_triangle_source_address': cell.tolist(), 'surface_discrepancy_mm': d.tolist(), 'source_FDI': fdi[ids].tolist() if ids else 'UNKNOWN_NO_SOURCE_ADDRESS', 'CBCT_root_FDI': 'UNKNOWN_NOT_ANNOTATED', 'independent_registration': 'ABSTAIN_METRIC_CALIBRATION_AND_LANDMARK_TRUTH_MISSING', 'uniform_registration_error_bound_mm': None, 'physical_prediction': 'UNKNOWN'}
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
