"""Independent registered-anatomy import demo, never joined to Teeth3DS tasks."""
from dental_release.paths import expand as _release_expand
import hashlib, json, time, zipfile
import numpy as np
from scipy.spatial import cKDTree
from .io import ROOT, DATA, dump, freeze, sha
from .checks.inherited import load
SOURCES = [('Bits2Bites', _release_expand('@DENTAL_DATA_ROOT@/geometry/Bits2Bites/Bits2Bites_v01.zip'), 2), ('Bite2Text', _release_expand('@DENTAL_DATA_ROOT@/geometry/Bite2Text/Bite2Text.zip'), 1)]

def selection():
    rows = []
    for (dataset, path, n) in SOURCES:
        with zipfile.ZipFile(path) as z:
            upper = sorted((s for s in z.namelist() if s.endswith('upper.stl')), key=lambda x: hashlib.sha256(x.encode()).hexdigest())[:n]
            for u in upper:
                l = u.replace('upper.stl', 'lower.stl')
                rows.append(dict(dataset=dataset, zip_path=path, upper_member=u, lower_member=l, upper_bytes=z.getinfo(u).file_size, lower_bytes=z.getinfo(l).file_size, upper_zip_crc=z.getinfo(u).CRC, lower_zip_crc=z.getinfo(l).CRC))
    return rows

def run():
    rows = selection()
    freeze(ROOT / 'data/REGISTERED_BITE_FROZEN.json', dict(selection=rows, metric='Symmetric nearest-vertex distance on deterministic subsamples, not clearance/contact/force', transforms='identity only; original paired coordinates', sample_count=20000, negative_control='Translate lower arch by [1000,0,0] mm; median nearest distance must increase by >=900 mm', geometry_not_linked_to_Teeth3DS=True))
    module = load('gencad_x2_import', 'LANE_X2_OCCLUSION_B2B/code/occlusion_operator.py')
    out = []
    tic = time.perf_counter()
    for (i, row) in enumerate(rows):
        with zipfile.ZipFile(row['zip_path']) as z:
            (U, uh) = module.read_stl(z, row['upper_member'])
            (L, lh) = module.read_stl(z, row['lower_member'])
        up = np.unique(U.reshape(-1, 3), axis=0)
        lo = np.unique(L.reshape(-1, 3), axis=0)
        up = up[np.linspace(0, len(up) - 1, min(20000, len(up)), dtype=int)]
        lo = lo[np.linspace(0, len(lo) - 1, min(20000, len(lo)), dtype=int)]
        distances = np.r_[cKDTree(lo).query(up, workers=4)[0], cKDTree(up).query(lo, workers=4)[0]]
        changed = lo + np.array([1000.0, 0, 0])
        bad = np.r_[cKDTree(changed).query(up, workers=4)[0], cKDTree(up).query(changed, workers=4)[0]]
        caught = float(np.median(bad) - np.median(distances)) >= 900
        assert caught
        case = hashlib.sha256((row['dataset'] + row['upper_member']).encode()).hexdigest()[:12]
        sampled = DATA / ('registered_' + case + '.npz')
        np.savez_compressed(sampled, upper=up, lower=lo)
        out.append(dict(dataset=row['dataset'], case_hash=case, upper_sha256=uh, lower_sha256=lh, triangles_upper=len(U), triangles_lower=len(L), native_pose_transform='identity', median_nearest_vertex_mm=float(np.median(distances)), minimum_nearest_vertex_mm=float(distances.min()), corrupt_pose_rejected=caught, corrupt_pose_median_mm=float(np.median(bad)), bound='Nearest-vertex distances overestimate true nearest-surface distance; no signed penetration conclusion', exported_samples=dict(path=str(sampled), sha256=sha(sampled)), reference_locator='https://ditto.ing.unimore.it/' + ('bits2bites/' if row['dataset'] == 'Bits2Bites' else 'bite2text/')))
    report = dict(cases=out, wall_s=time.perf_counter() - tic, registered_pairs_read=len(out), no_patient_join=True, physical_contact_and_force='UNKNOWN', external_referent=dict(kind='published_dataset', locator='https://ditto.ing.unimore.it/bits2bites/; https://ditto.ing.unimore.it/bite2text/', compared_quantity='native registered paired IOS geometry and pose integrity, not force or optimal crown', refutes_us=True))
    dump(ROOT / 'raw/registered_bite_results.json', report)
    return report
