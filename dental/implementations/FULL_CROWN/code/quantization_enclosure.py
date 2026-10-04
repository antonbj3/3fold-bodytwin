"""Supplementary enclosure: ideal digital film versus the actual selected indexed exports."""
from dental_release.paths import expand as _release_expand
from fc_common import *
import trimesh
sys.path.insert(0, _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V5/deps'))
from prescan_solid import cavity_triangles, norm_interval
from export_lab import parse

def transform_scale_bound(R):
    lo = np.zeros((3, 3))
    hi = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            pp = R[:, i] * R[:, j]
            pl = np.nextafter(pp, -np.inf)
            ph = np.nextafter(pp, np.inf)
            l = np.nextafter(np.nextafter(pl[0] + pl[1], -np.inf) + pl[2], -np.inf)
            h = np.nextafter(np.nextafter(ph[0] + ph[1], np.inf) + ph[2], np.inf)
            lo[i, j] = np.nextafter(l - (i == j), -np.inf)
            hi[i, j] = np.nextafter(h - (i == j), np.inf)
    b = np.maximum(abs(lo), abs(hi))
    s = np.nextafter(np.nextafter(b[:, 0] + b[:, 1], np.inf) + b[:, 2], np.inf)
    delta = float(s.max())
    return (float(np.nextafter(np.sqrt(1 - delta), -np.inf)), float(np.nextafter(np.sqrt(1 + delta), np.inf)), delta)

def run():
    rows = []
    for case in read(ROOT / 'FROZEN_LAB_EXPORTS.json')['selected']:
        if 'path' not in case:
            continue
        path = Path(case['path'])
        c = case['certificates']
        p = npz(DATA / 'public' / (case['key'] + '.npz'))
        die = cavity_triangles(np.array(c['axis_top']), c['die_radius_mm'], float(p['margin_z']) - 1.0, closed=True)
        expected = die.vertices @ p['source_R'].T + p['source_base']
        blob = (path / 'die.stl').read_bytes()
        dt = np.dtype([('normal', '<f4', 3), ('v', '<f4', (3, 3)), ('attr', '<u2')])
        raw = np.frombuffer(blob, dtype=dt, offset=84)['v'].astype(float)
        assert raw.shape == expected[die.faces].shape
        eps_file = float(norm_interval((raw - expected[die.faces]).reshape(-1, 3))[1].max())
        (low, high, delta) = transform_scale_bound(p['source_R'])
        (cv, cf) = parse(path / 'crown.3mf')
        (dv, df) = parse(path / 'die.3mf')
        u = 2.0 ** (-53)
        gamma = 64 * u / (1 - 64 * u)
        B = max(np.max(abs(cv)), np.max(abs(dv)), 100.0) * 10
        rounding = float(np.nextafter(2 * np.sqrt(3) * gamma * B, np.inf))
        budget = float(np.nextafter(eps_file + rounding, np.inf))
        nom = c['normal_film_interval_mm']
        interval = [float(np.nextafter(low * nom[0] - budget, -np.inf)), float(np.nextafter(high * nom[1] + budget, np.inf))]
        rows.append(dict(key=case['key'], family=case['family'], resolution='PER_SURFACE_REGION', ideal_common_facet_film_mm=nom, actual_export_set_distance_film_enclosure_mm=interval, die_STL_quantization_upper_mm=eps_file, frame_RtR_departure_upper=delta, arithmetic_and_frame_rounding_allowance_mm=rounding, scope='common nonbasal paired surfaces in native modeled cavity domain; exact binary64 primitive/input model, IEEE754 arithmetic and frozen face correspondences; no physical cement', clinical_status='UNKNOWN', upper_expands_ideal=interval[1] > 0.08))
    result = dict(rows=rows, claim_type='capability', scope='This supplement qualifies the nominal film interval in earlier frozen certificates; actual die3MF inherits float32STL quantization. No earlier file or failed gate overwritten.', external_referent=dict(kind='closed_form', locator='code/quantization_enclosure.py; DECOMPOSITION_R4.json; standard1-Lipschitz distance and Hausdorff perturbation inequalities', compared_quantity='explicit outer bound for geometric gap after coordinate quantization and affine frame conversion', refutes_us=True))
    dump(ROOT / 'LAB_EXPORT_ENCLOSURES.json', result)
    return result
if __name__ == '__main__':
    run()
