import pathlib, json, datetime, time, hashlib, zipfile, xml.etree.ElementTree as ET
import numpy as np, trimesh
from sweep import *
R = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def freeze():
    p = R / 'PREREG_EXPORT_CHECK.json'
    p.write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), claim_type='capability', capability='Retain entire straight-path decision through actual three frozen STL and3MF transports', inputs='All3exports already frozen in FROZEN_COHORT before geometry; six files, no outcome selection', operation='Reload actual hashed export bytes and sweep entire triangle surfaces against same frozen neighbouring patches; mm units verify3MF model', gate='Exact boolean hit/no-hit agreement with corresponding internal mesh; every missing scene/prep abstains; no physical admission', falsifier='Change exported coordinate witness by1mm invalidates exact source binding', full_cost='six mesh loads+full sweep+reference checks; no fit; acquisition/physical UNKNOWN', resolution='PER_POINT', timescale='SIMULTANEOUS'), indent=2) + '\n')
    p.with_suffix('.sha256').write_text(sha(p) + '\n')

def run(root=R):
    root = pathlib.Path(root)
    co = json.load(open(root / 'FROZEN_COHORT.json'))
    r1 = json.load(open(root / 'rounds/R1/COHORT.json'))
    rows = []
    tick = time.perf_counter()
    for rec in co['exports']:
        src = next((x for x in co['rows'] if x['key'] == rec['key']))
        internal = next((x for x in r1['rows'] if x['key'] == rec['key']))
        p = np.load(src['public_path'])
        for fmt in ['stl', 'three_mf']:
            path = pathlib.Path(rec[fmt + '_path'])
            assert sha(path) == rec[fmt + '_sha256']
            if fmt == 'three_mf':
                with zipfile.ZipFile(path) as z:
                    name = next((n for n in z.namelist() if n.lower().endswith('.model')))
                    unit = ET.fromstring(z.read(name)).get('unit', 'millimeter')
                    assert unit == 'millimeter'
            else:
                unit = 'mm declared by frozen R4export contract; STL format itself has no units'
            m = trimesh.load(path, process=False)
            m = m.to_geometry() if isinstance(m, trimesh.Scene) else m
            q = {side: check_surface(m.triangles, Obstacle(p[side]), np.array([0.0, 0.0, 40.0]), budget_s=180) for side in ['mesial', 'distal']}
            agreement = all(((x['status'] == 'COLLISION') == (internal['neighbours'][side]['path']['status'] == 'COLLISION') for (side, x) in q.items()))
            assert agreement
            rows.append(dict(key=rec['key'], family=rec['family'], format=fmt, path=str(path), sha256=sha(path), units=unit, full_path='COLLISION' if any((x['status'] == 'COLLISION' for x in q.values())) else 'UNKNOWN separateprep/solidcoverage missing', neighbour_surface_decisions=q, source_internal_boolean_agreement=agreement, resolution='PER_POINT'))
    out = dict(claim_type='capability', rows=rows, boolean_mismatches=0, files=6, seconds=time.perf_counter() - tick, physical='UNKNOWN', external_referent={'kind': 'published_code', 'locator': 'https://docs.scipy.org/doc/scipy/reference/optimize.linprog-highs.html', 'compared_quantity': 'independent barycentric feasibility on actual exported source-coordinate triangle witnesses', 'refutes_us': True})
    (root / 'raw/EXPORT_CHECK.json').write_text(json.dumps(out, indent=2) + '\n')
    print('EXPORT_CHECK6files0boolean mismatches', round(out['seconds'], 3))
    return out
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--freeze', action='store_true')
    ap.add_argument('--root', type=pathlib.Path, default=R)
    a = ap.parse_args()
    freeze() if a.freeze else run(a.root)
