from common_r3 import *
import time
sys.path.insert(0, str(V6 / 'code'))
from contact import compare
start = time.perf_counter()
out = DATA / 'inputs_C'
out.mkdir(exist_ok=True)
rows = []
for rec in read(DATA / 'inputs/RECORDS.json'):
    p = npz(DATA / 'inputs' / (rec['key'] + '.npz'))
    raw = p['contact_reference_gap_unclipped']
    ff = p['contact_faces']
    g = raw[ff]
    active = np.isfinite(g).all(1) & (g.min(1) <= 0.1) & (g.max(1) >= 0)
    ids = np.unique(ff[active])
    old = p['contact_gap_mm'].copy()
    p['contact_gap_mm'][ids] = raw[ids]
    c = compare(p['contact_xy'], ff, p['contact_gap_mm'], raw)
    np.savez_compressed(out / (rec['key'] + '.npz'), **p)
    row = dict(rec, augmented_nodes=int(np.count_nonzero(p['contact_gap_mm'][ids] != old[ids])), crossing_nodes=len(ids), acquisition_contact_comparison=c)
    rows.append(row)
dump(out / 'RECORDS.json', rows)
freeze(ROOT / 'FROZEN_ACQUISITION_C.json', dict(files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()}, rows=rows, seconds=time.perf_counter() - start, source_sha256=sha(ROOT / 'FROZEN_ACQUISITION_A.json'), role='Explicitly added raw endpoints on contact-crossing faces; no hidden source triangles given'))
xy = np.array([[0, 0], [1, 0], [0, 1]], float)
ff = np.array([[0, 1, 2]])
a = np.array([0.0625, 0.5, 0.5])
b = np.array([0.0625, 5, 5])
ca = np.clip(a, -0.25, 0.35)
cb = np.clip(b, -0.25, 0.35)
assert np.array_equal(ca, cb)
aa = compare(xy, ff, a, a)
bb = compare(xy, ff, b, a)
w = dict(identity_error=0.0, stored_gap_mm=ca, source_a_mm=a, source_b_mm=b, source_area_a_mm2=aa['predicted']['area_mm2'], source_area_b_mm2=bb['predicted']['area_mm2'], downstream_symdiff_mm2=bb['symdiff_mm2'], minimal_extension_for_pair='one crossing fraction on either symmetric edge; general PL contact requires crossing fractions on all affected edges', resolution='PER_POINT endpoints -> PER_SURFACE_REGION band', external_referent=dict(kind='our_own_fixture', locator='code/acquire_c.py', compared_quantity='same stored sheet, different source band', refutes_us=True))
dump(ROOT / 'raw/CLIPPING_SUFFICIENCY.json', w)
print(json.dumps(clean(dict(rows=len(rows), max_acquisition_symdiff_mm2=max((r['acquisition_contact_comparison']['symdiff_mm2'] for r in rows)), witness=w)), indent=2))
