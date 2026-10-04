from common_r3 import *
sys.path.insert(0, str(V6 / 'code'))
from contact import compare, gates
rows = []
for r in read(DATA / 'inputs/RECORDS.json'):
    p = npz(DATA / 'inputs' / (r['key'] + '.npz'))
    g = p['contact_reference_gap_unclipped']
    c = compare(p['contact_xy'], p['contact_faces'], g, g)
    rows.append(dict(key=r['key'], family=r['family'], native_self=c, native_self_gates=gates(c, read(ROOT / 'PREREG_A.json')['metrics']['v6'])))
dump(ROOT / 'raw/NATIVE_SELF_CONTACT.json', dict(rows=rows, scope='read existing frozen acquisition; native against itself tests whether original digital conjunction is satisfiable by copying native'))
for r in rows:
    print(r['key'], r['native_self_gates'], r['native_self'].get('negative_gap_area_mm2'))
