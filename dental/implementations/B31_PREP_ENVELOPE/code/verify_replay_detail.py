"""Compare substantive bytes and values, including negative-family checkpoints."""
from common import *
import sys

def compare(audit):
    audit = Path(audit)
    checks = []
    for round_ in ['R2', 'R3', 'R4']:
        src = ROOT / 'raw' / (round_ + '_RESULTS.json')
        dst = audit / 'raw' / (round_ + '_RESULTS.json')
        if not dst.exists():
            continue
        aa = {r['key']: r for r in read(src)}
        bb = {r['key']: r for r in read(dst)}
        assert set(aa) == set(bb)
        for (key, a) in aa.items():
            b = bb[key]
            q = {'round': round_, 'key': key}
            if round_ == 'R2':
                if 'field_path' not in b:
                    p = audit / 'raw/r2' / (key + '_FIELD_CHECKPOINT.json')
                    if p.exists():
                        b = {**b, **read(p)}
                        q['metadata_lookup'] = 'preserved negative checkpoint'
                q.update(status_same=a['status'] == b['status'], field_bytes_same=a.get('field_sha256') == b.get('field_sha256'), candidate_bytes_same=a.get('candidate_sha256') == b.get('candidate_sha256'), STL_bytes_same=a.get('stl_sha256') == b.get('stl_sha256'), exact_volume_same=a.get('canonical_intaglio', {}).get('exact_volume_mm3') == b.get('canonical_intaglio', {}).get('exact_volume_mm3'), whole_wall_values_same=a.get('wall') == b.get('wall'))
            elif round_ == 'R3':
                q.update(status_same=a['whole_path_status'] == b['whole_path_status'], flag_bytes_same=a.get('flag_sha256') == b.get('flag_sha256'), collision_count_same=a.get('collision_pairs') == b.get('collision_pairs'), LP_count_same=a.get('LP_controls') == b.get('LP_controls'))
            else:
                q.update(status_same=a['status'] == b['status'], witness_bytes_same=a.get('exact_witness_sha256') == b.get('exact_witness_sha256'), proper_crossing_count_same=a.get('proper_crossing_pairs') == b.get('proper_crossing_pairs'))
            checks.append(q)
    allpass = all((v for r in checks for (k, v) in r.items() if isinstance(v, bool)))
    out = {'all_pass': allpass, 'rows': checks, 'audit_result_sha256': sha(audit / 'results.json'), 'scope': 'exact geometry/field/STL/flags/witness bytes, exact-volume strings and full wall-result equality. V1 missing negative-row metadata resolved only by its already preserved checkpoint; original strict failed comparator remains archived.'}
    name = 'DETAIL_' + audit.parent.name + '.json'
    dump(ROOT / 'raw' / name, out)
    if not allpass:
        raise RuntimeError('Detailed replay changed; preserve ' + name)
    print(name, len(checks), 'all_pass', allpass)
if __name__ == '__main__':
    compare(sys.argv[1])
