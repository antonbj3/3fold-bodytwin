from fc_common import *
import trimesh, resource
sys.path.insert(0, str(V4 / 'code'))
sys.path.insert(0, str(V4 / 'payload/participant_code'))
sys.path.insert(0, str(V4 / 'payload/participant_code/legacy/vendor'))
from source_ops import pair

def run():
    st = time.perf_counter()
    rows = []
    pub = DATA / 'public'
    pub.mkdir(exist_ok=True)
    meta = read(V4 / 'payload/whole_inputs/RECORDS.json')
    pr = read(ROOT / 'PREREG_R2.json')
    for sel in read(V4 / 'PREREG_R3.json')['selection']:
        (arches, sources) = pair(sel)
        low = arches['lower']
        for rec in meta:
            if rec['case_key'] != sel['case_key']:
                continue
            key = rec['key']
            p = npz(V4 / 'payload/whole_inputs' / key / 'preparation.npz')
            a = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
            tri = a['source_triangles']
            m = float(p['margin_z'])
            mask = (tri[:, :, 2].max(1) <= m + 0.5) & (tri[:, :, 2].min(1) >= m - 1.0)
            collar = tri[mask]
            other = rec['source_fdi'] + 10
            local = (low['v'] - p['source_base']) @ p['source_R']
            donor = local[low['f'][low['owner'] == other]]
            out = dict(p, collar_triangles=collar, homolog_triangles=donor, mesial_triangles=a['mesial_triangles'], distal_triangles=a['distal_triangles'], antagonist_triangles=a['antagonist_triangles'])
            dest = pub / (key + '.npz')
            np.savez_compressed(dest, **out)
            rows.append(dict(rec, homolog_fdi=other, collar_faces=len(collar), homolog_faces=len(donor), public_file_sha256=sha(dest), sources=sources, collar_discarded_fraction=1 - float(mask.mean()), collar_upper_mm=m + 0.5, collar_lower_mm=m - 1.0))
            print(key, 'public collar', len(collar), 'homolog', len(donor), flush=True)
        del arches, low
    dump(pub / 'RECORDS.json', rows)
    freeze(ROOT / 'FROZEN_PUBLIC_INPUTS.json', dict(rows=rows, files={p.name: {'sha256': sha(p), 'bytes': p.stat().st_size} for p in pub.glob('*')}, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
if __name__ == '__main__':
    run()
