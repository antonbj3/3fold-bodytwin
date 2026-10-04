"""Export candidate aligned crown STL plus provenance sidecar; source unit assumed."""
from registration import *

def main():
    r = json.loads((ROOT / 'raw/R2_RESULTS.json').read_text())
    dst = DATA / 'exports'
    dst.mkdir(exist_ok=True)
    out = []
    for j in r['jaws']:
        pieces = []
        owners = []
        T = np.array(j['source_to_target'])
        for s in j['source_labels']:
            if sha(s['source']) != s['sha256']:
                raise ValueError('STALE_IOS_SOURCE')
            a = stl(s['source'])
            pieces.append(transform(a.reshape(-1, 3), T).reshape(-1, 3, 3))
            owners.extend([s['fdi']] * len(a))
        tri = np.concatenate(pieces)
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        norm = np.linalg.norm(n, axis=1)
        n[norm > 0] /= norm[norm > 0, None]
        dt = np.dtype([('normal', '<f4', (3,)), ('vertices', '<f4', (3, 3)), ('attribute', '<u2')])
        arr = np.zeros(len(tri), dtype=dt)
        arr['normal'] = n
        arr['vertices'] = tri
        path = dst / f"Demo1_{j['jaw']}_IOS_crowns_in_CBCT.stl"
        path.write_bytes(b'X69 candidate source-paired IOS crowns; assumed mm; anatomy accuracy UNKNOWN'.ljust(80, b' ') + np.array([len(tri)], dtype='<u4').tobytes() + arr.tobytes())
        owner = dst / f"Demo1_{j['jaw']}_facet_FDI.npy"
        np.save(owner, np.array(owners, dtype=np.int16))
        out.append({'jaw': j['jaw'], 'path': str(path), 'sha256': sha(path), 'facet_owner_path': str(owner), 'facet_owner_sha256': sha(owner), 'source_to_target': T.tolist(), 'source_labels': j['source_labels'], 'coordinate_frame': 'CBCT_DEMO1', 'metric_units': 'assumed mm; STL does not encode units', 'resolution': 'PER_POINT', 'time_scale': 'SIMULTANEOUS', 'acquisition_time_relation': 'UNKNOWN', 'role': 'Candidate aligned crown geometry for research inspection; no clinical/manufacture release', 'license': 'CC BY4.0', 'attribution': 'Jiaxiang Liu et al., dataset10.5281/zenodo.8027553; Patterns10.1016/j.patter.2023.100825; rigid transformations by X69'})
    dump(ROOT / 'EXPORT_SIDECAR.json', {'exports': out, 'CBCT_source': 'See raw/HAO_ACQUISITION_MANIFEST.json', 'clinical_release': 'ABSTAIN', 'source_and_CBCT_anatomical_accuracy': 'UNKNOWN', 'registration_trueness': 'UNKNOWN'})
if __name__ == '__main__':
    main()
