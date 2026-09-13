"""Measure a mesh before designing or applying the ingest adapter."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import trimesh


def measure(path):
    raw=Path(path).read_bytes()
    m=trimesh.load_mesh(path,process=False)
    if not isinstance(m,trimesh.Trimesh):raise ValueError('One triangle mesh required')
    edges=np.sort(m.edges,axis=1);_,counts=np.unique(edges,axis=0,return_counts=True)
    return {'source_sha256':hashlib.sha256(raw).hexdigest(),'vertices':len(m.vertices),
            'faces':len(m.faces),'watertight':bool(m.is_watertight),
            'winding_consistent':bool(m.is_winding_consistent),
            'unpaired_edges':int(np.count_nonzero(counts!=2)),
            'degenerate_faces':int(np.count_nonzero(m.area_faces<=0)),
            'bounds_mm':m.bounds.tolist(),'volume_mm3':float(m.volume),
            'surface_area_mm2':float(m.area),'finite':bool(np.isfinite(m.vertices).all()),
            'array_sha256':hashlib.sha256(m.vertices.tobytes()+m.faces.tobytes()).hexdigest()}


def main():
    p=argparse.ArgumentParser();p.add_argument('mesh');p.add_argument('--units',required=True,choices=['mm'])
    p.add_argument('--out',required=True);args=p.parse_args()
    rows=[measure(args.mesh) for _ in range(2)]
    gates={'exact_repeats':rows[0]==rows[1],
           'closed_oriented':all(r['watertight'] and r['winding_consistent'] and r['unpaired_edges']==0 for r in rows),
           'valid_geometry':all(r['finite'] and r['degenerate_faces']==0 and r['volume_mm3']>0 for r in rows)}
    report={'units':args.units,'legs':rows,'gates':gates}
    Path(args.out).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
