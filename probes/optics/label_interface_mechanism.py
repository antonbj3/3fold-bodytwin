"""Count exact synthetic label interfaces before a shared-face mesh design."""
import hashlib
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]


def fixtures():
    layer=np.ones((2,2,2),np.uint8);layer[1:]=2
    cavity=np.ones((3,3,3),np.uint8);cavity[1,1,1]=2
    edge=np.zeros((2,2,1),np.uint8);edge[0,0,0]=edge[1,1,0]=1
    return dict(block=np.ones((2,2,2),np.uint8),layers=layer,cavity=cavity,edge_contact=edge)


def measure():
    rows=[]
    for name,labels in fixtures().items():
        volume={str(int(k)):int(np.count_nonzero(labels==k)) for k in np.unique(labels) if k}
        pairs={};boundaries={k:0 for k in volume};padded=np.pad(labels,1)
        for axis in range(3):
            left=[slice(1,-1)]*3;right=left.copy();left[axis]=slice(0,-1);right[axis]=slice(1,None)
            a=padded[tuple(left)];b=padded[tuple(right)];mask=a!=b
            for x,y in zip(a[mask],b[mask]):
                key=f'{min(x,y)}:{max(x,y)}';pairs[key]=pairs.get(key,0)+1
                for label in (x,y):
                    if label:boundaries[str(int(label))]+=1
        rows.append(dict(case=name,sha256=hashlib.sha256(labels.tobytes()).hexdigest(),shape=list(labels.shape),
                         region_voxels=volume,pair_faces=pairs,region_boundary_faces=boundaries))
    return rows


if __name__=='__main__':
    a,b=measure(),measure();gates=dict(full_repeat=a==b,nonzero_boundaries=all(n>0 for r in a for n in r['region_boundary_faces'].values()))
    report=dict(rows=a,gates=gates,scope='Synthetic face incidence only, no manifold certificate.')
    (ROOT/'reports/label_interface_mechanism.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));raise SystemExit(0 if all(gates.values()) else 1)
