"""Measure face multiplicity and sidedness before a tetrahedral adapter."""
import hashlib
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]


def fixtures():
    vertices=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1],[.25,.25,2]],float)
    pair=np.array([[0,1,2,3],[0,2,1,4]],np.int32)
    return dict(single=(vertices,pair[:1],np.array([1],np.int32)),
                same=(vertices,pair,np.array([1,1],np.int32)),
                layers=(vertices,pair,np.array([1,2],np.int32)),
                inverted=(vertices,pair[:,::-1][:,[1,0,2,3]],np.array([1,2],np.int32)),
                duplicate=(vertices,np.repeat(pair[:1],2,axis=0),np.array([1,2],np.int32)),
                triple=(vertices,np.vstack((pair,[0,1,2,5])),np.array([1,2,3],np.int32)))


def measure():
    rows=[]
    for name,(v,t,labels) in fixtures().items():
        p=v[t];six=np.einsum('ij,ij->i',p[:,1]-p[:,0],np.cross(p[:,2]-p[:,0],p[:,3]-p[:,0]))
        faces={}
        for i,tet in enumerate(t):
            for opposite in range(4):
                face=tuple(sorted(np.delete(tet,opposite)));faces.setdefault(face,[]).append((i,int(tet[opposite])))
        same_side=0;internal=0
        for face,owners in faces.items():
            if len(owners)==2:
                a,b,c=v[list(face)];normal=np.cross(b-a,c-a)
                signed=[float(np.dot(v[o]-a,normal)) for _,o in owners]
                same_side+=int(signed[0]*signed[1]>=0)
                internal+=int(labels[owners[0][0]]!=labels[owners[1][0]])
        rows.append(dict(case=name,cell_six_volumes=six.tolist(),unique_faces=len(faces),
                         exterior_faces=sum(len(x)==1 for x in faces.values()),shared_label_faces=internal,
                         maximum_face_incidence=max(map(len,faces.values())),same_side_pairs=same_side,
                         hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in zip(('vertices','cells','labels'),(v,t,labels))}))
    return rows


if __name__=='__main__':
    a,b=measure(),measure();gates=dict(full_repeat=a==b,finite_volumes=all(np.isfinite(r['cell_six_volumes']).all() for r in a))
    report=dict(rows=a,gates=gates,scope='Synthetic tetrahedron face inventory, not a boundary converter certificate.')
    (ROOT/'reports/tetra_interface_mechanism.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));raise SystemExit(0 if all(gates.values()) else 2)
