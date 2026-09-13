"""Distinguish volume-partition topology from per-material boundary topology."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.tetra_interfaces_v1 import tetra_interfaces


def edge_histogram(faces,nvertices):
    edges=faces[:,[0,1,1,2,2,0]].reshape(-1,2);edges=np.sort(edges,axis=1)
    _,count=np.unique(edges[:,0].astype(np.int64)*nvertices+edges[:,1],return_counts=True)
    values,counts=np.unique(count,return_counts=True)
    return {str(int(k)):int(v) for k,v in zip(values,counts)}


def measure(v,t,labels):
    mesh=tetra_interfaces(v,t,labels)
    all_domain=tetra_interfaces(v,t,np.ones(len(t),np.int32))
    outer=all_domain['faces'];points=v[t]
    six=np.einsum('ij,ij->i',points[:,1]-points[:,0],np.cross(points[:,2]-points[:,0],points[:,3]-points[:,0]))
    tri=v[outer];surface_volume=float(np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2])).sum()/6)
    volume=float(np.abs(six).sum()/6)
    faces=np.concatenate([np.sort(np.delete(t,i,axis=1),axis=1) for i in range(4)])
    owners=np.tile(np.arange(len(t)),4)
    _,first,inverse,counts=np.unique(faces,axis=0,return_index=True,return_inverse=True,return_counts=True)
    sums=np.zeros(len(counts),np.int64);np.add.at(sums,inverse,np.arange(len(faces)))
    shared=counts==2;a=owners[first[shared]];b=owners[sums[shared]-first[shared]]
    graph=coo_matrix((np.ones(2*len(a),np.int8),(np.r_[a,b],np.r_[b,a])),shape=(len(t),len(t))).tocsr()
    parts=[]
    for label in np.unique(labels):
        f=np.concatenate((mesh['faces'][mesh['back']==label],mesh['faces'][mesh['front']==label][:,::-1]))
        parts.append(dict(label=int(label),edge_incidence=edge_histogram(f,len(v))))
    return dict(cells=len(t),exterior_faces=len(outer),maximum_face_incidence=int(counts.max()),cell_components=int(connected_components(graph,directed=False,return_labels=False)),total_cell_volume=volume,outer_signed_volume=surface_volume,relative_volume_error=abs(surface_volume-volume)/volume,outer_edge_incidence=edge_histogram(outer,len(v)),regions=parts,source_hashes={k:hashlib.sha256(np.asarray(a).tobytes()).hexdigest() for k,a in [('vertices',v),('cells',t),('labels',labels)]})


def fixtures():
    v=np.array([[0,0,-1],[0,0,1],[1,0,0],[0,1,0],[-1,0,0],[0,-1,0]],float)
    t=np.array([[0,1,2+i,2+(i+1)%4] for i in range(4)],np.int32)
    return [('one_material',v,t,np.ones(4,np.int32)),('alternating',v,t,np.array([1,2,1,2],np.int32))]


def main():
    first=[dict(case=name,**measure(v,t,labels)) for name,v,t,labels in fixtures()]
    second=[dict(case=name,**measure(v,t,labels)) for name,v,t,labels in fixtures()]
    gates=dict(repeat=first==second,closed_partition=all(set(r['outer_edge_incidence'])=={'2'} and r['maximum_face_incidence']==2 and r['cell_components']==1 and r['relative_volume_error']<=1e-10 for r in first),material_defect_detected=all('4' in r['edge_incidence'] for r in first[1]['regions']))
    report=dict(rows=first,gates=gates,scope='Synthetic partition-versus-boundary observer, not repair or transport acceptance.')
    (ROOT/'reports/tetra_partition_mechanism.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
