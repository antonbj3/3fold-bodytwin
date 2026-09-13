"""Prepared straight-ray traversal of a conforming closed tetrahedral partition.

First exterior exit only. Exact singular launches/exits are refused, not shifted.
Nonadjacent cell intersections are outside this contract. No implicit file writes.
"""
import numpy as np
from bodytwin.geometry.tetra_interfaces_v1 import tetra_interfaces


class TetraRayWalk:
    def __init__(self,vertices,cells,labels):
        mesh=tetra_interfaces(vertices,cells,labels)
        self.vertices=mesh['vertices'];self.cells=np.asarray(cells,dtype=np.int64).copy();self.labels=np.asarray(labels,dtype=np.int32).copy()
        v=self.vertices;t=self.cells
        outer=tetra_interfaces(v,t,np.ones(len(t),np.int32))['faces']
        edges=np.sort(outer[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1)
        _,incidence=np.unique(edges,axis=0,return_counts=True)
        if np.any(incidence!=2):raise ValueError('Nonmanifold outer partition boundary')
        patterns=np.array([[1,2,3],[0,2,3],[0,1,3],[0,1,2]])
        faces=np.sort(t[:,patterns].reshape(-1,3),axis=1)
        _,first,inverse,counts=np.unique(faces,axis=0,return_index=True,return_inverse=True,return_counts=True)
        sums=np.zeros(len(counts),np.int64);np.add.at(sums,inverse,np.arange(len(faces)))
        self.neighbors=np.full((len(t),4),-1,np.int64);self.entry_faces=np.full((len(t),4),-1,np.int64)
        for left,right in zip(first[counts==2],sums[counts==2]-first[counts==2]):
            a,i=divmod(int(left),4);b,j=divmod(int(right),4)
            self.neighbors[a,i]=b;self.neighbors[b,j]=a;self.entry_faces[a,i]=j;self.entry_faces[b,j]=i
        p=v[t];self.origins=p[:,0].copy();self.inverse=np.linalg.inv(np.stack((p[:,1]-p[:,0],p[:,2]-p[:,0],p[:,3]-p[:,0]),axis=2))
        for array in (self.vertices,self.cells,self.labels,self.neighbors,self.entry_faces,self.origins,self.inverse):array.flags.writeable=False

    def trace(self,origin,direction):
        origin=np.asarray(origin,dtype=np.float64);direction=np.asarray(direction,dtype=np.float64)
        if origin.shape!=(3,) or direction.shape!=(3,) or not np.isfinite([origin,direction]).all():raise ValueError('Finite three-component ray required')
        length=np.linalg.norm(direction)
        if not np.isfinite(length) or length==0:raise ValueError('Nonzero finite ray direction required')
        direction=direction/length
        local=np.einsum('nij,nj->ni',self.inverse,origin-self.origins)
        bary=np.column_stack((1-local.sum(axis=1),local));inside=np.flatnonzero(np.min(bary,axis=1)>0)
        if len(inside)!=1:raise ValueError('Source must be strictly inside exactly one cell')
        cell=int(inside[0]);entry=-1;parameter=0.;records=[];visited=set()
        for _ in range(len(self.cells)+1):
            if cell in visited:raise ValueError('Traversal revisited a convex cell')
            visited.add(cell)
            q=self.inverse[cell]@direction;delta=np.r_[-q.sum(),q]
            candidates=[]
            for face in range(4):
                if face!=entry and delta[face]<0:
                    value=-bary[cell,face]/delta[face]
                    if value<parameter:raise ValueError('Negative cell interval')
                    candidates.append((float(value),face))
            if not candidates:raise ValueError('No forward exit from bounded cell')
            candidates.sort();exit_parameter,face=candidates[0]
            if len(candidates)>1 and candidates[1][0]==exit_parameter:raise ValueError('Ambiguous edge or vertex exit')
            if exit_parameter<=parameter:raise ValueError('Zero-length cell interval')
            records.append((cell,parameter,exit_parameter))
            neighbor=int(self.neighbors[cell,face])
            if neighbor<0:return np.asarray(records,dtype=np.float64)
            entry=int(self.entry_faces[cell,face]);cell=neighbor;parameter=exit_parameter
        raise ValueError('Cell traversal bound exceeded')


if __name__=='__main__':
    import runpy
    from pathlib import Path
    runpy.run_path(str(Path(__file__).resolve().parents[3]/'probes/optics/tetra_ray_walk_probe.py'),run_name='__main__')
