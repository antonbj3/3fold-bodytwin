"""Strict array seam for closed surfaces and positive tetrahedral complexes.

Coordinates must explicitly use mm. Optional welding merges identical float64
coordinate bytes only. No self-intersection or anatomical-validity certificate.
"""
from dataclasses import dataclass
import numpy as np
import trimesh


@dataclass(frozen=True)
class MeshArrays:
    vertices_mm: np.ndarray
    faces: np.ndarray
    welded_vertices: int = 0

    def field_arrays(self):
        return self.vertices_mm.copy(),self.faces.copy()

    def motion_arrays(self):
        """Explicit SI conversion; does not integrate a triangle contact backend."""
        return self.vertices_mm*.001,self.faces.copy()

    @property
    def volume_mm3(self):
        return float(trimesh.Trimesh(self.vertices_mm,self.faces,process=False).volume)


def _arrays(vertices,indices,width,units):
    if units!='mm':raise ValueError('Explicit mm coordinates required')
    v=np.asarray(vertices,dtype=np.float64)
    i=np.asarray(indices)
    if v.ndim!=2 or v.shape[1]!=3 or len(v)<4 or not np.isfinite(v).all():
        raise ValueError('Finite (N,3) vertices required')
    if i.ndim!=2 or i.shape[1]!=width or len(i)==0 or not np.issubdtype(i.dtype,np.integer):
        raise ValueError('Nonempty integer connectivity required')
    if i.min()<0 or i.max()>=len(v) or i.max()>np.iinfo(np.int32).max:
        raise ValueError('Connectivity outside vertex range')
    return np.ascontiguousarray(v).copy(),np.ascontiguousarray(i,dtype=np.int32)


def surface_mesh(vertices,faces,*,units,weld_exact=False):
    v,f=_arrays(vertices,faces,3,units);count=0
    if weld_exact:
        _,first,inverse=np.unique(v.view(np.dtype((np.void,24))).ravel(),return_index=True,return_inverse=True)
        count=len(v)-len(first);v=v[first].copy();f=inverse[f].astype(np.int32)
    if np.any(np.diff(np.sort(f,axis=1),axis=1)==0):raise ValueError('Repeated face vertex')
    if len(np.unique(np.sort(f,axis=1),axis=0))!=len(f):raise ValueError('Duplicate face')
    m=trimesh.Trimesh(v,f,process=False)
    if not m.is_watertight or not m.is_winding_consistent:raise ValueError('Closed consistently oriented surface required')
    if np.any(m.area_faces<=0) or not np.isfinite(m.volume) or m.volume<=0:
        raise ValueError('Positive volume and nondegenerate triangles required')
    v.setflags(write=False);f.setflags(write=False)
    return MeshArrays(v,f,count)


def tetrahedral_mesh(vertices,tetrahedra,*,units):
    v,t=_arrays(vertices,tetrahedra,4,units)
    if len(np.unique(np.sort(t,axis=1),axis=0))!=len(t):raise ValueError('Duplicate tetrahedron')
    p=v[t];det=np.linalg.det(np.stack([p[:,1]-p[:,0],p[:,2]-p[:,0],p[:,3]-p[:,0]],axis=2))
    if not np.isfinite(det).all() or np.any(det<=0):raise ValueError('Positive tetrahedron orientation required')
    oriented=t[:,[[1,2,3],[0,3,2],[0,1,3],[0,2,1]]].reshape(-1,3)
    keys=np.sort(oriented,axis=1)
    _,first,inverse,counts=np.unique(keys,axis=0,return_index=True,return_inverse=True,return_counts=True)
    if np.any(counts>2):raise ValueError('Nonmanifold tetrahedron face')
    # Opposite parity on each interior face is required; no overlapping-side gluing.
    parity=np.where(((oriented[:,0]>oriented[:,1]).astype(int)+(oriented[:,0]>oriented[:,2])+(oriented[:,1]>oriented[:,2]))%2, -1,1)
    orientation_sum=np.bincount(inverse,weights=parity,minlength=len(counts))
    if np.any(orientation_sum[counts==2]!=0):raise ValueError('Inconsistent interior face orientation')
    return surface_mesh(v,oriented[first[counts==1]],units=units)
