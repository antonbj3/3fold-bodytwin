"""Oriented shared material boundaries of a conforming tetrahedral mesh.

No repair, smoothing, compute initialization or implicit artifact writes. Exterior
label is zero; input cells have positive labels. Nonadjacent cell intersections
are not detected: callers must supply a conforming mesh, not a triangle soup.
"""
import numpy as np


def tetra_interfaces(vertices, cells, labels):
    v=np.asarray(vertices,dtype=np.float64)
    t=np.asarray(cells);regions=np.asarray(labels)
    if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all() or len(v)>np.iinfo(np.int32).max:
        raise ValueError('expected finite vertices within int32 index capacity')
    if t.ndim!=2 or t.shape[1]!=4 or not np.issubdtype(t.dtype,np.integer) or len(t)==0:
        raise ValueError('expected nonempty integral four-node cells')
    if t.min()<0 or t.max()>=len(v):
        raise ValueError('cell vertex index outside vertices')
    if regions.shape!=(len(t),) or not np.issubdtype(regions.dtype,np.integer) or regions.min()<1 or regions.max()>np.iinfo(np.int32).max:
        raise ValueError('expected one positive int32-range material per cell')
    t=t.astype(np.int64);regions=regions.astype(np.int32)
    p=v[t]
    six=np.einsum('ij,ij->i',p[:,1]-p[:,0],np.cross(p[:,2]-p[:,0],p[:,3]-p[:,0]))
    if not np.isfinite(six).all() or np.any(six==0):
        raise ValueError('zero or nonfinite cell volume')
    patterns=np.array([[1,2,3],[0,2,3],[0,1,3],[0,1,2]])
    faces=t[:,patterns].reshape(-1,3)
    unique,first,inverse,count=np.unique(np.sort(faces,axis=1),axis=0,return_index=True,return_inverse=True,return_counts=True)
    if np.any(count>2):
        raise ValueError('face has more than two incident cells')
    owners=np.repeat(np.arange(len(t)),4)
    opposites=t.reshape(-1)
    points=v[unique]
    normal=np.cross(points[:,1]-points[:,0],points[:,2]-points[:,0])
    side=np.einsum('ij,ij->i',normal,v[opposites[first]]-points[:,0])
    index_sum=np.zeros(len(unique),np.int64)
    np.add.at(index_sum,inverse,np.arange(len(faces)))
    shared=count==2;second=index_sum[shared]-first[shared]
    other_side=np.einsum('ij,ij->i',normal[shared],v[opposites[second]]-points[shared,0])
    if np.any(np.signbit(side[shared])==np.signbit(other_side)) or np.any(other_side==0):
        raise ValueError('shared face has cells on the same side')
    back=regions[owners[first]];front=np.zeros(len(unique),np.int32)
    front[shared]=regions[owners[second]]
    keep=front!=back
    triangles=unique.copy();flip=side>0
    triangles[flip]=triangles[flip][:,[0,2,1]]
    return dict(vertices=v.copy(),faces=triangles[keep].astype(np.int32),front=front[keep],back=back[keep])


if __name__=='__main__':
    from pathlib import Path
    import runpy
    import sys
    probe=Path(__file__).resolve().parents[3]/'probes/optics/tetra_interface_probe.py'
    sys.path.insert(0,str(probe.parent))
    runpy.run_path(str(probe),run_name='__main__')
