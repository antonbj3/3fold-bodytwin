"""Separate CUDA straight-ray cell walk with explicit capacity/error status."""
import numpy as np
import warp as wp
from bodytwin.geometry.tetra_ray_walk_v1 import TetraRayWalk


@wp.kernel
def _walk(inverse:wp.array3d(dtype=wp.float64),bary:wp.array2d(dtype=wp.float64),
          neighbors:wp.array2d(dtype=wp.int32),entry_faces:wp.array2d(dtype=wp.int32),
          directions:wp.array(dtype=wp.vec3d),initial:int,capacity:int,
          cells:wp.array2d(dtype=wp.int32),intervals:wp.array3d(dtype=wp.float64),
          counts:wp.array(dtype=wp.int32),status:wp.array(dtype=wp.int32)):
    ray=wp.tid();direction=directions[ray];cell=initial;entry=int(-1);parameter=wp.float64(0.0)
    finished=int(0)
    for slot in range(capacity):
        q0=inverse[cell,0,0]*direction[0]+inverse[cell,0,1]*direction[1]+inverse[cell,0,2]*direction[2]
        q1=inverse[cell,1,0]*direction[0]+inverse[cell,1,1]*direction[1]+inverse[cell,1,2]*direction[2]
        q2=inverse[cell,2,0]*direction[0]+inverse[cell,2,1]*direction[1]+inverse[cell,2,2]*direction[2]
        delta=wp.vec4d(-(q0+q1+q2),q0,q1,q2)
        best=wp.float64(1.0e300);face=int(-1);ties=int(0);negative=int(0)
        for i in range(4):
            if i!=entry and delta[i]<wp.float64(0.0):
                value=-bary[cell,i]/delta[i]
                if value<parameter:negative=1
                if value<best:best=value;face=i;ties=1
                elif value==best:ties+=1
        if negative!=0:status[ray]=1;finished=1;break
        if ties>1:status[ray]=2;finished=1;break
        if face<0 or best<=parameter:status[ray]=3;finished=1;break
        cells[ray,slot]=cell;intervals[ray,slot,0]=parameter;intervals[ray,slot,1]=best;counts[ray]=slot+1
        neighbor=neighbors[cell,face]
        if neighbor<0:finished=1;break
        entry=entry_faces[cell,face];cell=neighbor;parameter=best
    if finished==0:status[ray]=4


def trace(prepared,origin,directions,*,capacity=32,device='cuda:0'):
    """Return owned cells/intervals/counts/status arrays; Warp manages kernel caching.

    Status0 complete,1 negative interval,2 ambiguous tie,3 missing/nonpositive
    exit,4 capacity exceeded. A failed ray retains only its explicitly counted
    prefix and is never a completed transport result.
    """
    if not isinstance(prepared,TetraRayWalk):raise ValueError('Prepared closed cell partition required')
    origin=np.asarray(origin,dtype=np.float64);directions=np.asarray(directions,dtype=np.float64)
    if origin.shape!=(3,) or not np.isfinite(origin).all() or directions.ndim!=2 or directions.shape[1]!=3 or len(directions)==0 or not np.isfinite(directions).all():raise ValueError('Finite source and nonempty direction rows required')
    lengths=np.array([np.linalg.norm(d) for d in directions])
    if np.any(lengths==0) or not np.isfinite(lengths).all():raise ValueError('Nonzero finite directions required')
    if type(capacity) is not int or not 1<=capacity<=4096 or len(prepared.cells)>np.iinfo(np.int32).max:raise ValueError('Unsupported capacity or cell count')
    directions=np.ascontiguousarray(directions/lengths[:,None])
    local=np.einsum('nij,nj->ni',prepared.inverse,origin-prepared.origins)
    bary=np.column_stack((1-local.sum(axis=1),local));inside=np.flatnonzero(np.min(bary,axis=1)>0)
    if len(inside)!=1:raise ValueError('Source must be strictly inside exactly one cell')
    cells=wp.full((len(directions),capacity),-1,dtype=wp.int32,device=device)
    distances=wp.zeros((len(directions),capacity,2),dtype=wp.float64,device=device)
    counts=wp.zeros(len(directions),dtype=wp.int32,device=device);status=wp.zeros(len(directions),dtype=wp.int32,device=device)
    wp.launch(_walk,len(directions),inputs=[wp.array(prepared.inverse.copy(),dtype=wp.float64,device=device),wp.array(bary,dtype=wp.float64,device=device),wp.array(prepared.neighbors.astype(np.int32),dtype=wp.int32,device=device),wp.array(prepared.entry_faces.astype(np.int32),dtype=wp.int32,device=device),wp.array(directions,dtype=wp.vec3d,device=device),int(inside[0]),capacity,cells,distances,counts,status],device=device)
    return dict(cells=cells.numpy(),intervals=distances.numpy(),counts=counts.numpy(),status=status.numpy())


if __name__=='__main__':
    import runpy
    from pathlib import Path
    runpy.run_path(str(Path(__file__).resolve().parents[3]/'probes/optics/tetra_ray_cuda_probe.py'),run_name='__main__')
