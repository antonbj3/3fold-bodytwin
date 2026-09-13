"""Integer-energy attenuation over complete certified cell traces, not scattering."""
import numpy as np
import warp as wp
PACKET=1<<30


@wp.kernel
def _score(cells:wp.array2d(dtype=wp.int32),intervals:wp.array3d(dtype=wp.float64),counts:wp.array(dtype=wp.int32),status:wp.array(dtype=wp.int32),labels:wp.array(dtype=wp.int32),mu:wp.array(dtype=wp.float64),dose:wp.array(dtype=wp.int64),terminal:wp.array2d(dtype=wp.int64)):
    ray=wp.tid()
    if status[ray]!=0:
        terminal[ray,2]=wp.int64(1073741824)
    else:
        continuous=wp.float64(1073741824.0);weight=wp.int64(1073741824)
        for slot in range(counts[ray]):
            cell=cells[ray,slot];step=mu[labels[cell]]*(intervals[ray,slot,1]-intervals[ray,slot,0])
            continuous*=wp.exp(-step);after=wp.int64(wp.floor(continuous+wp.float64(.5)))
            wp.atomic_add(dose,cell,weight-after);weight=after
        terminal[ray,0]=wp.int64(1073741824)-weight;terminal[ray,1]=weight


def attenuate(prepared,traversal,mu,*,device='cuda:0'):
    """Return owned dose/terminal arrays; failed traces remain fully unscored."""
    from bodytwin.geometry.tetra_ray_walk_v1 import TetraRayWalk
    if not isinstance(prepared,TetraRayWalk):raise ValueError('Prepared partition required')
    cells=np.asarray(traversal['cells']);intervals=np.asarray(traversal['intervals'],dtype=np.float64);counts=np.asarray(traversal['counts']);status=np.asarray(traversal['status']);mu=np.asarray(mu,dtype=np.float64)
    if cells.ndim!=2 or intervals.shape!=cells.shape+(2,) or counts.shape!=(len(cells),) or status.shape!=counts.shape or any(a.dtype.kind not in 'iu' for a in (cells,counts,status)) or len(cells)==0 or len(cells)*PACKET>np.iinfo(np.int64).max:raise ValueError('Invalid trace array contract')
    if np.any(counts<0) or np.any(counts>cells.shape[1]) or np.any(status<0) or np.any(status>4):raise ValueError('Invalid counts/status')
    if mu.ndim!=1 or len(mu)<=int(prepared.labels.max()) or not np.isfinite(mu).all() or np.any(mu<0):raise ValueError('Finite nonnegative attenuation per label required')
    for ray in np.flatnonzero(status==0):
        n=int(counts[ray]);ids=cells[ray,:n];d=intervals[ray,:n]
        if n==0 or ids.min()<0 or ids.max()>=len(prepared.cells) or not np.isfinite(d).all() or d[0,0]!=0 or np.any(d[:,1]<=d[:,0]) or np.any(d[1:,0]!=d[:-1,1]):raise ValueError('Invalid complete cell intervals')
    dose=wp.zeros(len(prepared.cells),dtype=wp.int64,device=device);terminal=wp.zeros((len(cells),3),dtype=wp.int64,device=device)
    wp.launch(_score,len(cells),inputs=[wp.array(cells.astype(np.int32),dtype=wp.int32,device=device),wp.array(intervals,dtype=wp.float64,device=device),wp.array(counts.astype(np.int32),dtype=wp.int32,device=device),wp.array(status.astype(np.int32),dtype=wp.int32,device=device),wp.array(prepared.labels.copy(),dtype=wp.int32,device=device),wp.array(mu,dtype=wp.float64,device=device),dose,terminal],device=device)
    return dict(dose=dose.numpy(),terminal=terminal.numpy())


if __name__=='__main__':
    import runpy
    from pathlib import Path
    runpy.run_path(str(Path(__file__).resolve().parents[4]/'probes/optics/tetra_attenuation_probe.py'),run_name='__main__')
