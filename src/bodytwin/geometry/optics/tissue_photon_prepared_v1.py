"""Owned device geometry for repeated exact segmented photon transport.

Construction snapshots geometry and optical properties. Each launch allocates fresh
outputs and uses the frozen transport kernel. No implicit artifact writes.
"""
import numpy as np
import warp as wp
from bodytwin.geometry.optics.tissue_photon_segmented_v1 import _trace, PACKET


class PreparedPhotonScene:
    def __init__(self, vertices, faces, properties, face_front, face_back, *, device='cuda:0'):
        v=np.asarray(vertices,dtype=np.float32);f=np.asarray(faces,dtype=np.int32)
        prop=np.asarray(properties,dtype=np.float32);front=np.asarray(face_front,dtype=np.int32);back=np.asarray(face_back,dtype=np.int32)
        if v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3:raise ValueError('Triangle arrays required')
        if not np.isfinite(v).all() or f.min()<0 or f.max()>=len(v):raise ValueError('Invalid geometry')
        if prop.ndim!=2 or prop.shape[1]!=4 or not np.isfinite(prop).all():raise ValueError('Optical rows must be mua,mus,g,n')
        if np.any(prop[:,:2]<0) or np.any(np.abs(prop[:,2])>1) or np.any(prop[:,3]<=0):raise ValueError('Invalid optical properties')
        if len(front)!=len(f) or len(back)!=len(f) or min(front.min(),back.min())<0 or max(front.max(),back.max())>=len(prop):raise ValueError('Invalid face region labels')
        normal=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]])
        lengths=np.linalg.norm(normal,axis=1)
        if np.any(lengths==0):raise ValueError('Degenerate triangle')
        normal/=lengths[:,None]
        plane_point=v[f[:,0]].astype(np.float64)
        plane_normal=np.cross(v[f[:,1]].astype(np.float64)-plane_point,v[f[:,2]].astype(np.float64)-plane_point)
        plane_normal/=np.linalg.norm(plane_normal,axis=1)[:,None]
        self._device=device;self._regions=len(prop)
        self._mesh=wp.Mesh(points=wp.array(v,dtype=wp.vec3,device=device),indices=wp.array(f.ravel(),dtype=int,device=device))
        self._front=wp.array(front,dtype=int,device=device);self._back=wp.array(back,dtype=int,device=device)
        self._normals=wp.array(normal,dtype=wp.vec3,device=device)
        self._plane_points=wp.array(plane_point,dtype=wp.vec3d,device=device)
        self._plane_normals=wp.array(plane_normal,dtype=wp.vec3d,device=device)
        self._properties=wp.array(prop,dtype=float,device=device)

    def simulate(self, *, source, direction, initial_region, origin, shape, pitch=1.,
                 photons=1000000, seed=20260912, reflect=False, window_ns=5.):
        if not 0<initial_region<self._regions or photons<1 or photons*PACKET>np.iinfo(np.int64).max:raise ValueError('Invalid launch/energy count')
        shape=tuple(int(x) for x in shape)
        if len(shape)!=3 or min(shape)<=0 or pitch<=0 or window_ns<=0:raise ValueError('Positive grid/window required')
        direction=np.array(direction,dtype=float,copy=True);direction/=np.linalg.norm(direction)
        device=self._device
        absorption=wp.zeros(int(np.prod(shape)),dtype=wp.int64,device=device)
        terminal=wp.zeros((photons,4),dtype=wp.int64,device=device);counters=wp.zeros(3,dtype=wp.int64,device=device)
        paths=wp.zeros((photons,self._regions-1),dtype=wp.float64,device=device)
        exits=wp.zeros((photons,7),dtype=wp.float64,device=device)
        wp.launch(_trace,photons,inputs=[self._mesh.id,self._front,self._back,
            self._normals,self._plane_points,self._plane_normals,self._properties,wp.vec3(*source),wp.vec3(*direction),
            initial_region,seed,int(reflect),wp.vec3(*origin),float(pitch),*shape,float(window_ns),absorption,terminal,counters,paths,exits],device=device)
        return {'absorption':absorption.numpy().reshape(shape),'terminal':terminal.numpy(),'counters':counters.numpy(),'paths':paths.numpy(),'exits':exits.numpy()}


if __name__=='__main__':
    import subprocess
    import sys
    from pathlib import Path
    root=Path(__file__).resolve().parents[4]
    raise SystemExit(subprocess.call([sys.executable,str(root/'probes/optics/photon_segmented_transport_probe.py'),'--prepared',*sys.argv[1:]]))
