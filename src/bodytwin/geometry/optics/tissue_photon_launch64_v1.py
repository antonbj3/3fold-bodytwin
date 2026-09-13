"""Prepared photon transport with a measured 64-thread launch configuration."""
import numpy as np
import warp as wp
from bodytwin.geometry.optics.tissue_photon_segmented_v1 import _trace, PACKET
from bodytwin.geometry.optics.tissue_photon_prepared_v1 import PreparedPhotonScene as PreparedBaseline


class PreparedPhotonScene(PreparedBaseline):
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
            initial_region,seed,int(reflect),wp.vec3(*origin),float(pitch),*shape,float(window_ns),absorption,terminal,counters,paths,exits],device=device,block_dim=64)
        return {'absorption':absorption.numpy().reshape(shape),'terminal':terminal.numpy(),'counters':counters.numpy(),'paths':paths.numpy(),'exits':exits.numpy()}


if __name__=='__main__':
    import subprocess
    import sys
    from pathlib import Path
    root=Path(__file__).resolve().parents[4]
    raise SystemExit(subprocess.call([sys.executable,str(root/'probes/optics/photon_segmented_transport_probe.py'),'--launch64',*sys.argv[1:]]))
