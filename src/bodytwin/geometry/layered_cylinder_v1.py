"""Synthetic nested curved two-region fixture, not measured anatomy."""
import numpy as np
import trimesh


def layered_cylinder():
    outer=trimesh.creation.cylinder(radius=30.,height=120.,sections=64)
    inner=trimesh.creation.cylinder(radius=27.,height=114.,sections=64)
    return dict(vertices=np.concatenate((outer.vertices,inner.vertices))+[32.,32.,62.],
                faces=np.concatenate((outer.faces,inner.faces+len(outer.vertices))).astype(np.int32),
                front=np.concatenate((np.zeros(len(outer.faces),np.int32),np.ones(len(inner.faces),np.int32))),
                back=np.concatenate((np.ones(len(outer.faces),np.int32),np.full(len(inner.faces),2,np.int32))))


if __name__=='__main__':
    from pathlib import Path
    import runpy
    runpy.run_path(str(Path(__file__).resolve().parents[3]/'probes/optics/photon_curved_layer_probe.py'),run_name='__main__')
