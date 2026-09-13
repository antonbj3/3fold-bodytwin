"""Measure Warp mesh point accessor indexing against the source triangle corners."""
import os
os.environ['CUDA_VISIBLE_DEVICES']=''
import json
from pathlib import Path
import numpy as np
import trimesh
import warp as wp
ROOT=Path(__file__).resolve().parents[2]


@wp.kernel
def read_corners(mesh:wp.uint64,direct:wp.array(dtype=wp.vec3),indirect:wp.array(dtype=wp.vec3)):
    i=wp.tid()
    direct[i]=wp.mesh_get_point(mesh,i)
    indirect[i]=wp.mesh_get_point(mesh,wp.mesh_get_index(mesh,i))


def main():
    wp.init();m=trimesh.creation.box(extents=[60]*3);m.apply_translation([30]*3)
    mesh=wp.Mesh(points=wp.array(m.vertices,dtype=wp.vec3,device='cpu'),indices=wp.array(m.faces.ravel(),dtype=int,device='cpu'))
    expected=np.asarray(m.vertices[m.faces.ravel()],dtype=np.float32);legs=[]
    for _ in range(2):
        direct=wp.zeros(36,dtype=wp.vec3,device='cpu');indirect=wp.zeros(36,dtype=wp.vec3,device='cpu')
        wp.launch(read_corners,36,inputs=[mesh.id,direct,indirect],device='cpu')
        a=direct.numpy();b=indirect.numpy()
        legs.append({'direct':a.tolist(),'extra_indirection':b.tolist(),'direct_mismatched_corners':int(np.any(a!=expected,axis=1).sum()),
                     'extra_indirection_mismatched_corners':int(np.any(b!=expected,axis=1).sum()),'direct_bytes_exact':a.tobytes()==expected.tobytes()})
    gates={'complete_repeat_exact':legs[0]==legs[1],'direct_source_corners_exact':all(r['direct_bytes_exact'] for r in legs)}
    report={'legs':legs,'gates':gates,'source_contract':'Warp1.13 mesh_get_point applies mesh.indices internally; argument is a face-corner index'}
    (ROOT/'reports/photon_vertex_accessor.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'gates':gates,'direct_mismatches':legs[0]['direct_mismatched_corners'],'extra_indirection_mismatches':legs[0]['extra_indirection_mismatched_corners']}))
    return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
