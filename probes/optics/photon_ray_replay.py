"""CPU replay of rounded diagnostic rays; not a full-path equivalence claim."""
import os
os.environ['CUDA_VISIBLE_DEVICES']=''
import json
from pathlib import Path
import numpy as np
import trimesh
import warp as wp
ROOT=Path(__file__).resolve().parents[2]


@wp.kernel
def queries(mesh:wp.uint64,positions:wp.array(dtype=wp.vec3),directions:wp.array(dtype=wp.vec3),limits:wp.array(dtype=float),out:wp.array2d(dtype=float)):
    i=wp.tid()
    q=wp.mesh_query_ray(mesh,positions[i],directions[i],limits[i])
    if q.result:out[i,0]=1.0;out[i,1]=q.t
    q=wp.mesh_query_ray(mesh,positions[i],directions[i],1.e20)
    if q.result:out[i,2]=1.0;out[i,3]=q.t


def exact_hit(p,d,vertices,faces):
    hits=[]
    for i,f in enumerate(faces):
        a,b,c=vertices[f];e1=b-a;e2=c-a;h=np.cross(d,e2);det=np.dot(e1,h)
        if det==0:continue
        q=p-a;u=np.dot(q,h)/det
        if not 0<=u<=1:continue
        k=np.cross(q,e1);v=np.dot(d,k)/det
        if v<0 or u+v>1:continue
        t=float(np.dot(e2,k)/det)
        if t>=0:hits.append((t,i))
    return min(hits) if hits else (None,None)


def main():
    wp.init();r=json.loads((ROOT/'reports/photon_boundary_failure_l4.json').read_text())
    samples=[s for c in r['cases'] if c['case']=='slab' for g in c['legs'][0]['groups'] if g['code'] in (2,3) for s in g['samples']]
    m=trimesh.creation.box(extents=[100]*3);m.apply_translation([50]*3)
    p=np.asarray([s['position'] for s in samples],np.float32);d=np.asarray([s['direction'] for s in samples],np.float32)
    limits=np.asarray([s['distance'] for s in samples],np.float32)
    mesh=wp.Mesh(points=wp.array(m.vertices,dtype=wp.vec3,device='cpu'),indices=wp.array(m.faces.ravel(),dtype=int,device='cpu'))
    legs=[]
    for _ in range(2):
        out=wp.zeros((len(samples),4),dtype=float,device='cpu')
        wp.launch(queries,len(samples),inputs=[mesh.id,wp.array(p,dtype=wp.vec3,device='cpu'),wp.array(d,dtype=wp.vec3,device='cpu'),wp.array(limits,dtype=float,device='cpu'),out],device='cpu')
        a=out.numpy();rows=[]
        for i,s in enumerate(samples):
            t,face=exact_hit(p[i].astype(float),d[i].astype(float),m.vertices,m.faces)
            rows.append({'photon':s['photon'],'position':p[i].astype(float).tolist(),'limit':float(limits[i]),
                'bounded_hit':bool(a[i,0]),'bounded_distance':float(a[i,1]),'unbounded_hit':bool(a[i,2]),'unbounded_distance':float(a[i,3]),
                'double_triangle_distance':t,'double_face':int(face) if face is not None else None,
                'double_within_limit':bool(t is not None and t<=float(limits[i]))})
        legs.append(rows)
    gates={'complete_repeat_exact':legs[0]==legs[1],'all_rays_have_double_hit':all(r['double_face'] is not None for r in legs[0])}
    report={'legs':legs,'gates':gates,'scope':'CPU replay of diagnostic coordinates truncated at1e-9; this is not exact reconstruction of every original float32 coordinate or a GPU path claim'}
    (ROOT/'reports/photon_ray_replay.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True);return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
