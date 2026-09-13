"""Two-region synthetic box seam with one shared oriented interface."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import trimesh


def layered_box(side=60.,split=10.):
    if not np.isfinite([side,split]).all() or not 0<split<side:
        raise ValueError('Require a finite positive side and interior split')
    vertices=[];faces=[];front=[];back=[];offset=0
    for region,lo,hi in ((1,0.,split),(2,split,side)):
        mesh=trimesh.creation.box(extents=[side,side,hi-lo]);mesh.apply_translation([side/2,side/2,(lo+hi)/2])
        interface=np.all(mesh.triangles[:,:,2]==split,axis=1)
        keep=np.ones(len(mesh.faces),dtype=bool) if region==1 else ~interface
        vertices.append(mesh.vertices);faces.append(mesh.faces[keep]+offset)
        labels=np.zeros(len(mesh.faces),dtype=np.int32)
        if region==1:labels[interface]=2
        front.append(labels[keep]);back.append(np.full(int(keep.sum()),region,dtype=np.int32));offset+=len(mesh.vertices)
    v=np.concatenate(vertices);f=np.concatenate(faces)
    v,inverse=np.unique(v,axis=0,return_inverse=True)
    return {'vertices':v,'faces':inverse[f].astype(np.int32),'front':np.concatenate(front),'back':np.concatenate(back)}


def main():
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
    from bodytwin.geometry.mesh_ingest_v1 import surface_mesh
    legs=[]
    for _ in range(2):
        out=layered_box();v=out['vertices'];f=out['faces'];volumes=[]
        for region in (1,2):
            outward=np.concatenate([f[out['back']==region],f[out['front']==region][:,::-1]])
            m=surface_mesh(v,outward,units='mm');volumes.append(m.volume_mm3)
        interface=(out['front']==2)&(out['back']==1)
        mesh=trimesh.Trimesh(v,f,process=False)
        legs.append({'hashes':{k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in out.items()},
            'vertices':len(v),'faces':len(f),'region_volumes':volumes,'strict_region_validation':True,
            'interface_faces':int(interface.sum()),'interface_area':float(mesh.area_faces[interface].sum())})
    gates={'full_array_repeat_exact':legs[0]==legs[1],
        'closed_oriented_region_boundaries':all(r['strict_region_validation'] for r in legs),
        'analytic_volumes_exact':all(r['region_volumes']==[36000.,180000.] for r in legs),
        'single_shared_interface':all(r['vertices']==12 and r['faces']==22 and r['interface_faces']==2 and r['interface_area']==3600. for r in legs)}
    report={'legs':legs,'gates':gates,'scope':'synthetic two-region topology/volume seam, not transport certification'}
    (Path(__file__).resolve().parents[3]/'reports/layered_box_v1.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
