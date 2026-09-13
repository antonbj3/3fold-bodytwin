"""Public atlas -> strict mesh seam -> unchanged sibling field CPU backend."""
import hashlib
import json
import os
from pathlib import Path
import sys
import numpy as np
import trimesh

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.mesh_ingest_v1 import surface_mesh,tetrahedral_mesh
PITCH=2.0


def digest(*arrays):
    h=hashlib.sha256()
    for a in arrays:
        a=np.ascontiguousarray(a);h.update(str(a.dtype).encode());h.update(str(a.shape).encode());h.update(a.tobytes())
    return h.hexdigest()


def contracts():
    v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]])
    t=np.array([[0,1,2,3],[0,2,1,4]])
    m=tetrahedral_mesh(v,t,units='mm')
    cases=[lambda:tetrahedral_mesh(v,t,units='m'),
           lambda:tetrahedral_mesh(v,t[:,[1,0,2,3]],units='mm'),
           lambda:tetrahedral_mesh(v,np.vstack([t,t[0]]),units='mm'),
           lambda:tetrahedral_mesh(v,np.array([[0,1,2,2]]),units='mm'),
           lambda:tetrahedral_mesh(v,np.array([[0,1,2,8]]),units='mm'),
           lambda:tetrahedral_mesh(v*np.nan,t,units='mm'),
           lambda:surface_mesh(m.vertices_mm,m.faces[:-1],units='mm'),
           lambda:surface_mesh(m.vertices_mm,m.faces.astype(float),units='mm')]
    rejects=[]
    for call in cases:
        try:call();rejects.append(False)
        except ValueError:rejects.append(True)
    mv,mf=m.motion_arrays()
    return {'tet_volume_mm3':m.volume_mm3,'tet_faces':len(m.faces),'rejections':rejects,
            'motion_scale_exact':bool(np.array_equal(mv,m.vertices_mm*.001) and np.array_equal(mf,m.faces)),
            'tet_sha256':digest(m.vertices_mm,m.faces)}


def boundary_bound(v,f,origin,shape):
    band=np.zeros(shape,dtype=bool)
    for tri in v[f]:
        low=np.ceil((tri.min(0)-origin)/PITCH-.5001).astype(int)
        high=np.floor((tri.max(0)-origin)/PITCH+.5001).astype(int)
        low=np.maximum(low,0);high=np.minimum(high,np.asarray(shape)-1)
        if np.all(high>=low):band[tuple(slice(a,b+1) for a,b in zip(low,high))]=True
    return int(band.sum())*PITCH**3


def run(wp,field):
    raw=trimesh.load_mesh(ROOT/'examples/anatomy/left_kidney_atlas.obj',process=False)
    try:surface_mesh(raw.vertices,raw.faces,units='mm');strict_reject=False
    except ValueError:strict_reject=True
    mesh=surface_mesh(raw.vertices,raw.faces,units='mm',weld_exact=True)
    v,f=mesh.field_arrays();lo=v.min(0)-4*PITCH
    result=field.mesh_to_sdf_del(wp,v,f,PITCH,0.,lo,'cpu',grind='flagga',metod='raypar_vindning',returnera_falt=True)
    origin=lo+np.asarray(result['gmin'])*PITCH
    solid=result['solid_final'];sd=result['sd'];volume=float(solid.sum())*PITCH**3
    bound=boundary_bound(v,f,origin,solid.shape)
    gate=result['vattentathet'];contract=contracts()
    report={'contract':contract,'strict_original_rejected':strict_reject,'welded_vertices':mesh.welded_vertices,
            'face_coordinates_unchanged':raw.vertices[raw.faces].tobytes()==v[f].tobytes(),
            'mesh_volume_mm3':mesh.volume_mm3,'voxel_volume_mm3':volume,
            'volume_absolute_error_mm3':abs(volume-mesh.volume_mm3),
            'volume_relative_error':abs(volume-mesh.volume_mm3)/mesh.volume_mm3,
            'raster_bound_mm3':bound,'shape':list(solid.shape),'origin_mm':origin.tolist(),'pitch_mm':PITCH,
            'field_gate_ok':gate['status']=='OK','field_mesh_closed':bool(gate['vattentat']),
            'field_parity_disagreement':gate.get('paritet_divergens_frac'),
            'field_deep_parity_disagreement':gate.get('paritet_divergens_frac_djup'),
            'sparse_dense_mismatch':result['n_diff_pre_margin_gpu_vs_cpu'],
            'mesh_sha256':digest(v,f),'field_sha256':digest(sd,solid,origin),
            'finite_field':bool(np.isfinite(sd).all())}
    return report,{'vertices_mm':v,'faces':f,'sdf_mm':sd,'solid':solid,'origin_mm':origin,'pitch_mm':np.array(PITCH)}


def main():
    siblings=Path(os.environ.get('THREEFOLD_ROOT',ROOT.parent))
    sys.path.insert(0,str(siblings/'3fold-field-engine'/'src'/'field_engine'))
    import faltkarna_v1_mesh_to_sdf as field
    os.environ['CUDA_VISIBLE_DEVICES']=''
    import warp as wp
    wp.init()
    results=[run(wp,field) for _ in range(2)];legs=[x[0] for x in results]
    gates={'two_full_results_identical':legs[0]==legs[1] and all(np.array_equal(results[0][1][k],results[1][1][k]) for k in results[0][1]),
           'original_rejected_and_exact_weld':all(r['strict_original_rejected'] and r['welded_vertices']==4 and r['face_coordinates_unchanged'] for r in legs),
           'tet_contract':all(abs(r['contract']['tet_volume_mm3']-1/3)<=1e-12/3 and r['contract']['tet_faces']==6 and all(r['contract']['rejections']) and r['contract']['motion_scale_exact'] for r in legs),
           'frozen_field_gate':all(r['field_gate_ok'] and r['field_mesh_closed'] and r['finite_field'] for r in legs),
           'sparse_dense_exact':all(r['sparse_dense_mismatch']==0 for r in legs),
           'volume_band':all(r['volume_relative_error']<=.05 and r['volume_absolute_error_mm3']<=r['raster_bound_mm3'] for r in legs)}
    out=ROOT/'examples/anatomy'
    report={'legs':legs,'gates':gates,'backend':'CPU','scope':'public atlas geometry seam; no clinical or motion contact integration claim'}
    (out/'mesh_to_field_result.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(out/'mesh_to_field_arrays.npz',**results[0][1])
    print(json.dumps(report),flush=True)
    return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
