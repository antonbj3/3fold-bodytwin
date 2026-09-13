"""Synthetic ocular stack tessellation and refraction measurement before design."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import numpy as np
import trimesh
ROOT=Path(__file__).resolve().parents[2]
Z=np.array([0.,.6,3.6,7.6,23.6]);C=np.array([1/8,1/7,1/10,-1/6,0.])


def sag(r2,c):return c*r2/(1+np.sqrt(1-c*c*r2))


def measure(radial,angular):
    radii=np.arange(1,radial+1)*3/radial;angles=np.arange(angular)*2*np.pi/angular
    xy=np.vstack((np.zeros((1,2)),np.stack((radii[:,None]*np.cos(angles),radii[:,None]*np.sin(angles)),axis=2).reshape(-1,2)))
    triangles=[]
    for j in range(angular):triangles.append([0,1+j,1+(j+1)%angular])
    for ring in range(radial-1):
        a=1+ring*angular;b=a+angular
        for j in range(angular):
            q=(j+1)%angular;triangles.extend(([a+j,b+j,b+q],[a+j,b+q,a+q]))
    triangles=np.array(triangles,np.int32);vertices=[];faces=[];front=[];back=[];n=len(xy);errors=[]
    for i,(z,c) in enumerate(zip(Z,C)):
        v=np.column_stack((xy,z+sag(np.sum(xy*xy,axis=1),c)));vertices.append(v)
        fv=v[triangles];cent=fv.mean(1);delta=cent[:,2]-(z+sag(np.sum(cent[:,:2]**2,axis=1),c))
        normal=np.cross(fv[:,1]-fv[:,0],fv[:,2]-fv[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
        slope=c*cent[:,:2]/np.sqrt(1-c*c*np.sum(cent[:,:2]**2,axis=1))[:,None]
        analytic=np.column_stack((-slope,np.ones(len(slope))));analytic/=np.linalg.norm(analytic,axis=1)[:,None]
        angle=np.arccos(np.clip(np.sum(normal*analytic,axis=1),-1,1));errors.append([float(np.abs(delta).max()),float(angle.max())])
        faces.extend((triangles+i*n).tolist());front.extend([i+1 if i<4 else 0]*len(triangles));back.extend([i]*len(triangles))
    outer=1+(radial-1)*angular
    for region in range(1,5):
        for j in range(angular):
            k=(j+1)%angular;a=(region-1)*n+outer+j;b=(region-1)*n+outer+k;c=region*n+outer+j;d=region*n+outer+k
            faces.extend(([a,b,d],[a,d,c]));front.extend([0,0]);back.extend([region,region])
    vertices=np.concatenate(vertices);faces=np.array(faces,np.int32);front=np.array(front);back=np.array(back)
    volumes=[];closed=True
    for region in range(1,5):
        f=np.concatenate((faces[back==region],faces[front==region][:,::-1]));mesh=trimesh.Trimesh(vertices,f,process=False)
        volumes.append(float(mesh.volume));closed &= bool(mesh.is_watertight and mesh.is_winding_consistent and mesh.volume>0)
    return dict(radial=radial,angular=angular,vertices=len(vertices),triangles=len(faces),volumes=volumes,closed=closed,max_sag_error=max(x[0] for x in errors),max_normal_angle=max(x[1] for x in errors),mesh_hashes=[hashlib.sha256(x.tobytes()).hexdigest() for x in (vertices,faces,front,back)])


def refraction(reference):
    rows=[]
    for n1,n2 in ((1.38,1.34),(1.34,1.42),(1.42,1.34),(1.34,1.)):
        a=float(np.float32(n1));b=float(np.float32(n2))
        for angle in (0.,10.,30.,50.):
            theta=np.deg2rad(angle);d=np.array([[np.sin(theta),0,np.cos(theta)]],np.float32)
            sysd=dict(z=np.array([0.]),c=np.zeros(1),k=np.zeros(1),poly=np.zeros((1,6)),n=np.array([a,b]),ca=np.array([100.]),istop=-1)
            _,expected,alive=reference.trace(sysd,np.array([[0.,0.,-1.]]),d.astype(np.float64))
            cosine=d[0,2];eta=np.float32(a)/np.float32(b);s2=eta*eta*(np.float32(1)-cosine*cosine)
            transmitted=bool(s2<1)
            if transmitted:
                ct=np.sqrt(np.float32(1)-s2);out=eta*d[0]+np.array([0,0,ct-eta*cosine],np.float32);out/=np.linalg.norm(out)
                error=float(np.max(np.abs(out.astype(np.float64)-expected[0])))
                rs=(a*float(cosine)-b*float(ct))/(a*float(cosine)+b*float(ct));rp=(b*float(cosine)-a*float(ct))/(b*float(cosine)+a*float(ct));fresnel=.5*(rs*rs+rp*rp)
            else:error=None;fresnel=1.
            rows.append(dict(n1=a,n2=b,angle=angle,transmitted=transmitted,reference_alive=bool(alive[0]),direction_error=error,fresnel=fresnel))
    return rows


def main():
    source=Path(os.environ['LENS_REFERENCE_SOURCE']);spec=importlib.util.spec_from_file_location('lens_reference',source);reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)
    pairs=[]
    for _ in range(2):pairs.append(dict(meshes=[measure(r,a) for r,a in ((8,32),(16,64),(32,128))],refraction=refraction(reference)))
    gates=dict(full_repeat=pairs[0]==pairs[1],finite=all(np.isfinite(r['volumes']).all() and np.isfinite(r['max_sag_error']) and np.isfinite(r['max_normal_angle']) for r in pairs[0]['meshes']),closed_regions=all(r['closed'] for r in pairs[0]['meshes']))
    report=dict(**pairs[0],gates=gates,reference_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),scope='Synthetic geometry and float32 refraction mechanism only; no photon transport or anatomical acceptance.')
    (ROOT/'reports/ocular_geometry_mechanism.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
