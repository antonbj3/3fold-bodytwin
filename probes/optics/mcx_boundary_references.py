"""Measure refractive cube and slab/diffusion anchors before transport extension."""
import importlib.metadata
import hashlib
import json
from pathlib import Path
import numpy as np
import pmcx
ROOT=Path(__file__).resolve().parents[2]


def slab_profile(field):
    xy=np.arange(100,dtype=float)+.5-50.
    xx,yy=np.meshgrid(xy,xy,indexing='ij');mask=xx*xx+yy*yy<=100.
    return np.array([field[:,:,z:z+2][mask].mean() for z in range(5,21,2)])


def diffusion_profile():
    # Semi-infinite Green function, with source/mirror depths of the MCX anchor.
    mua=.01;musp=1.;diff=1/(3*(mua+musp));z0=1/(mua+musp);zb=2*diff
    x=np.arange(100,dtype=float)+.5-50.
    xx,yy,zz=np.meshgrid(x,x,np.arange(100,dtype=float)+.5,indexing='ij')
    real=np.sqrt(xx*xx+yy*yy+(zz-z0)**2);mirror=np.sqrt(xx*xx+yy*yy+(zz+z0+2*zb)**2)
    attenuation=np.sqrt(3*mua*musp)
    phi=(np.exp(-attenuation*real)/real-np.exp(-attenuation*mirror)/mirror)/(4*np.pi*diff)
    return slab_profile(phi)


def main():
    rows=[];saved={};analytic=diffusion_profile()
    for case in ('cube60b','slab'):
        fields=[];case_rows=[]
        for leg in range(2):
            cfg=pmcx.mcxcreate('cube60b' if case=='cube60b' else 'cube60')
            if case=='slab':
                cfg.update(vol=np.ones((100,100,100),dtype=np.uint8),srcpos=[50,50,0],
                           prop=[[0,0,1,1],[.01,10,.9,1]])
                cfg.pop("detpos",None)
            out=pmcx.run(cfg);fluence=np.asarray(out['flux'])[...,0]*cfg['tstep'];fields.append(fluence.copy())
            row={'case':case,'leg':leg+1,'absorbed_fraction':float(np.sum(fluence,dtype=float)*cfg['prop'][1][0]),
                 'finite_nonnegative':bool(np.isfinite(fluence).all() and fluence.min()>=0),
                 'field_sha256':hashlib.sha256(fluence.tobytes()).hexdigest()}
            if case=='slab':
                profile=slab_profile(fluence);row.update(profile=profile.tolist(),diffusion_relative_errors=(np.abs(profile-analytic)/analytic).tolist())
            saved[f'{case}_{leg}']=fluence;case_rows.append(row)
        rows.append({'case':case,'legs':case_rows,'bit_identical':fields[0].tobytes()==fields[1].tobytes(),
                     'max_repeat_difference':float(np.max(np.abs(fields[0]-fields[1])))})
    gates={'finite_nonnegative':all(x['finite_nonnegative'] for c in rows for x in c['legs']),
           'slab_diffusion_band':all(max(x['diffusion_relative_errors'])<=.05 for x in rows[1]['legs']),
           'external_exact_repeat':all(c['bit_identical'] for c in rows)}
    report={'cases':rows,'gates':gates,'diffusion_profile':analytic.tolist(),'pmcx_version':importlib.metadata.version('pmcx'),
            'slab_depth_intervals_mm':[[z,z+2] for z in range(5,21,2)],'slab_radius_mm':10.,
            'scope':'fresh external boundary and slab anchors;5% engineering diffusion band, no clinical calibration'}
    (ROOT/'reports/mcx_boundary_references.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/mcx_boundary_reference_fields.npz',**saved)
    print(json.dumps(report),flush=True);return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
