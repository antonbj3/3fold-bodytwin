"""Synthetic wavelength-resolved retinal receiver with an external refraction oracle."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))


def launch_grid():
    x,y=np.meshgrid(np.linspace(-.9,.9,7),np.linspace(-.9,.9,7),indexing='ij');keep=x*x+y*y<=.9**2+1e-14
    p=np.column_stack((x[keep],y[keep],np.full(keep.sum(),-.5))).astype(np.float32)
    return p,np.tile(np.array([0,0,1],np.float32),(len(p),1))


def properties(wavelength):
    indices=[1.,1.,1.382,1.342,1.425,1.342] if wavelength==532 else [1.,1.,1.376,1.336,1.415,1.336]
    absorption=[0,0,.005,.0002,.01,.0002] if wavelength==532 else [0,0,.003,.0001,.006,.0001]
    return np.array([[a,0,0,n] for a,n in zip(absorption,indices)],np.float32)


def oracle(reference,p,d,props):
    sysd=dict(z=np.array([0.,.6,3.6,7.6]),c=np.array([1/8,1/7,1/10,-1/6]),k=np.zeros(4),poly=np.zeros((4,6)),n=props[1:,3].astype(np.float64),ca=np.full(4,3.),istop=-1)
    record=[];end,direction,alive=reference.trace(sysd,p.astype(np.float64),d.astype(np.float64),record=record)
    receiver=reference.to_plane(end,direction,23.6);prob=np.ones(len(p));previous=p.astype(np.float64);optical_depth=np.zeros(len(p))
    for i,hit,valid in record:
        segment=hit-previous;length=np.linalg.norm(segment,axis=1);incoming=segment/length[:,None]
        radius=np.linalg.norm(hit[:,:2],axis=1);slope=sysd['c'][i]*hit[:,:2]/np.sqrt(1-sysd['c'][i]**2*radius**2)[:,None]
        normal=np.column_stack((-slope,np.ones(len(p))));normal/=np.linalg.norm(normal,axis=1)[:,None]
        ci=np.sum(incoming*normal,axis=1);n1=sysd['n'][i];n2=sysd['n'][i+1];ct=np.sqrt(1-(n1/n2)**2*(1-ci*ci))
        rs=(n1*ci-n2*ct)/(n1*ci+n2*ct);rp=(n2*ci-n1*ct)/(n2*ci+n1*ct);prob*=1-.5*(rs*rs+rp*rp)
        optical_depth+=float(props[i+1,0])*length;previous=hit
    optical_depth+=float(props[5,0])*np.linalg.norm(receiver-end,axis=1)
    return receiver,direction,prob,np.exp(-optical_depth),alive


def main():
    from bodytwin.geometry.optics.ocular_photon_v1 import simulate,PACKET
    source=Path(os.environ['LENS_REFERENCE_SOURCE']);spec=importlib.util.spec_from_file_location('reference',source);reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)
    photons=int(os.environ.get('OCULAR_PHOTONS','10000'));p,d=launch_grid();cohort=np.arange(photons)%len(p);launched=np.bincount(cohort,minlength=len(p))
    rows=[];arrays={};edges=np.linspace(-.25,.25,65)
    for wavelength in (532,650):
        props=properties(wavelength);expected,expected_d,prob,attenuation,alive=oracle(reference,p,d,props)
        for leg in range(2):
            out=simulate(props,sources=p,directions=d,photons=photons)
            terminal=out['terminal'];sums=terminal[:,:3].sum(axis=0,dtype=np.int64);retina=out['history'][:,1]==5;direct=retina & (out['history'][:,0]==0)
            observed=np.bincount(cohort[direct],minlength=len(p));expected_count=launched*prob;bound=6*np.sqrt(launched*prob*(1-prob))+1
            position_error=float(np.max(np.abs(out['exits'][direct,:3]-expected[cohort[direct]]))) if direct.any() else float('inf')
            direction_error=float(np.max(np.abs(out['exits'][direct,3:6]-expected_d[cohort[direct]]))) if direct.any() else float('inf')
            # Integer histogram of escaped packet energy, with explicit overflow bin.
            x=np.searchsorted(edges,out['exits'][retina,0],side='right')-1;y=np.searchsorted(edges,out['exits'][retina,1],side='right')-1;inside=(x>=0)&(x<64)&(y>=0)&(y<64)
            histogram=np.zeros((64,64),np.int64);np.add.at(histogram,(x[inside],y[inside]),terminal[retina,1][inside]);overflow=int(terminal[retina,1][~inside].sum(dtype=np.int64))
            hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in out.items()};hashes['retinal_histogram']=hashlib.sha256(histogram.tobytes()).hexdigest()
            rows.append(dict(wavelength=wavelength,leg=leg,photons=photons,launches=len(p),absorbed=int(sums[0]),escaped=int(sums[1]),residual=int(sums[2]),leaks=int(terminal[:,3].sum()),caps=int(out['counters'][2]),retinal_packets=int(retina.sum()),retinal_energy=int(terminal[retina,1].sum(dtype=np.int64)),histogram_energy=int(histogram.sum(dtype=np.int64)),histogram_overflow_energy=overflow,position_error=position_error,direction_error=direction_error,max_fresnel_count_error=float(np.max(np.abs(observed-expected_count))),max_count_error_over_bound=float(np.max(np.abs(observed-expected_count)/bound)),finite_nonnegative=bool(all(np.isfinite(out[k]).all() and np.min(out[k])>=0 for k in ('absorption','terminal','paths'))),all_reference_alive=bool(alive.all()),hashes=hashes))
            for k,a in out.items():arrays[f'{wavelength}_{leg}_{k}']=a
            arrays[f'{wavelength}_{leg}_retinal_histogram']=histogram
        arrays[f'{wavelength}_reference_positions']=expected;arrays[f'{wavelength}_reference_directions']=expected_d;arrays[f'{wavelength}_direct_probability']=prob;arrays[f'{wavelength}_direct_attenuation']=attenuation
    gates=dict(exact_energy=all(r['absorbed']+r['escaped']==photons*PACKET for r in rows),zero_failed=all(r['residual']==r['leaks']==r['caps']==0 for r in rows),full_repeat=all(rows[i]['hashes']==rows[i+1]['hashes'] for i in (0,2)),finite_scoring=all(r['finite_nonnegative'] and r['histogram_energy']+r['histogram_overflow_energy']==r['retinal_energy'] for r in rows),refraction=all(r['all_reference_alive'] and r['position_error']<=3e-4 and r['direction_error']<=2e-6 for r in rows),fresnel_statistics=all(r['max_count_error_over_bound']<=1 for r in rows))
    report=dict(rows=rows,gates=gates,photons_per_arm=photons,retinal_bin_edges=edges.tolist(),reference_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),transport_sha256=hashlib.sha256((ROOT/'src/bodytwin/geometry/optics/ocular_photon_v1.py').read_bytes()).hexdigest(),scope='Synthetic spherical ocular stack, illustrative spectral indices/absorption; ideal retinal receiver. Ballistic Fresnel sampling, no anatomical/clinical calibration or scattering claim.')
    tag='ocular_transport_'+str(photons);(ROOT/'reports'/f'{tag}.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports'/f'{tag}_arrays.npz',**arrays);print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
