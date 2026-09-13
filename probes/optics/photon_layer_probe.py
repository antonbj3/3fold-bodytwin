"""Certify frozen v6 on the measured two-region mesh and external fluence."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import warp as wp
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.layered_box_v1 import layered_box
from bodytwin.geometry.optics.tissue_photon_mc_v6 import simulate,PACKET
PROPERTIES=np.array([[0.,0.,1.,1.],[.005,1.,.01,1.37],[.01,2.,.8,1.4]])


def coarse(a):return a.reshape(12,5,12,5,12,5).sum(axis=(1,3,5))


def main():
    wp.init();mesh=layered_box();references=np.load(ROOT/'reports/mcx_layer_reference_fields_l4.npz')
    labels=np.ones((60,60,60),dtype=int);labels[:,:,10:]=2;mua=PROPERTIES[labels,0]
    rows=[];outputs=[]
    for _ in range(2):
        out=simulate(mesh['vertices'],mesh['faces'],PROPERTIES,mesh['front'],mesh['back'],source=[29,29,.0001],
            direction=[0,0,1],initial_region=1,origin=[0,0,0],shape=[60]*3,reflect=True)
        outputs.append(out);sums=out['terminal'][:,:3].sum(axis=0,dtype=np.int64)
        absorbed=out['absorption'].astype(float)/(1000000*PACKET);fluence=absorbed/mua
        refs=[references[f'field_{leg}'].astype(float) for leg in (0,1)]
        layer_abs=[float(absorbed[labels==i].sum()) for i in (1,2)]
        errors=[]
        for a in refs:
            refabs=a*mua
            errors.append({'total_absorption':float(abs(absorbed.sum()-refabs.sum())/refabs.sum()),
                'layer_absorption':[float(abs(layer_abs[i-1]-refabs[labels==i].sum())/refabs[labels==i].sum()) for i in (1,2)],
                'binned_fluence':float(np.abs(coarse(fluence)-coarse(a)).sum()/a.sum())})
        row={'hashes':{k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in out.items()},
            'absorbed_integer':int(sums[0]),'escaped_integer':int(sums[1]),'residual_integer':int(sums[2]),
            'leaks':int(out['terminal'][:,3].sum()),'event_caps':int(out['counters'][2]),'events':int(out['counters'][0]),'boundary_hits':int(out['counters'][1]),
            'window_absorbed_fraction':float(absorbed.sum()),'layer_absorbed_fractions':layer_abs,'reference_errors':errors,
            'finite_nonnegative':bool(np.isfinite(fluence).all() and fluence.min()>=0)}
        rows.append(row)
    gates={'full_arrays_repeat_exact':all(outputs[0][k].tobytes()==outputs[1][k].tobytes() for k in outputs[0]),
        'exact_energy':all(r['absorbed_integer']+r['escaped_integer']==1000000*PACKET for r in rows),
        'no_failed_packets':all(r['residual_integer']==r['leaks']==r['event_caps']==0 for r in rows),
        'finite_nonnegative':all(r['finite_nonnegative'] for r in rows),
        'total_absorption':all(e['total_absorption']<=.02 for r in rows for e in r['reference_errors']),
        'per_layer_absorption':all(max(e['layer_absorption'])<=.02 for r in rows for e in r['reference_errors']),
        'binned_fluence':all(e['binned_fluence']<=.05 for r in rows for e in r['reference_errors'])}
    report={'legs':rows,'gates':gates,'properties':PROPERTIES.tolist(),
        'backend_sha256':hashlib.sha256((ROOT/'src/bodytwin/geometry/optics/tissue_photon_mc_v6.py').read_bytes()).hexdigest(),
        'scope':'fixed synthetic two-region transport; detector and physiology integration still separate'}
    (ROOT/'reports/photon_layer.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
