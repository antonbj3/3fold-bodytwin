"""Certify unperturbed paths and a fixed detector against external path statistics."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import warp as wp
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.layered_box_v1 import layered_box
from bodytwin.geometry.optics.tissue_photon_paths_v1 import simulate,PACKET
PROPERTIES=np.array([[0.,0.,1.,1.],[.005,1.,.01,1.37],[.01,2.,.8,1.4]])
N=1000000


def statistics(weights,paths):
    total=weights.sum();fraction=total/N;mean=np.sum(weights[:,None]*paths,axis=0)/total
    return {'count':len(weights),'fraction':float(fraction),'fraction_variance':float((np.sum(weights**2)/N-fraction**2)/(N-1)),
        'layer_path_mean':mean.tolist(),'layer_path_variance':(np.sum(weights[:,None]**2*(paths-mean)**2,axis=0)/total**2).tolist()}


def main():
    wp.init();m=layered_box();base=json.loads((ROOT/'reports/photon_layer_l4.json').read_text())['legs'][0]['hashes']
    ref=np.load(ROOT/'reports/mcx_layer_reference_fields_l4.npz');refs=[]
    for leg in (0,1):
        p=ref[f'paths_{leg}'];w=np.exp(-(p*PROPERTIES[1:,0]).sum(axis=1));refs.append(statistics(w,p))
    rows=[];saved={}
    for leg in range(2):
        out=simulate(m['vertices'],m['faces'],PROPERTIES,m['front'],m['back'],source=[29,29,.0001],direction=[0,0,1],
            initial_region=1,origin=[0,0,0],shape=[60]*3,reflect=True)
        hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in out.items()}
        canonical={k:hashes[k] for k in base};escaped=out['terminal'][:,1]>0
        reconstructed=np.floor(PACKET*np.exp(-(out['paths']*PROPERTIES.astype(np.float32).astype(float)[1:,0]).sum(axis=1))+.5).astype(np.int64)
        energy_error=int(np.max(np.abs(reconstructed[escaped]-out['terminal'][escaped,1])))
        distance2=np.sum((out['exits'][:,:3]-[39,29,0])**2,axis=1)
        detected=escaped&(out['exits'][:,6]<=5.)&(distance2<=4.)
        p=out['paths'][detected];w=out['terminal'][detected,1].astype(float)/PACKET
        stat=statistics(w,p);comparisons=[]
        for reference in refs:
            delta=abs(stat['fraction']-reference['fraction']);bound=5*np.sqrt(stat['fraction_variance']+reference['fraction_variance'])
            path_delta=np.abs(np.array(stat['layer_path_mean'])-reference['layer_path_mean'])
            path_bound=5*np.sqrt(np.array(stat['layer_path_variance'])+reference['layer_path_variance'])
            comparisons.append({'fraction_difference':float(delta),'fraction_bound':float(bound),'fraction_pass':bool(delta<=bound),
                'path_differences':path_delta.tolist(),'path_bounds':path_bound.tolist(),'paths_pass':bool(np.all(path_delta<=path_bound))})
        row={'hashes':hashes,'canonical_identical':canonical==base,'reconstructed_energy_max_integer_error':energy_error,
            'detector':stat,'exit_z_max_abs':float(np.max(np.abs(out['exits'][detected,2]))),'comparisons':comparisons}
        rows.append(row);saved[f'paths_{leg}']=p;saved[f'weights_{leg}']=w;saved[f'exits_{leg}']=out['exits'][detected]
    gates={'canonical_transport_exact':all(r['canonical_identical'] for r in rows),'all_full_arrays_repeat_exact':rows[0]['hashes']==rows[1]['hashes'],
        'nonempty_detector':all(r['detector']['count']>0 for r in rows),
        'path_energy_reconstruction':all(r['reconstructed_energy_max_integer_error']<=1 for r in rows),
        'detector_on_exit_face':all(r['exit_z_max_abs']<=1e-8 for r in rows),
        'external_fraction_five_standard_errors':all(c['fraction_pass'] for r in rows for c in r['comparisons']),
        'external_paths_five_standard_errors':all(c['paths_pass'] for r in rows for c in r['comparisons'])}
    report={'legs':rows,'references':refs,'gates':gates,'scope':'fixed layered detector and nonperturbing metadata, not NIRS/clinical calibration'}
    (ROOT/'reports/photon_detector.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/photon_detector_paths.npz',**saved)
    print(json.dumps(report));return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
