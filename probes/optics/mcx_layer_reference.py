"""Measure a fixed two-layer external field and detector-path reference."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import numpy as np
import pmcx
ROOT=Path(__file__).resolve().parents[2]
PROPERTIES=np.array([[0.,0.,1.,1.],[.005,1.,.01,1.37],[.01,2.,.8,1.4]])


def main():
    fields=[];rows=[];saved={};detectors=[[39.,29.,0.,2.]]
    for leg in range(2):
        cfg=pmcx.mcxcreate('cube60b');vol=np.ones((60,60,60),dtype=np.uint8);vol[:,:,10:]=2
        cfg.update(vol=vol,prop=PROPERTIES,detpos=detectors,issavedet=1,issaveexit=1,savedetflag='dpxv',maxdetphoton=1000000)
        out=pmcx.run(cfg);field=np.asarray(out['flux'])[...,0]*cfg['tstep'];fields.append(field.copy())
        detp=pmcx.detphoton(out['detp'],2,'dpxv');paths=np.asarray(detp['ppath'],dtype=float)
        weights=np.exp(-(paths*PROPERTIES[1:,0]).sum(axis=1))
        helper=np.asarray(pmcx.detweight(detp,PROPERTIES),dtype=float).reshape(-1)
        ids=np.asarray(detp['detid']).reshape(-1);positions=np.asarray(detp['p'],dtype=float)
        mua=PROPERTIES[vol,0];weighted=field.astype(float)*mua
        detector_rows=[]
        for detid in range(1,len(detectors)+1):
            keep=ids==detid;w=weights[keep];pp=paths[keep]
            detector_rows.append({'id':detid,'photons':int(keep.sum()),'detected_fraction':float(w.sum()/cfg['nphoton']),
                'weighted_layer_path_mm':(np.sum(w[:,None]*pp,axis=0)/w.sum()).tolist() if w.sum()>0 else [],
                'second_layer_visitors':int(np.count_nonzero(pp[:,1]>0))})
        row={'absorbed_fraction':float(weighted.sum()),'layer_absorbed_fractions':[float(weighted[vol==i].sum()) for i in (1,2)],
             'field_sha256':hashlib.sha256(field.tobytes()).hexdigest(),'detector_path_sha256':hashlib.sha256(paths.tobytes()).hexdigest(),
             'detectors':detector_rows,'weight_helper_max_difference':float(np.max(np.abs(weights-helper))),
             'finite_nonnegative':bool(np.isfinite(field).all() and field.min()>=0 and np.isfinite(paths).all() and paths.min()>=0),
             'detector_exit_z_max_abs':float(np.max(np.abs(positions[:,2])))}
        rows.append(row);saved[f'field_{leg}']=field;saved[f'paths_{leg}']=paths;saved[f'detector_ids_{leg}']=ids;saved[f'exit_positions_{leg}']=positions
    gates={'finite_nonnegative':all(r['finite_nonnegative'] for r in rows),
           'physical_absorbed_fraction':all(0<r['absorbed_fraction']<1 for r in rows),
           'nonempty_detector':all(r['detectors'][0]['photons']>0 for r in rows),
           'weight_formula_matches_external_helper':all(r['weight_helper_max_difference']<=1e-12 for r in rows),
           'external_full_field_repeat_exact':fields[0].tobytes()==fields[1].tobytes()}
    report={'legs':rows,'gates':gates,'max_field_repeat_difference':float(np.max(np.abs(fields[0]-fields[1]))),
        'properties':PROPERTIES.tolist(),'layer_split_mm':10.,'shape':[60,60,60],'source':[29,29,0],'detectors':detectors,
        'photons':1000000,'window_ns':5.,'seed':1648335518,'pmcx_version':importlib.metadata.version('pmcx'),
        'source':'https://github.com/fangq/mcx/blob/master/pmcx/pmcx/utils.py',
        'scope':'synthetic layered external reference, inside-tissue launch; no clinical calibration or throughput claim'}
    (ROOT/'reports/mcx_layer_reference.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/mcx_layer_reference_fields.npz',**saved)
    print(json.dumps(report),flush=True);return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
