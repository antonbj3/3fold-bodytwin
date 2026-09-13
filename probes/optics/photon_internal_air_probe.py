"""Frozen propagation through a measured synthetic internal transparent region."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.label_interfaces_v1 import label_interfaces


def main():
    from bodytwin.geometry.optics.tissue_photon_paths_v1 import simulate,PACKET
    labels=np.ones((3,3,3),np.uint8);labels[1,1,1]=2
    mesh=label_interfaces(labels);rows=[]
    for mu in (0.,.2):
        properties=np.array([[0.,0.,1.,1.],[mu,0.,1.,1.],[0.,0.,1.,1.]])
        outputs=[simulate(mesh['vertices'],mesh['faces'],properties,mesh['front'],mesh['back'],
                          source=[1.25,1.375,.0001],direction=[0,0,1],initial_region=1,
                          origin=[0,0,0],shape=[3,3,3]) for _ in range(2)]
        out=outputs[0];terminal=out['terminal'];sums=terminal[:,:3].sum(axis=0,dtype=np.int64)
        prediction=np.floor(PACKET*np.exp(-mu*out['paths'][:,0])+.5).astype(np.int64)
        rows.append(dict(mu=mu,full_repeat=all(out[k].tobytes()==outputs[1][k].tobytes() for k in out),
                         hashes={k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in out.items()},
                         absorbed=int(sums[0]),escaped=int(sums[1]),residual=int(sums[2]),
                         leaks=int(terminal[:,3].sum()),caps=int(out['counters'][2]),
                         exit_z_max_error=float(np.max(np.abs(out['exits'][:,2]-3))),
                         region_length_max_errors=np.max(np.abs(out['paths']-[1.9999,1.]),axis=0).tolist(),
                         all_regions_visited=bool(np.all(out['paths']>0)),
                         packet_attenuation_max_error=int(np.max(np.abs(prediction-terminal[:,1])))))
    gates=dict(full_repeat=all(r['full_repeat'] for r in rows),
               exact_energy=all(r['absorbed']+r['escaped']==1000000*PACKET for r in rows),
               zero_failed_packets=all(r['residual']==r['leaks']==r['caps']==0 for r in rows),
               outer_exit=all(r['exit_z_max_error']==0 for r in rows),
               internal_region_visited=all(r['all_regions_visited'] for r in rows),
               path_bound=all(max(r['region_length_max_errors'])<=1e-4 for r in rows),
               attenuation_bound=all(r['packet_attenuation_max_error']<=1 for r in rows),
               transparent_no_absorption=rows[0]['absorbed']==0)
    report=dict(rows=rows,gates=gates,scope='Synthetic ballistic internal-air control, not anatomical transport validation.')
    (ROOT/'reports/photon_internal_air.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
