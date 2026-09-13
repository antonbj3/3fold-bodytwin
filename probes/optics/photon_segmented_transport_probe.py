"""Full photon transport parity and synchronized API timing on a fixed fixture."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import warp as wp
from bodytwin.geometry.optics import tissue_photon_paths_v1 as baseline
from bodytwin.geometry.optics import tissue_photon_segmented_v1 as candidate


def run(fixture, photons=128, seed=20260912, prepared=False, launch64=False, window2=False, cell_certified=False):
    if type(photons) is not int or photons<1 or type(seed) is not int:
        raise ValueError("Explicit positive photon count and integer seed required")
    data=np.load(fixture)
    before={k:hashlib.sha256(data[k].tobytes()).hexdigest() for k in data.files}
    v,f,front,back=(data[k] for k in ('vertices','faces','front','back'))
    prop=np.array([[0,0,1,1],[.019,7.8182,.89,1.37],[.019,7.8182,.89,1.37],[.0004,.009,.89,1.37],[.02,9,.89,1.37],[.08,40.9,.84,1.37],[0,0,1,1]])
    kwargs=dict(source=[75,67.38,167.5],direction=[.1636,.4569,-.8743],initial_region=1,origin=[0,0,0],shape=tuple(data['shape']),reflect=True,photons=photons,seed=seed)
    reference_module=baseline
    candidate_module=candidate
    if prepared:
        from bodytwin.geometry.optics import tissue_photon_prepared_v1 as candidate_module
        reference_module=candidate
    if launch64:
        from bodytwin.geometry.optics import tissue_photon_prepared_v1 as reference_module
        from bodytwin.geometry.optics import tissue_photon_launch64_v1 as candidate_module
    if window2:
        from bodytwin.geometry.optics import tissue_photon_launch64_v1 as reference_module
        from bodytwin.geometry.optics import tissue_photon_window2_v1 as candidate_module
    if cell_certified:
        from bodytwin.geometry.optics import tissue_photon_window2_v1 as reference_module
        from bodytwin.geometry.optics import tissue_photon_cell_v1 as candidate_module
    construction_ms={};captures={};rows=[];preparation_ms=None;cold_candidate_ms=None
    for name,module in (('baseline',reference_module),('candidate',candidate_module)):
        if launch64 or window2 or cell_certified:
            wp.synchronize();construct_start=time.perf_counter()
            scene=module.PreparedPhotonScene(v,f,prop,front,back)
            wp.synchronize();construction_ms[name]=(time.perf_counter()-construct_start)*1000
            execute=lambda:scene.simulate(**kwargs)
            execute();wp.synchronize()
        elif prepared and name=='candidate':
            preparation_ms=[];cold_candidate_ms=[]
            for cold_leg in range(2):
                wp.synchronize();cold_start=time.perf_counter()
                scene=module.PreparedPhotonScene(v,f,prop,front,back);wp.synchronize()
                preparation_ms.append((time.perf_counter()-cold_start)*1000)
                execute=lambda:scene.simulate(**kwargs)
                execute();wp.synchronize();cold_candidate_ms.append((time.perf_counter()-cold_start)*1000)
        else:
            execute=lambda:module.simulate(v,f,prop,front,back,**kwargs)
            execute();wp.synchronize()
        for leg in range(2):
            wp.synchronize();start=time.perf_counter()
            output=execute()
            wp.synchronize();elapsed=(time.perf_counter()-start)*1000
            captures[name,leg]=output
            rows.append(dict(name=name,leg=leg,milliseconds=elapsed,hashes={k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in output.items()},
                             events=int(output['counters'][0]),caps=int(output['counters'][2]),leaks=int(output['terminal'][:,3].sum()),
                             energy=int(output['terminal'][:,:3].sum())))
    reference=captures['baseline',0]
    gates=dict(full_transport_parity=all(reference[k].tobytes()==out[k].tobytes() for out in captures.values() for k in reference),
               complete_repeat=all(captures[name,0][k].tobytes()==captures[name,1][k].tobytes() for name in ('baseline','candidate') for k in reference),
               energy_and_leaks=all(row['energy']==photons*baseline.PACKET and row['leaks']==0 and row['caps']==0 for row in rows),
               faster_both=max(row['milliseconds'] for row in rows if row['name']=='candidate')<min(row['milliseconds'] for row in rows if row['name']=='baseline'),
               input_immutable=before=={k:hashlib.sha256(data[k].tobytes()).hexdigest() for k in data.files})
    return dict(photons=photons,seed=seed,prepared=prepared,launch64=launch64,window2=window2,cell_certified=cell_certified,construction_ms=construction_ms,preparation_ms=preparation_ms,cold_candidate_ms=cold_candidate_ms,rows=rows,gates=gates,overall_pass=all(gates.values()),source_hashes={name:hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() for name,module in (('baseline',reference_module),('candidate',candidate_module))}),captures


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fixture',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--photons',type=int,default=128);p.add_argument('--seed',type=int,default=20260912);p.add_argument('--prepared',action='store_true');p.add_argument('--launch64',action='store_true');p.add_argument('--window2',action='store_true');p.add_argument('--cell-certified',action='store_true');a=p.parse_args()
    report,_=run(a.fixture,a.photons,a.seed,a.prepared,a.launch64,a.window2,a.cell_certified);a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));raise SystemExit(0 if report['overall_pass'] else 2)
