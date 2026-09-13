"""Check spec-generated chains against immutable propagation and observer records."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import os
import subprocess
import tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.chains.nirs_spec_v1 import Chain,generate
from bodytwin.chains.nirs_chain_v1 import leg,propagate


def original_gates(legs):
    return dict(full_chain_repeat_exact=legs[0]==legs[1],finite_physical=all(r['finite_physical'] for r in legs),
        zero_band=all(r['zero_band_exact'] for r in legs),zero_blood=all(r['zero_blood_ratio']==1 for r in legs),
        oxygen_direction=all(r['oxygen_direction'] for r in legs),
        sample_composition=all(r['joint_span']<=r['sum_isolated_spans'] for r in legs),
        grid_composition=all(r['grid_span']<=r['grid_isolated_span_sum'] for r in legs),
        samples_inside_measured_grid=all(r['random_within_grid'] for r in legs))


def main():
    rows=[];arrays={}
    for name in ('box','curved','narrow'):
        raw=(ROOT/f'examples/chains/nirs_{name}.json').read_bytes();spec=json.loads(raw);chain=Chain(spec)
        filename='photon_detector_paths_l4.npz' if spec['fixture']=='box' else 'photon_curved_detector_paths_l4.npz'
        paths=np.load(ROOT/'reports'/filename)['paths_0'];legs=[chain.observe(paths),chain.observe(paths)]
        inputs=chain.nominal+chain.bands*np.random.default_rng(spec['seed']).uniform(-1.,1.,(spec['draws'],3))
        values=chain.propagate(paths,inputs);arrays[name]=values
        row=dict(name=name,spec_sha256=hashlib.sha256(raw).hexdigest(),legs=legs,gates=original_gates(legs),
                 full_records_reference_exact=values.tobytes()==propagate(paths,inputs).tobytes())
        if name!='narrow':row['frozen_observer_bytes_exact']=json.dumps(legs[0],sort_keys=True).encode()==json.dumps(leg(paths),sort_keys=True).encode()
        generated=generate(spec)
        row['generated_bytes_exact']=generated.encode()==generate(json.loads(raw)).encode()
        with tempfile.TemporaryDirectory() as folder:
            script=Path(folder)/'chain.py';script.write_text(generated)
            env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'src')
            children=[subprocess.run([sys.executable,str(script),str(ROOT/'reports'/filename)],env=env,capture_output=True,text=True,check=True,timeout=120).stdout for _ in range(2)]
        row['generated_execution_exact']=children[0]==children[1] and json.loads(children[0])==legs[0]
        rows.append(row)
    bad=[];base=json.loads((ROOT/'examples/chains/nirs_box.json').read_text())
    x=copy.deepcopy(base);x['nodes'][0]['module']='unknown';bad.append(x)
    x=copy.deepcopy(base);x['nodes'][0]['inputs']={};bad.append(x)
    x=copy.deepcopy(base);x['nodes'][0]['inputs']['po2']='optical';bad.append(x)
    x=copy.deepcopy(base);x['bands']=[-1,0,0];bad.append(x)
    x=copy.deepcopy(base);x['nodes'][0]['inputs']['po2']='input.paths';bad.append(x)
    rejected=0
    for spec in bad:
        try:Chain(spec)
        except ValueError:rejected+=1
    gates=dict(frozen_two_exact=all(r.get('frozen_observer_bytes_exact',True) and r['full_records_reference_exact'] for r in rows),
               all_repeat=all(r['legs'][0]==r['legs'][1] for r in rows),original_gates=all(all(r['gates'].values()) for r in rows),
               third_narrower=rows[2]['legs'][0]['joint_span']<rows[1]['legs'][0]['joint_span'],invalid_rejected=rejected==5,generated=all(r['generated_bytes_exact'] and r['generated_execution_exact'] for r in rows))
    report=dict(rows=rows,gates=gates,rejected_controls=rejected,scope='Fixed registry NIRS data graphs; conditional synthetic paths; no clinical calibration or arbitrary program generation.')
    (ROOT/'reports/nirs_spec_probe.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(ROOT/'reports/nirs_spec_records.npz',**arrays)
    print(json.dumps(dict(gates=gates,rows=[dict(name=r['name'],gates=r['gates'],joint_span=r['legs'][0]['joint_span'],frozen_exact=r.get('frozen_observer_bytes_exact')) for r in rows])))
    return 0 if all(gates.values()) else 2



def reduced_compare():
    import time
    from bodytwin.chains import nirs_chain_v1 as baseline
    from bodytwin.chains import nirs_reduced_chain_v1 as candidate
    rows=[];arrays={};all_exact=True;all_gates=True;unchanged=True
    for filename in ('photon_detector_paths_l4.npz','photon_curved_detector_paths_l4.npz'):
        paths=np.load(ROOT/'reports'/filename)['paths_0'];before=paths.tobytes()
        inputs=baseline.NOMINAL+baseline.BANDS*np.random.default_rng(baseline.SEED).uniform(-1.,1.,(baseline.DRAWS,3))
        baseline.leg(paths);candidate.leg(paths);results={}
        for repeat,order in enumerate((('baseline','candidate'),('candidate','baseline'))):
            for name in order:
                module=baseline if name=='baseline' else candidate
                start=time.perf_counter();record=module.leg(paths);seconds=time.perf_counter()-start
                values=module.propagate(paths,inputs);arrays[f'{filename}_{name}_{repeat}']=values
                results[name,repeat]=record
                rows.append(dict(fixture=filename,name=name,repeat=repeat,seconds=seconds,record=record,propagation_sha256=hashlib.sha256(values.tobytes()).hexdigest()))
        all_exact &= all(record==results['baseline',0] for record in results.values())
        all_exact &= all(arrays[f'{filename}_{name}_{rep}'].tobytes()==arrays[f'{filename}_baseline_0'].tobytes() for name in ('baseline','candidate') for rep in (0,1))
        all_gates &= all(original_gates([results['candidate',0],results['candidate',1]]).values())
        unchanged &= paths.tobytes()==before
    # Positive attenuation weights: zeros, subnormals, large dynamic range, and odd/even row counts.
    controls=[]
    for count in (1,2,3,17,4095,4096,4119):
        weights=np.random.default_rng(20260913).uniform(size=(count,2))
        controls.append(weights)
        controls.append(np.resize(np.array([[0.,1.],[np.nextafter(0.,1.),1e-300],[1.,1e-16],[1e-16,1.]]),(count,2)))
    reduction_exact=all(a.sum(axis=0).tobytes()==np.add.accumulate(a,axis=0)[-1].tobytes() for a in controls)
    speed=all(max(r['seconds'] for r in rows if r['fixture']==f and r['name']=='candidate')<min(r['seconds'] for r in rows if r['fixture']==f and r['name']=='baseline') for f in ('photon_detector_paths_l4.npz','photon_curved_detector_paths_l4.npz'))
    gates=dict(complete_exact=bool(all_exact),original_eight=bool(all_gates),inputs_unchanged=bool(unchanged),reduction_controls=reduction_exact,faster_both=bool(speed))
    report=dict(rows=rows,gates=gates,reduction_controls=len(controls),source_sha256={name:hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest() for name,m in (('baseline',baseline),('candidate',candidate))},scope='Existing synthetic NIRS chains only; no clinical or anatomy promotion.')
    (ROOT/'reports/nirs_reduction.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/nirs_reduction_arrays.npz',**arrays)
    print(json.dumps(report));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(reduced_compare() if '--reduced' in sys.argv[1:] else main())
