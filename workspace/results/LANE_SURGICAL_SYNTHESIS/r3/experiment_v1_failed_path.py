import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import argparse,copy,hashlib,json,time,resource
from pathlib import Path
import numpy as np
from surgical_chain.chain import write,run,resume
from surgical_chain.scenario import scenario_config
from surgical_chain.shared import late_response
from surgical_chain.vendor import response_r1 as r
from diagnose_r2 import score,load_reference
ROOT=Path(__file__).resolve().parent

def all_r3_config():
    cfg=scenario_config()
    p=json.loads((ROOT/'sources/LANE_SURGICAL_RESPONSE/r1/native_chain/native_geometry_response.json').read_text())['parameters']
    for key,pk in {'gap':'gap_m','biological_width':'width_m','perfusion_width':'width_m','face_conductance':'face_conductance','k_deposit':'k_deposit_day','k_mature_legacy':'k_mature_day'}.items():
        cfg['parameters'][key]['value']=p[pk]
        cfg['parameters'][key]['source']='Frozen RESPONSE R1 native config inherited by R3; early strength fitted k_deposit/k_mature, not measured collagen'
    return cfg

def main(out):
    assert hashlib.sha256((ROOT/'PREREG_R3.json').read_bytes()).hexdigest()==(ROOT/'PREREG_R3.sha256').read_text().strip()
    out.mkdir(exist_ok=False,parents=True)
    st=time.perf_counter();cpu=time.process_time();rows=[]
    for name,cfg in [('nominal',scenario_config()),('all_R3',all_r3_config()),('half_dt',all_r3_config()),('half_dx',all_r3_config())]:
        if name=='half_dt':cfg['parameters']['late_dt']['value']=.05;cfg['parameters']['early_dt']['value']=.05
        if name=='half_dx':cfg['parameters']['grid_dx']['value']=2.5e-6
        d=run(cfg,mode='scenario',out=out/name,inventory='shared')
        with np.load(out/name/'fields.npz') as f:days=f['days']; pred=f['strength_pct']
        with np.load(out/name/'shared_reference.npz') as f:
            chem=f['chemical_strength_pct']; C=f['C']; u=f['UIM']
        row=dict(case=name,score=score(days,pred),chemical_score=score(days,chem),shared_reference=d['shared_reference'],costs=d['costs'],strength90_pct=float(pred[-1]),C90=float(C[-1]),UIM90=u[-1].tolist())
        rows.append(row);print(name,row['score']['RMSE_pp'],row['chemical_score']['RMSE_pp'],d['shared_reference']['max_pointwise_trace_error'],flush=True)
    cfg=all_r3_config();params={k:v['value'] for k,v in cfg['parameters'].items()}
    meta=json.loads((out/'all_R3/checkpoint21.json').read_text());p=r.Params(**meta['parameters']);g=r.Grid(p,dx=params['grid_dx'])
    with np.load(out/'all_R3/fields.npz') as z:early=dict(state=z['early_state'].copy(),hazard=z['hazard'].copy())
    for name in ['noformation','noremove','pathway_block_k1','pathway_block_k2','cohort_control']:
        pp=p;pa=copy.deepcopy(params)
        if name=='noformation':pp=r.replace(p,k_deposit_day=0)
        if name=='pathway_block_k1':pa['k_U_I']=0
        if name=='pathway_block_k2':pa['k_I_M']=0
        s=time.perf_counter();l=late_response(g,pp,early,pa,90.,.1,removal=name!='noremove',compare_cohorts=True)
        raw=out/(name+'.npz');np.savez_compressed(raw,days=l['days'],C=l['C'],UIM=l['UIM'],strength_pct=l['strength_pct'],chemical_marks=l['chemical_marks_series'],chemical_strength_pct=l['chemical_strength_pct'],ledger=l['ledger'],Q=l['Q_final'])
        row=dict(case=name,score=score(l['days'],l['strength_pct']),shared_reference=l['shared_reference'],wall_s=time.perf_counter()-s,strength90_pct=float(l['strength_pct'][-1]),C90=float(l['C'][-1]))
        rows.append(row);print(name,row['score']['RMSE_pp'],row['strength90_pct'],flush=True)
    s=time.perf_counter();l=resume(out/'all_R3/checkpoint21.npz',out/'restart')
    with np.load(out/'all_R3/fields.npz') as z:
        mask=z['days']>21+1e-8
        restart=dict(max_strength_error_pp=float(abs(l['strength_pct']-z['strength_pct'][mask]).max()),max_tensor_error=float(abs(l['Q_final']-z['Q_final']).max()),wall_s=time.perf_counter()-s,prefix_replayed=False,global_suffix_recomputed=True)
    with np.load(out/'all_R3/fields.npz') as z:base=z['strength_pct'];bd=z['days']
    refinement={}
    for name in ['half_dt','half_dx']:
        with np.load(out/name/'fields.npz') as z:refinement[name]=float(abs(base-np.interp(bd,z['days'],z['strength_pct'])).max())
    ref=load_reference()
    with np.load(out/'all_R3/shared_reference.npz') as z:
        observed=75*z['C']*np.interp(z['days'],ref['days'],ref['HP_wound_reference_fraction'])
        observed_score=score(z['days'],observed)
    # Reference uses sparse FV quantity interpolation and external HP; preserve exact older diagnostic.
    reference_score=score(ref['days'],ref['strength_pct'])
    result=dict(review_state='PENDING_INDEPENDENT_REVIEW',scientific_admission=False,rows=rows,reference_score=reference_score,external_HP_with_new_C_diagnostic=observed_score,external_HP_fed_into_candidate=False,late_strength_fit=False,refinement_max_strength_pp=refinement,restart=restart,engineering_gate='PASS' if all(x['shared_reference']['max_pointwise_trace_error']<1e-10 and x['shared_reference']['max_collagen_ledger_error']<1e-10 for x in rows) and max(refinement.values())<1 and restart['max_strength_error_pp']<1e-9 else 'FAIL',autonomous_component_reproduction_gate='PASS' if rows[1]['score']['RMSE_pp']<=8.406367094671637 else 'FAIL',strongest_control_gate='TIE',empirical_joint_gate='UNKNOWN',total_wall_s=time.perf_counter()-st,total_CPU_s=time.process_time()-cpu,max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    write(out/'EXPERIMENT.json',result)
    print('engineering',result['engineering_gate'],'autonomous',result['autonomous_component_reproduction_gate'],'refinement',refinement,'restart',restart,flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--out',type=Path,required=True);main(a.parse_args().out)
