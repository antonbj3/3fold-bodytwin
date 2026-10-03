"""Frozen-port diagnosis; never changes the R1 chain or fits a target."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[k]='2'
import argparse, copy, hashlib, itertools, json, time
from pathlib import Path
import numpy as np
from surgical_chain.chain import run, write, late_response, UNITS
from surgical_chain.scenario import scenario_config
from surgical_chain.vendor import response_r1 as r

ROOT=Path(__file__).resolve().parent
PEER=ROOT.parent/'LANE_SURGICAL_RESPONSE'

def load_reference():
    with np.load(ROOT/'sources/LANE_SURGICAL_RESPONSE/r3/crosslink_port/native.npz') as z:
        return {k:z[k].copy() for k in z.files}

def targets():
    return json.loads((ROOT/'r1/INTEGRATED_STRENGTH_DIAGNOSTIC.json').read_text())['rows']

def score(days, prediction):
    rows=[]
    for row in targets():
        value=float(np.interp(row['day'],days,prediction))
        rows.append(dict(day=row['day'],target_pct=row['target_pct'],prediction_pct=value,
                         error_pp=value-row['target_pct'],independent_day=row['independent_day']))
    return dict(RMSE_pp=float(np.sqrt(np.mean([x['error_pp']**2 for x in rows if x['independent_day']]))),rows=rows)

def factors(z):
    abundance=z['UIM'].sum(axis=1)
    mature=np.divide(z['UIM'][:,2],abundance,out=np.zeros_like(abundance),where=abundance>0)
    return abundance,mature

def main(out):
    out.mkdir(exist_ok=False,parents=True)
    start=time.perf_counter();cpu=time.process_time()
    pre=ROOT/'PREREG_R2.json'
    assert hashlib.sha256(pre.read_bytes()).hexdigest()==pre.with_suffix('.sha256').read_text().strip()
    ref=load_reference()
    with np.load(ROOT/'r1/nominal_complete_v4/fields.npz') as f:
        z={k:f[k].copy() for k in f.files}
    days=z['days'];A,m=factors(z)
    Ar=np.interp(days,ref['days'],ref['C']);mr=np.interp(days,ref['days'],ref['HP_wound_reference_fraction'])
    baseline=score(days,z['strength_pct'])['RMSE_pp']
    factor_error=float(np.max(abs(75*A*m-z['strength_pct'])))
    cases={};curves={}
    def add(name,pred,scope,port):
        a=score(days,pred);a.update(case=name,port=port,delta_RMSE_pp=baseline-a['RMSE_pp'],scope=scope)
        cases[name]=a;curves[name]=pred
    add('baseline',z['strength_pct'],'Preserved R1 integrated proxy','none')
    for a,b in itertools.product([0,1],repeat=2):
        add(f'factor_A{a}_m{b}',75*(Ar if a else A)*(mr if b else m),
            'Diagnostic matched observation formula; R3 quantity inherits early strength calibration; HP observed input',
            f'abundance={"R3_C" if a else "chain_UIM"},maturity={"R3_HP" if b else "chain_M_fraction"}')
    cf=z['legacy_fields'][:,3]
    monitor=np.maximum(0,np.minimum(z['edges'][1:],.00025)-z['edges'][:-1])/.00025
    C=np.interp(days,z['legacy_field_days'],cf@monitor)
    add('current_FV_C_only',75*C*m,'Substitute preserved FV C (sparse saved times interpolated), retain UIM mature fraction','abundance')
    hp=np.interp(days,z['chemistry_days'],z['chemistry_D_H_A'][1])
    mhp=np.maximum(0,(hp-.003)/.040)
    add('DHA_HP_only',75*A*mhp,'Illustrative R3 HP observation normalization on isolated R4 H, held constant below day10; not identified native map','maturity')
    add('bridge_reference_identity',z['strength_pct'].copy(),'R1 and R3 both bridge=1; no differing reference input','bridge')
    add('angle_reference_identity',z['strength_pct'].copy(),'Isotropic UIM tensors make direction immaterial; R3 isotropic reference','load_direction')
    original=json.loads((PEER/'r1/native_chain/native_geometry_response.json').read_text())
    pp=original['parameters']
    changes={'gap':'gap_m','biological_width':'width_m','perfusion_width':'width_m',
             'face_conductance':'face_conductance','k_deposit':'k_deposit_day','k_mature_legacy':'k_mature_day'}
    base=scenario_config();compute_costs={}
    for key,sourcekey in changes.items():
        cfg=copy.deepcopy(base);cfg['parameters'][key]['value']=pp[sourcekey]
        d=run(cfg,mode='scenario',out=out/('upstream_'+key))
        pred=np.interp(days,d['strength_time']['value']['days'],d['strength_time']['value']['pct'])
        add('upstream_'+key,pred,'Single R3 port substitution; all other R1 ports retained',key)
        compute_costs[key]=d['costs']
    cfg=copy.deepcopy(base)
    for key,sourcekey in changes.items():cfg['parameters'][key]['value']=pp[sourcekey]
    d=run(cfg,mode='scenario',out=out/'upstream_all_R3')
    add('upstream_all_R3',np.asarray(d['strength_time']['value']['pct']),
        'All six R3 geometry/FV rate/face ports; chain UIM observer unchanged','six_upstream_ports')
    compute_costs['all_R3']=d['costs']
    # Swap early state/hazard, interpolating from its own registered spatial grid.
    with np.load(PEER/'r1/native_chain/native_geometry_response.npz') as f:
        early={'state':np.array([np.interp(z['x'],f['x'],v) for v in f['early_state']]),
               'hazard':np.interp(z['x'],f['x'],f['hazard'])}
    meta=json.loads((ROOT/'r1/nominal_complete_v4/checkpoint21.json').read_text())
    p=r.Params(**meta['parameters']);g=r.Grid(p,dx=5e-6)
    params={k:v['value'] for k,v in base['parameters'].items()}
    s=time.perf_counter();l=late_response(g,p,early,params,90.,.1)
    compute_costs['early_R3']={'wall_s':time.perf_counter()-s}
    add('early_state_hazard_R3',l['strength_pct'],
        'R3 early state and hazard mapped in physical x; R1 geometry, FV law, seal-boundary and observer retained','early_state_and_hazard')
    np.savez_compressed(out/'ATTRIBUTION_CURVES.npz',days=days,chain_abundance=A,chain_maturity=m,R3_abundance=Ar,
                        R3_maturity=mr,current_FV_C=C,**curves)
    f=lambda a,b:cases[f'factor_A{a}_m{b}']['RMSE_pp']
    shapley_A=.5*((f(0,0)-f(1,0))+(f(0,1)-f(1,1)))
    shapley_m=.5*((f(0,0)-f(0,1))+(f(1,0)-f(1,1)))
    interactions=f(0,0)-f(1,0)-f(0,1)+f(1,1)
    result={'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'baseline_RMSE_pp':baseline,
            'reference_RMSE_pp':score(ref['days'],ref['strength_pct'])['RMSE_pp'],
            'factorization_max_error_pp':factor_error,'factorial_interaction_RMSE_pp':interactions,
            'two_factor_Shapley_RMSE_reduction_pp':{'abundance':shapley_A,'maturity':shapley_m},
            'factor_check_days':[{ 'day':t,'chain_A':float(np.interp(t,days,A)),'chain_m':float(np.interp(t,days,m)),
               'R3_C':float(np.interp(t,days,Ar)),'R3_HP_fraction':float(np.interp(t,days,mr))} for t in [7,14,21,28,42,90]],
            'rows':list(cases.values()),'costs':compute_costs,'wall_s':time.perf_counter()-start,'CPU_s':time.process_time()-cpu,
            'strongest_control':'TIE: same conventional port interventions and factorial algebra',
            'lineage_correction':{'R3_k_deposit_day':pp['k_deposit_day'],'chain_k_deposit_day':.1,
              'R3_rate_fit':'response_r1.py main: least_squares to legacy strength3% at7 and20% at21, not new R3 strength fit',
              'R3_C_is_measured_collagen':False,'same_chemistry_input_claim':False,
              'chemistry_strength_coupling_in_R1':False},
            'gate':'PASS' if factor_error<1e-10 and abs(f(1,1)-8.306367094671637)<1e-8 else 'FAIL',
            'integrated_empirical_gate':'UNKNOWN','reference_proxy_gate':'FAIL'}
    write(out/'ATTRIBUTION.json',result)
    for row in result['rows']:print(row['case'],round(row['RMSE_pp'],8),round(row['delta_RMSE_pp'],8),flush=True)
    print('factor Shapley',result['two_factor_Shapley_RMSE_reduction_pp'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args().out)
