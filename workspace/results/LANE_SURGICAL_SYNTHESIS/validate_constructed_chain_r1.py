import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import copy,json,time,hashlib
import numpy as np
from surgical_chain.chain import run,write
from surgical_chain.scenario import scenario_config
R=Path(__file__).resolve().parent
pre={'operation':'Refine constructed chain and test causal seal ablation with shared scenario directions','frozen_before_compute':True,'cases':['nominal final code','half early/late steps','half front dx','seal-to-O2 coupling0','shared direction-1','shared direction+1'],'joint_scenario_directions':'Not a distribution or posterior: sign jointly modifies demand±10%,radius±5%,injury/perfusion width±40um; no probabilities or CI','budget':'pressure.5Torr; strength1pp for numerical refinements; no empirical bound','strongest_control':'Identical conventional FV/cohort/chemical/plug history operations:TIE','checkpoint':'Sufficient state includes full late FV state,earlyhazard,threeQ,ledger totals,DHA,grid,absoluteclock,law/rates/reference'}
p=R/'PREREG_R1_REFINEMENT.json'
write(p,pre);p.with_suffix('.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest()+'\n')
base=scenario_config();cases={}
for name in ['nominal_complete_v4','refine_time_v1','refine_space_v1','no_seal_oxygen_v1','joint_minus_v1','joint_plus_v1']:
 cfg=copy.deepcopy(base)
 if name=='refine_time_v1':
  for k in ['early_dt','late_dt']:cfg['parameters'][k]['value']/=2
 if name=='refine_space_v1':cfg['parameters']['grid_dx']['value']/=2
 if name=='no_seal_oxygen_v1':cfg['parameters']['seal_oxygen_coupling']['value']=0.
 if name.startswith('joint_'):
  sign=-1 if 'minus' in name else 1
  cfg['parameters']['oxygen_consumption']['value']*=1+.1*sign
  cfg['parameters']['radius_bins']['value']=[x*(1+.05*sign) for x in cfg['parameters']['radius_bins']['value']]
  for k in ['biological_width','perfusion_width']:cfg['parameters'][k]['value']+=sign*40e-6
  for x in cfg['parameters'].values():x['joint_id']=name
 start=time.perf_counter();d=run(cfg,mode='scenario',out=R/'r1'/name);d['costs']['runner_wall_s']=time.perf_counter()-start
 cases[name]=d
 print(name,d['strength_time']['value']['pct'][-1],flush=True)
nom=cases['nominal_complete_v4'];raw=np.load(R/'r1/nominal_complete_v4/fields.npz');checks=[]
for name in ['refine_time_v1','refine_space_v1']:
 d=cases[name];a=np.load(R/'r1'/name/'fields.npz')
 err=max(abs(np.interp(raw['days'],a['days'],a['strength_pct'])-raw['strength_pct']))
 p_err=abs(d['oxygen_edge_pressure']['value']-nom['oxygen_edge_pressure']['value'])
 checks.append({'case':name,'strength_pp':float(err),'pressure_Torr':p_err,'gate':'PASS' if err<=1 and p_err<=.5 else 'FAIL'})
abl=np.load(R/'r1/no_seal_oxygen_v1/fields.npz')
pre_mask=raw['early_series'][:,0]<=4.1
causal_error=float(abs(raw['early_series'][pre_mask,:6]-abl['early_series'][pre_mask,:6]).max())
assert causal_error<1e-9
r4=json.loads((R/'sources/LANE_SURGICAL_RESPONSE/PREREG_R4_GLYCO_RESERVOIR.json').read_text())
res={'review_state':'PENDING_INDEPENDENT_REVIEW','numerical_refinements':checks,'causal_preseal_max_error':causal_error,
 'seal_oxygen_ablation_strength90_pp':nom['strength_time']['value']['pct'][-1]-cases['no_seal_oxygen_v1']['strength_time']['value']['pct'][-1],
 'joint_scenario_strength90_values':{k:cases[k]['strength_time']['value']['pct'][-1] for k in ['joint_minus_v1','nominal_complete_v4','joint_plus_v1']},'scenario_range_is_CI':False,
 'stage_solver_wall_sum_s':sum(d['costs']['full_wall_s'] for d in cases.values()),'stage_solver_CPU_sum_s':sum(d['costs']['CPU_s'] for d in cases.values()),'case_costs':{k:d['costs'] for k,d in cases.items()},'gate':'PASS' if all(x['gate']=='PASS' for x in checks) else 'FAIL','empirical_joint_gate':'UNKNOWN','strongest_control_gate':'TIE'}
write(R/'REFINEMENT_R1.json',res)
print(json.dumps(res,indent=2))
