from __future__ import annotations
from pathlib import Path
from dataclasses import asdict
import hashlib,json,time,resource
import numpy as np
from .ports import Port,from_dict,closure,unknown
from .path_preference import incision_path_gate
from .hemostasis import hydraulic,bleeding
from .inventories import JointInventory,query_strength,reservoir
from .vendor import response_r1 as r
from .vendor.storage_r1 import early_storage
from .closure_manifest import manifest

UNITS={'chemistry_D10':'mol/molreference collagen','chemistry_H10':'mol/molreference collagen','chemistry_A10':'mol/molreference collagen','chemistry_rate':'1/day','chemistry_yield':'1','seal_rate':'1/min','gap':'m','biological_width':'m','perfusion_width':'m','radius_bins':'m','density_bins':'1/m²','driving_pressure':'Pa','blood_viscosity':'Pa s','vascular_path':'m','domain':'m','grid_dx':'m','early_dt':'min','early_until':'min','late_dt':'day','late_until':'day','oxygen_storage':'mlO2/(m³ Torr)','oxygen_permeability':'mlO2/(m min Torr)','oxygen_consumption':'mlO2/(m³ min)','face_conductance':'mlO2/(m² min Torr)','face_pressure':'Torr','seal_oxygen_coupling':'1','aggregation_factor':'1','k_deposit':'1/day','k_mature_legacy':'1/day','k_U_I':'1/day','k_I_M':'1/day','turnover_U_I':'1/day','birth_scale':'1/day','bridge_fraction':'1','load_angle':'degree'}
REVIEW='PENDING_INDEPENDENT_REVIEW'

def write(path,obj):
 with path.open('x') as f:json.dump(obj,f,indent=2,ensure_ascii=False,allow_nan=False,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item())

def check_config(config):
 ports={k:from_dict(config['parameters'][k]) for k in UNITS}
 for k,u in UNITS.items():
  if ports[k].unit!=u:raise ValueError('Unit mismatch for '+k)
 for k in ['path','cut_depth']:ports[k]=from_dict(config[k])
 if ports['path'].unit!='m' or ports['cut_depth'].unit!='m':raise ValueError('Invalid geometry units')
 ports['edge_radius']=from_dict(config['tool']['edge_radius']);ports['speed']=from_dict(config['tool']['speed'])
 for k,u in [('edge_radius','m'),('speed','m/s')]:
  if ports[k].unit!=u:raise ValueError('Tool units')
 for j,layer in enumerate(config['layers']):
  for name,u in [('thickness','m'),('gamma_cut','J/m²')]:
   p=from_dict(layer[name]);ports[f'layer{j}.{name}']=p
   if p.unit!=u:raise ValueError('Layer units')
 return ports

def incision(config,ports):
 path=np.asarray(ports['path'].require('m'),float)
 if path.ndim!=2 or path.shape[1]!=3 or path.shape[0]<2:raise ValueError('Polyline needs at least two xyz points')
 length=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum());depth=float(ports['cut_depth'].require('m',nonnegative=True))
 force=0.;remaining=depth;parents=[ports['path'],ports['cut_depth']];layerwork=[]
 for j,layer in enumerate(config['layers']):
  t=float(ports[f'layer{j}.thickness'].require('m',positive=True));cut=min(t,remaining);remaining-=cut
  gp=ports[f'layer{j}.gamma_cut'];gamma=float(gp.require('J/m²',nonnegative=True));force+=gamma*cut
  parents.extend([ports[f'layer{j}.thickness'],gp]);layerwork.append({'layer':layer['name'],'cut_thickness_m':cut,'work_J':gamma*cut*length})
 if remaining>1e-12:raise ValueError('Cut deeper than tissue stack')
 port=lambda v,u,n:closure(v,u,'incision single crack-plane ledger',n,parents).dict()
 return dict(length=port(length,'m','Polyline length'),depth=port(depth,'m','Cut depth'),area=port(length*depth,'m²','One crack-plane reference area; not two faces'),force=port(force if length>0 else 0.,'N','Cleavage component only; total tool load/contact/release UNKNOWN'),work=port(force*length,'J','Cleavage work only'),layer_work=layerwork,
   total_tool_work=unknown('J','INCISION final port','Contact, prestress release and pre-tip work unknown').dict(),
   gap=ports['gap'].dict(),biological_width=ports['biological_width'].dict(),perfusion_width=ports['perfusion_width'].dict(),
   tool_law='Edge radius/speed are recorded but have no calibrated mapping to gamma/gap/biological injury; those are supplied ports')

def late_response(g,p,e,params,until,dt,*,checkpoint=None):
 kwargs=dict(k1=params['k_U_I'],k2=params['k_I_M'],death=params['turnover_U_I'])
 inv=JointInventory(g.n,initial=None if checkpoint is None else checkpoint['Q'],**kwargs)
 start=0. if checkpoint is None else checkpoint['absolute_day']
 if checkpoint is not None:
  inv.born=checkpoint['born'].copy();inv.removed=checkpoint['removed'].copy();inv.initial=checkpoint['initial'].copy()
 late_days=[];strength=[];inventory=[];saved={}
 def observe(day,Y,P):
  gate=np.minimum(P/(10+P)/(p.healthy_Torr/(10+p.healthy_Torr)),1.)
  birth=params['birth_scale']*Y[2]*gate*(1-Y[3])
  inv.step(dt,birth,gate)
  monitor=np.maximum(0,np.minimum(g.edges[1:],p.width_m)-g.edges[:-1])/p.width_m
  strength.append(float(query_strength(inv.q[2],params['load_angle'],params['bridge_fraction'])@monitor))
  late_days.append(day);inventory.append([float(np.trace(x,axis1=0,axis2=1)@monitor) for x in inv.q])
  if abs(day-21)<dt/3:
   saved.update(absolute_day=day,Q=inv.q.copy(),born=inv.born.copy(),removed=inv.removed.copy(),initial=inv.initial.copy(),long_state=Y.copy())
 late=r.long(g,p,e,until=until,initial=None if checkpoint is None else checkpoint['long_state'],start_day=start,observer=observe)
 return dict(days=np.asarray(late_days),strength_pct=np.asarray(strength),UIM=np.asarray(inventory),Q_final=inv.q,born=inv.born,removed=inv.removed,initial=inv.initial,checkpoint=saved,legacy=late,max_mass_balance=inv.maximum_balance)

def run(config,*,mode='evidence',out=None,grid_kind='graded',shear=True,inventory='legacy'):
 if inventory not in ['legacy','shared']:raise ValueError('Unknown inventory model')
 if mode not in ['evidence','scenario']:raise ValueError('Mode must be evidence or scenario')
 st=time.perf_counter();cpu=time.process_time();ports=check_config(config)
 path_gate=incision_path_gate(config.get('avlankning'),mode=mode)
 path_blocked=path_gate is not None and path_gate['blocks_straight_chain']
 unavailable=[k for k,p in ports.items() if p.status!="MEASURED"]
 unknown_names=[k for k,p in ports.items() if p.status=="UNKNOWN"]
 if mode=='evidence' or unknown_names or path_blocked:
  reason='Synthetic/unknown inputs cannot populate an empirical prediction' if mode=='evidence' else 'Scenario must explicitly fill all consumed unknowns'
  if path_blocked and mode=='scenario' and not unknown_names:reason='Requested preference does not license straight continuation; a branched geometry and growth solver is required'
  result={'review_state':REVIEW,'scientific_admission':False,'mode':mode,'engineering_gate':'PASS_UNKNOWN_PROPAGATION','empirical_joint_gate':'UNKNOWN','missing_or_nonempirical_ports':unavailable,'native_strength':unknown('%intact','Final source ports',reason).dict(),'stages':{s:{'status':"UNKNOWN",'value':None,'reason':reason} for s in ['injury_gap','bleeding_hemostasis','oxygen','healing_inventories','strength']}}
  if path_gate is not None:
   result['incision_path_gate']=path_gate
   result['stages']['path_continuation']={'status':"UNKNOWN",'value':None,'reason':"Straight continuation requires REVIEW in every requested box; deflection needs a branched geometry solver"}
   if path_blocked and mode=='scenario' and not unknown_names:
    result['engineering_gate']='PASS_PATH_GUARD'
    result['missing_or_nonempirical_ports']=unavailable+['straight_path_preference']
  if out is not None:out.mkdir(parents=True,exist_ok=False);write(out/'summary.json',result)
  return result
 positive={'domain','grid_dx','early_dt','early_until','late_dt','late_until','oxygen_storage','oxygen_permeability','oxygen_consumption','blood_viscosity','vascular_path','biological_width'}
 params={k:p.require(UNITS[k],nonnegative=k!='load_angle',positive=k in positive) for k,p in ports.items() if k in UNITS}
 for k in ['seal_oxygen_coupling','bridge_fraction','chemistry_yield']:
  if not 0<=params[k]<=1:raise ValueError(k+' outside [0,1]')
 if params['perfusion_width']>params['domain'] or params['biological_width']>params['domain']:raise ValueError('width exceeds domain')
 # Enforce integer-step horizon; no different inventory step hidden behind FV rounding.
 for h,d in [('early_until','early_dt'),('late_until','late_dt')]:
  if not np.isclose(params[h]/params[d],round(params[h]/params[d]),atol=1e-9):raise ValueError('Horizon must be a multiple of dt')
 mech=incision(config,ports);length=mech['length']['value'];depth=mech['depth']['value'];area=mech['area']['value']
 if path_gate is not None:mech['path_gate']=path_gate
 hist=hydraulic(params['radius_bins'],params['density_bins'],area,params['driving_pressure'],params['blood_viscosity'],params['vascular_path'])
 bstart=time.perf_counter();blood=bleeding(hist,minutes=params['early_until'],shear=shear,aggregation_factor=params['aggregation_factor'],seal_rate=params['seal_rate']);blood_s=time.perf_counter()-bstart
 p=r.Params(domain_m=params['domain'],width_m=params['perfusion_width'],gap_m=params['gap'],cut_length_m=length,cut_depth_m=depth,
  permeability=params['oxygen_permeability'],consumption=params['oxygen_consumption'],face_Torr=params['face_pressure'],face_conductance=0.,early_dt_min=params['early_dt'],early_minutes=params['early_until'],dt_day=params['late_dt'],k_deposit_day=params['k_deposit'],k_mature_day=params['k_mature_legacy'])
 if p.width_m<=0:raise ValueError('Positive perfusion_width needed by normalized FV wound monitor')
 g=r.Grid(p,grid_kind,dx=params['grid_dx'])
 def face(t):
  seal=np.interp(t,blood['t_min'],blood['seal'])
  return params['face_conductance']*(1-params['seal_oxygen_coupling']*seal)
 early=early_storage(g,p,params['oxygen_storage'],face_history=face)
 # Biological observation window is supplied independently of vascular width.
 pp=r.replace(p,width_m=params['biological_width'])
 # Legacy late boundary: sealed terminal face conductance, not re-opened naked wound.
 pp=r.replace(pp,face_conductance=face(params['early_until']))
 solver=late_response
 if inventory=='shared':
  from .shared import late_response as solver
 ls=time.perf_counter();late=solver(g,pp,early,params,params['late_until'],params['late_dt']);late_s=time.perf_counter()-ls
 allparents=list(ports.values());output=lambda v,u,s:closure(v,u,'surgical_chain coupled conditional closure',s,allparents).dict()
 chem_days=late['days'][late['days']>=10]
 chem=reservoir(params['chemistry_D10'],params['chemistry_H10'],params['chemistry_A10'],chem_days-10,params['chemistry_rate'],params['chemistry_yield'])
 result={'inherited_closures':manifest(),'chemistry_D_H_A':output({'days':chem_days.tolist(),'D_H_A':chem.tolist()},'mol/molreference collagen','Isolated RESPONSE R4 autonomous branch; fixedC10reference; no measured native chemistry→rupture coupling'),'review_state':REVIEW,'scientific_admission':False,'mode':'scenario','empirical_joint_gate':'UNKNOWN','strongest_control_gate':'TIE','cost_gain_claimed':False,
  'incision':mech,'bleeding':{k:output(v,'m³/s' if 'flow' in k else 'µL' if 'volume' in k else '1', 'Pressure/shear-driven synthetic bleed') for k,v in blood['summary'].items() if k in ['initial_flow_m3_s','blood_volume_uL','seal_final']},
  'hemostasis_time':output(blood['summary']['hemostasis_min'],'min','Joint coverage+thrombin; finite horizon censoring') if blood['summary']['hemostasis_min'] is not None else unknown('min','hemostasis closure','Not closed in chosen horizon; right-censored').dict(),
  'oxygen_edge_pressure':output(float(early['pressure'][0]),'Torr','At120min, FV surface-adjacent cell'),
  'hypoxic_width_proxy':output(float(np.sum(g.dx*(early['pressure']<10))),'m','10Torr threshold is synthetic; not viability/clinical necrosis'),
  'inventory_UIM_final':output(late['UIM'][-1].tolist(),'normalized reference mass','U/I/M observer driven by F and O2; not measured collagen mass'),
  'strength_time':output({'days':late['days'].tolist(),'pct':late['strength_pct'].tolist()},'%intact','Affine synthetic75%reference; chemistry-direction-bridge moments, not native rupture law'),
  'native_rupture_strength':unknown('Pa','RESPONSE final ports','Native chemistry→traction→rupture law absent').dict(),
  'numerics':{'cells':g.n,'oxygen_inventory':early['inventory'],'early_oxygen':early['oxygen'],'late_oxygen':late['legacy']['oxygen'],'UIM_balance':late['max_mass_balance'],'fibrinogen_error_nM':blood['summary']['max_fibrinogen_inventory_error_nM']},
  'blood_diagnostics':blood['summary'],'defaults_consumed':config,
  'coupling_scope':'Seal changes atmospheric face conductance, not cut perfusion. Early O2/Q036 hazard feeds FV healing; FV C feeds O2 permeability; new U/I/M directional marks are downstream observers and do not replace FV C feedback. No tool→injury law or unique numerical gain claimed.',
  'costs':{'blood_wall_s':blood_s,'late_wall_s':late_s,'full_wall_s':time.perf_counter()-st,'CPU_s':time.process_time()-cpu,'max_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'acquisition_wall_s':None,'setup_import_included_in_full_wall':False}}
 if inventory=='shared':
  result['shared_reference']=late['shared_reference']
  result['numerics']['shared_reference']=late['shared_reference']
  result['chemical_strength_diagnostic']=output({'days':late['days'].tolist(),'pct':late['chemical_strength_pct'].tolist()},'%intact','Shared birth/removal chemistry; synthetic HP-to-strength scale, native rupture UNKNOWN')
  result['coupling_scope']='Same C=trace(U+I+M) controls O2 permeability and strength marks; shared formation/removal. D/H/A marks share that collagen birth/removal; original isolated fixedC10 chemistry remains separate. No dynamic external HP or abundance inputs.'
 if out is not None:
  out.mkdir(parents=True,exist_ok=False);write(out/'summary.json',result)
  np.savez_compressed(out/'fields.npz',x=g.x,edges=g.edges,early_state=early['state'],hazard=early['hazard'],early_pressure=early['pressure'],early_series=early['series'],chemistry_days=chem_days,chemistry_D_H_A=chem,blood_t_min=blood['t_min'],coverage=blood['coverage'],seal=blood['seal'],seal_bins=blood['seal_bins'],blood_m3=blood['blood_m3'],days=late['days'],strength_pct=late['strength_pct'],UIM=late['UIM'],Q_final=late['Q_final'],legacy_days=late['legacy']['days'],legacy_strength_pct=late['legacy']['strength_pct'],legacy_fields=late['legacy']['fields'],legacy_field_days=late['legacy']['saved_days'])
  if inventory=='shared':np.savez_compressed(out/'shared_reference.npz',days=late['days'],C=late['C'],UIM=late['UIM'],born_removed=late['ledger'],chemical_marks=late['chemical_marks_series'],chemical_strength_pct=late['chemical_strength_pct'],born=late['born'],removed=late['removed'],initial=late['initial'])
  if late['checkpoint']:
   ck=late['checkpoint'];np.savez_compressed(out/'checkpoint21.npz',**ck,early_state=early['state'],hazard=early['hazard'],edges=g.edges,x=g.x,chemistry_D_H_A=reservoir(params['chemistry_D10'],params['chemistry_H10'],params['chemistry_A10'],11.,params['chemistry_rate'],params['chemistry_yield']))
   write(out/'checkpoint21.json',{'inventory_model':inventory,'review_state':REVIEW,'scientific_admission':False,'absolute_day':ck['absolute_day'],'parameters':asdict(pp),'inventory_rates':{k:params[k] for k in ['k_U_I','k_I_M','turnover_U_I','birth_scale','bridge_fraction','load_angle']},'full_config':config,'sufficient_for':'This same synthetic late FV+observer closure; requires all stored state/hazard, absolute clock, rates/boundaries, grid and reference; not an empirical history certificate'})
 return result

def resume(checkpoint_path,out,until=90.):
 cp=Path(checkpoint_path);meta=json.loads(cp.with_suffix('.json').read_text());cfg=meta['full_config'];params={k:from_dict(cfg['parameters'][k]).value for k in UNITS};p=r.Params(**meta['parameters']);g=r.Grid(p,dx=params['grid_dx'])
 with np.load(cp) as data:ck={k:data[k].copy() for k in data.files}
 if not np.array_equal(g.edges,ck['edges']):raise ValueError('Checkpoint grid mismatch')
 ck['absolute_day']=float(ck['absolute_day']);e={'hazard':ck['hazard'],'state':ck['early_state']}
 solver=late_response
 if meta.get('inventory_model','legacy')=='shared':
  from .shared import late_response as solver
 late=solver(g,p,e,params,until,p.dt_day,checkpoint=ck)
 chem_days=late['days'][late['days']>=ck['absolute_day']]
 if 'chemistry_D_H_A' in ck:
  chem=reservoir(*ck['chemistry_D_H_A'],chem_days-ck['absolute_day'],params['chemistry_rate'],params['chemistry_yield'])
 else:chem=np.empty((3,0))
 out.mkdir(parents=True,exist_ok=False);np.savez_compressed(out/'resumed.npz',days=late['days'],chemistry_days=chem_days,chemistry_D_H_A=chem,strength_pct=late['strength_pct'],UIM=late['UIM'],Q_final=late['Q_final'],long_state=late['legacy']['state'])
 if 'shared_reference' in late:np.savez_compressed(out/'shared_reference_resumed.npz',days=late['days'],C=late['C'],chemical_marks=late['chemical_marks_series'],chemical_strength_pct=late['chemical_strength_pct'],born_removed=late['ledger'])
 write(out/'summary.json',{'review_state':REVIEW,'empirical_joint_gate':'UNKNOWN','prefix_replayed':False,'global_spatial_suffix_recomputed':True,'mass_balance':late['max_mass_balance'],'wall_s':late['legacy']['wall_s']})
 return late
