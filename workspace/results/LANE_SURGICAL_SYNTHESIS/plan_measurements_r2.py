"""Noisy acquisition planning on the unchanged nonlinear synthetic R1 model.

Numbers are a designed prior and optimal-linear prediction MSE estimates,
not a native posterior, biological precision, or laboratory cost estimate.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[k]='2'
import argparse, copy, hashlib, json, time, resource
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from surgical_chain.chain import run, write
from surgical_chain.scenario import scenario_config

ROOT=Path(__file__).resolve().parent
STORE=Path('/mnt/games-240/research/bodytwin_solnight/LANE_SURGICAL_SYNTHESIS/r2')
RANGES={
 'gap':(100e-6,400e-6),'biological_width':(100e-6,700e-6),'perfusion_width':(100e-6,700e-6),
 'radius_factor':(.8,1.2),'density_factor':(.7,1.3),'driving_pressure':(2000.,6000.),
 'blood_viscosity':(.0028,.0042),'vascular_path':(.0007,.0013),
 'oxygen_consumption':(1000.,2000.),'oxygen_permeability':(.8e-6,1.8e-6),'oxygen_storage':(10.,100.),
 'face_conductance':(.001,.005),'face_pressure':(100.,160.),'seal_oxygen_coupling':(.5,1.),
 'aggregation_factor':(.6,1.4),'seal_rate':(10.,30.),'k_deposit':(.07,.13),
 'k_U_I':(.06,.14),'k_I_M':(.025,.075),'turnover_U_I':(.015,.045),'birth_scale':(.07,.13),
 'bridge_fraction':(.5,1.),'strength_reference_factor':(.7,1.3),
 'chemistry_D10':(.04,.09),'chemistry_H10':(.014,.022),'chemistry_rate':(.04,.1),
 'chemistry_yield':(.15,.45),'gamma_cut':(150.,380.)}
OUTPUTS=['strength28_pp','strength90_pp','restricted_hemostasis120_min','closed_by120_indicator',
         'injury_proxy_width_um','hypoxia_proxy_width_um']
# Cost units are transparent normalized workload assumptions, not quoted prices.
GROUPS={
 'M1':{'title':'Registered gap, independent injury/perfusion widths, radius/density/pressure/flow', 'cost':8.,
       'obs':['gap','biological_width','perfusion_width','radius_factor','density_factor','driving_pressure','blood_viscosity','vascular_path'],
       'sigma':[10e-6,25e-6,25e-6,.05,.1,200.,.000175,.00005], 'common_sigma':[0,10e-6,10e-6,.02,.02,80.,0,0]},
 'M2':{'title':'Oxygen/storage/boundary flux with fixed reference volume', 'cost':3.,
       'obs':['pressure20','pressure100','pressure300','face_flux','oxygen_storage','face_pressure'],
       'sigma':[2.,2.,2.,.02,6.,2.], 'common_sigma':[1.,1.,1.,.005,0,1.]},
 'M3':{'title':'Early radius-matched plug and persistent-seal dynamics', 'cost':5.,
       'obs':['coverage8_small','coverage8_large','seal8','initial_flow','aggregation_factor','seal_rate'],
       'sigma':[.03,.03,.05,1e-10,.1,2.], 'common_sigma':[.01,.01,.02,0,0,0]},
 'M4':{'title':'Absolute precursor/HP chemistry and reference collagen (missing observer map)', 'cost':12.,
       'obs':['chemistry_D10','chemistry_H10','chemistry_rate','chemistry_yield'],
       'sigma':[.005,.002,.007,.03], 'common_sigma':[.002,.001,0,0],
       'disconnected':True,'unrepresented':['GG/G/free subdivision','catalytic LOX map','native collagen mass to dimensionless UIM abundance']},
 'M5':{'title':'Chemistry-age-connected bridge and independent traction/ultimate-load port', 'cost':16.,
       'obs':['trace_M21','bridge_fraction','strength_reference_factor','strength21'],
       'sigma':[.015,.05,.05,1.], 'common_sigma':[.005,0,.02,.5],
       'scope':'Conditional upper-information design: assumes assays directly identify existing synthetic observer ports; native map UNKNOWN'},
 'M6':{'title':'Full force/work/new-area/interface and terminal-work ledger', 'cost':20.,
       'obs':['gamma_cut'],'sigma':[30.],'common_sigma':[0.],'disconnected':True,
       'unrepresented':['native tool-work to cell injury/perfusion transfer','mode-specific newly detached area','contact/release/pre-tip work']}}

def design():
    p=ROOT/'MEASUREMENT_DESIGN_R2.json'
    u=qmc.Sobol(d=len(RANGES)+1,scramble=True,seed=6122).random_base2(8)
    params=np.empty((256,len(RANGES)))
    names=list(RANGES)
    for j,(key,(lo,hi)) in enumerate(RANGES.items()):
        s=u[:,j+1]
        if key in ['biological_width','perfusion_width']:s=.65*u[:,0]+.35*s
        params[:,j]=lo+(hi-lo)*s
    description={'schema':'synthetic-noisy-measurement-planning-r2','frozen_before_ensemble':True,
      'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'ranges':RANGES,'parameter_names':names,
      'prior':'Uniform Sobol drivers on displayed bounds; shared width driver .65; all other drivers independent BY ASSUMPTION',
      'seed':6122,'n':256,'measurements':GROUPS,'outputs':OUTPUTS,
      'cost_unit':'Assumed relative combined preparation/instrument/analysis workload; nominal and0.5x..2x; actualhours/money UNKNOWN',
      'measurement_noise':'Independent Gaussian precision budgets plus one shared calibration error per package: R=diag(sigma^2)+common_sigma outer product',
      'rank_score':'mean fractional optimal-linear variance reduction of strength28,strength90,restrictedhemostasis,injuryproxy / workload',
      'error_floor':'Optional independent unresolved observation-map SD20pp in each strength output and50um in injury; sensitivity scenario, not identified discrepancy',
      'native_prediction_uncertainty':{'strength28':None,'strength90':None,'uncensored_hemostasis':None,'edge_necrosis':None},
      'M4_M6_policy':'Zero gain refers only to absent implemented output paths; assay value after future transfer-law identification remains UNKNOWN',
      'M5_policy':'Direct observation of synthetic Q/bridge/reference is a planning scenario, not an available calibrated native assay',
      'unmapped_subobservations':'Every missing subobservable explicit; no fabricated Hb advection, native necrosis or tool-injury law'}
    write(p,description)
    with (ROOT/'MEASUREMENT_PRIOR_R2.npz').open('xb') as f:np.savez_compressed(f,parameters=params,drivers=u,names=np.array(names))
    p.with_suffix('.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest()+'\n')
    print('design frozen',p,flush=True)

def get_params(i):
    with np.load(ROOT/'MEASUREMENT_PRIOR_R2.npz') as f:return dict(zip(f['names'].tolist(),f['parameters'][i].tolist()))

def config(par):
    c=scenario_config()
    for key,value in par.items():
        if key in c['parameters']:c['parameters'][key]['value']=value
    c['parameters']['radius_bins']['value']=[x*par['radius_factor'] for x in c['parameters']['radius_bins']['value']]
    c['parameters']['density_bins']['value']=[x*par['density_factor'] for x in c['parameters']['density_bins']['value']]
    for layer in c['layers']:layer['gamma_cut']['value']=par['gamma_cut']
    return c

def extract(d,fields,par):
    with np.load(fields) as f:
        y=f['strength_pct']*par['strength_reference_factor'] # bridge already in R1 query
        T=d['hemostasis_time']['value'];closed=T is not None
        injury=float(np.sum(np.diff(f['edges'])*(1-np.exp(-f['hazard'])>.20))*1e6)
        outputs=[float(np.interp(t,f['days'],y)) for t in [28,90]]+[float(T if closed else 120.),float(closed),injury,
                   float(d['hypoxic_width_proxy']['value']*1e6)]
        obs=dict(par)
        for t in [20,100,300]:obs['pressure'+str(t)]=float(np.interp(t*1e-6,f['x'],f['early_pressure']))
        obs['face_flux']=d['numerics']['oxygen_inventory']['terminal_face_flux_ml_m2_min']
        obs['coverage8_small']=float(np.interp(8,f['blood_t_min'],f['coverage'][0]))
        obs['coverage8_large']=float(np.interp(8,f['blood_t_min'],f['coverage'][1]))
        obs['seal8']=float(np.interp(8,f['blood_t_min'],f['seal']))
        obs['initial_flow']=d['blood_diagnostics']['initial_flow_m3_s']
        obs['trace_M21']=float(np.interp(21,f['days'],f['UIM'][:,2]))
        obs['strength21']=float(np.interp(21,f['days'],y))
    return outputs,obs

def batch(start,stop,out):
    out.mkdir(exist_ok=True,parents=True)
    begin=time.perf_counter();cpu=time.process_time();rows=[]
    for i in range(start,stop):
        par=get_params(i);c=config(par);dest=STORE/'ensemble_v1'/f'case{i:03d}'
        d=run(c,mode='scenario',out=dest)
        y,obs=extract(d,dest/'fields.npz',par)
        row={'sample':i,'parameters':par,'outputs':y,'measurements':obs,'solver_costs':d['costs'],
             'raw_fields':str(dest/'fields.npz'),'max_UIM_balance':d['numerics']['UIM_balance'],
             'oxygen_conservation':d['numerics']['oxygen_inventory']['relative_conservation_error']}
        write(out/f'sample{i:03d}.json',row);rows.append(row)
        if (i-start)%8==0:print('sample',i,'outputs',y,flush=True)
    write(out/f'BATCH_{start:03d}_{stop:03d}.json',{'samples':[r['sample'] for r in rows],
          'wall_s':time.perf_counter()-begin,'CPU_s':time.process_time()-cpu,'max_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})

def contraction(y,z,sigma,common,*,disconnected=False,noise_factor=1.,floor=None):
    cy=np.atleast_2d(np.cov(y,rowvar=False,ddof=1));cz=np.atleast_2d(np.cov(z,rowvar=False,ddof=1))
    yz=(y-y.mean(0)).T@(z-z.mean(0))/(len(y)-1)
    if disconnected:yz[:]=0.
    R=noise_factor**2*(np.diag(np.array(sigma)**2)+np.outer(common,common))
    # Scale measurement coordinates to avoid conditioning on heterogeneous units.
    scale=np.sqrt(np.diag(cz+R));zs=np.diag(1/scale)
    gain=(yz@zs)@np.linalg.solve(zs@(cz+R)@zs,(yz@zs).T)
    var=np.maximum(0,np.diag(cy));post=np.maximum(0,var-np.diag(gain))
    if floor is not None:var=var+np.array(floor)**2;post=post+np.array(floor)**2
    shrink=np.divide(var-post,var,out=np.zeros_like(var),where=var>1e-20)
    return {'prior_SD':np.sqrt(var).tolist(),'conditional_linear_SD':np.sqrt(post).tolist(),
            'variance_reduction_fraction':shrink.tolist(),
            'score':float(np.mean(shrink[[i for i in [0,1,2,4] if i<len(shrink)]])), 'B':np.linalg.solve(zs@(cz+R)@zs,(yz@zs).T).T@zs,
            'R':R}

def summary(out):
    rows=[json.loads((out/f'sample{i:03d}.json').read_text()) for i in range(256)]
    y=np.array([r['outputs'] for r in rows]);measurements={k:np.array([r['measurements'][k] for r in rows]) for k in rows[0]['measurements']}
    variants={};details={}
    for n,noise,floor in [(256,1.,None),(128,1.,None),(256,.5,None),(256,2.,None),(256,1.,[20,20,0,0,50,0])]:
        label=f'n{n}_noise{noise}_floor{floor is not None}';reports=[]
        for key,g in GROUPS.items():
            z=np.column_stack([measurements[k] for k in g['obs']])[:n]
            q=contraction(y[:n],z,g['sigma'],g['common_sigma'],disconnected=g.get('disconnected',False),noise_factor=noise,floor=floor)
            B,R=q.pop('B'),q.pop('R')
            # Heldout linear estimator diagnostic: train first128, test next128, include future noise.
            if n==256 and noise==1. and floor is None:
                train=contraction(y[:128],z[:128],g['sigma'],g['common_sigma'],disconnected=g.get('disconnected',False))
                pred=y[:128].mean(0)+(z[128:]-z[:128].mean(0))@train['B'].T
                mse=np.mean((y[128:]-pred)**2,axis=0)+np.diag(train['B']@train['R']@train['B'].T)
                base=np.mean((y[128:]-y[:128].mean(0))**2,axis=0)
                q['heldout_train128_test128_MSE_reduction_fraction']=np.divide(base-mse,base,out=np.zeros_like(base),where=base>1e-20).tolist()
                details[key]=[]
                for j,k in enumerate(g['obs']):
                    sub=contraction(y,measurements[k][:,None],[g['sigma'][j]],[g['common_sigma'][j]],disconnected=g.get('disconnected',False))
                    sub.pop('B');sub.pop('R');sub.update(observation=k,sigma=g['sigma'][j],common_sigma=g['common_sigma'][j]);details[key].append(sub)
            q.update(id=key,title=g['title'],workload_units=g['cost'],value_per_cost=q['score']/g['cost'],
                     value_per_cost_range=[q['score']/(2*g['cost']),q['score']/(.5*g['cost'])])
            reports.append(q)
        variants[label]=sorted(reports,key=lambda a:a['value_per_cost'],reverse=True)
    result={'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,
      'scope':'Optimal-linear MSE planning on a nonlinear synthetic ensemble; posteriorvariance interpretation requires Gaussian moment approximation',
      'native_outcomes':{'strength28_uncertainty':None,'strength90_uncertainty':None,'uncensored_hemostasis_uncertainty':None,'edge_necrosis_uncertainty':None},
      'output_names':OUTPUTS,'sample_count':256,'prior_mean':y.mean(0).tolist(),'prior_SD':y.std(0,ddof=1).tolist(),
      'closure_by120_probability_under_assumed_prior':float(y[:,3].mean()),'censored_samples':int(np.sum(y[:,3]==0)),
      'variants':variants,'individual_observations':details,'same_information_control_gate':'TIE',
      'costs':{'solver_wall_s':sum(r['solver_costs']['full_wall_s'] for r in rows),
               'solver_CPU_s':sum(r['solver_costs']['CPU_s'] for r in rows),
               'max_RSS_KiB':max(r['solver_costs']['max_RSS_KiB'] for r in rows),
               'batch_end_to_end':[json.loads(p.read_text()) for p in sorted(out.glob('BATCH_*.json'))],
               'laboratory_cost':None,'development_wall_s':None},
      'max_mass_balance':max(r['max_UIM_balance'] for r in rows),
      'max_oxygen_conservation':max(r['oxygen_conservation'] for r in rows)}
    write(out/'MEASUREMENT_VALUE.json',result)
    for key,rr in variants.items():print(key,[(q['id'],round(q['value_per_cost'],5)) for q in rr],flush=True)

def refinement(out):
    rows=[]
    for i in [0,63,127,255]:
        par=get_params(i);c=config(par)
        for key in ['early_dt','late_dt']:c['parameters'][key]['value']/=2
        dest=STORE/'refinement_v1'/f'case{i:03d}'
        d=run(c,mode='scenario',out=dest);y,_=extract(d,dest/'fields.npz',par)
        old=json.loads((out/f'sample{i:03d}.json').read_text())
        delta=abs(np.array(y)-old['outputs'])
        rows.append({'sample':i,'absolute_output_differences':delta.tolist(),'costs':d['costs'],
            'gate':'PASS' if max(delta[:2])<=1 and delta[2]<=1 and delta[4]<=25 else 'FAIL'})
    write(out/'REFINEMENT.json',{'rows':rows,'gate':'PASS' if all(x['gate']=='PASS' for x in rows) else 'FAIL'})
    print(rows,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['design','batch','summary','refinement'])
    p.add_argument('--start',type=int,default=0);p.add_argument('--stop',type=int,default=32)
    p.add_argument('--store',type=Path,default=STORE)
    p.add_argument('--out',type=Path,default=ROOT/'r2/measurement_v1');a=p.parse_args()
    allowed=Path('/mnt/games-240/research/bodytwin_solnight/LANE_SURGICAL_SYNTHESIS').resolve()
    if not (a.store.resolve().is_relative_to(allowed) or a.store.resolve().is_relative_to(ROOT)):
        raise ValueError('Field output must remain inside the authorized lane')
    STORE=a.store.resolve()
    if a.action=='design':design()
    elif a.action=='batch':batch(a.start,a.stop,a.out)
    elif a.action=='summary':summary(a.out)
    else:refinement(a.out)
