"""Nonlinear noisy-port prediction pivot; same synthetic physics, new validation."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='2'
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from plan_measurements_r2 import ROOT, RANGES, GROUPS, config, extract, contraction
from surgical_chain.chain import run, write

OUT=ROOT/'r2/nonlinear_information_v1'
STORE=Path('external_mount')
BANDWIDTHS=[.25,.5,1.,2.,4.]

def prior():
    pre=ROOT/'PREREG_R2_NONLINEAR_INFO.json'
    assert hashlib.sha256(pre.read_bytes()).hexdigest()==pre.with_suffix('.sha256').read_text().strip()
    u=qmc.Sobol(d=len(RANGES)+1,scramble=True,seed=6124).random_base2(6)
    x=np.empty((64,len(RANGES)))
    for j,(key,(lo,hi)) in enumerate(RANGES.items()):
        s=u[:,j+1]
        if key in ['biological_width','perfusion_width']:s=.65*u[:,0]+.35*s
        x[:,j]=lo+(hi-lo)*s
    with (ROOT/'NONLINEAR_INFORMATION_PRIOR_R2.npz').open('xb') as f:
        np.savez_compressed(f,drivers=u,parameters=x,names=np.array(list(RANGES)))
    OUT.mkdir(parents=True,exist_ok=False)
    write(OUT/'DESIGN.json',{'new_samples':64,'seed':6124,'measurement_noise_seed':6125,
        'prior':'Identical frozen synthetic R2 ranges and shared-width drivers, fresh Sobol scramble',
        'bandwidths':BANDWIDTHS,'training_noisy_LOO_draws':8,'test_noise_draws':32,
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False})

def build(start,stop):
    st=time.perf_counter();cpu=time.process_time()
    with np.load(ROOT/'NONLINEAR_INFORMATION_PRIOR_R2.npz') as f:
        pars=f['parameters'].copy();names=f['names'].tolist()
    for i in range(start,stop):
        p=dict(zip(names,pars[i].tolist()));dest=STORE/f'case{i:03d}'
        d=run(config(p),mode='scenario',out=dest);y,z=extract(d,dest/'fields.npz',p)
        write(OUT/f'sample{i:03d}.json',{'sample':i,'parameters':p,'outputs':y,'measurements':z,
            'solver_costs':d['costs'],'raw_fields':str(dest/'fields.npz')})
        if (i-start)%8==0:print('new validation',i,y,flush=True)
    write(OUT/f'BATCH_{start}_{stop}.json',{'wall_s':time.perf_counter()-st,'CPU_s':time.process_time()-cpu})

def distances(x,t):
    return np.maximum(0,np.sum(x*x,1)[:,None]+np.sum(t*t,1)[None,:]-2*x@t.T)

def kernel_predict(d,y,h):
    a=-d/(2*h*h);a-=np.max(a,axis=1)[:,None]
    weights=np.exp(a);weights/=weights.sum(1)[:,None]
    return weights@y,weights

def fit():
    st=time.perf_counter();cpu=time.process_time()
    train=[json.loads((ROOT/f'r2/measurement_v1/sample{i:03d}.json').read_text()) for i in range(256)]
    y=np.array([r['outputs'] for r in train]);sd=y.std(0,ddof=1);rng=np.random.default_rng(6126)
    results=[];models={}
    for key,g in GROUPS.items():
        z=np.array([[r['measurements'][k] for k in g['obs']] for r in train])
        q=contraction(y,z,g['sigma'],g['common_sigma'],disconnected=g.get('disconnected',False))
        R=q['R'];mu=z.mean(0)
        W=np.linalg.inv(np.linalg.cholesky(np.atleast_2d(np.cov(z,rowvar=False))+R))
        x=(z-mu)@W.T
        noisy=z[:,None,:]+rng.multivariate_normal(np.zeros(len(mu)),R,size=(len(z),8))
        queries=(noisy.reshape(-1,len(mu))-mu)@W.T
        d=distances(queries,x);exclude=np.repeat(np.arange(len(z)),8)
        d[np.arange(len(d)),exclude]=np.inf
        labels=np.repeat(y,8,axis=0);cols=[0,1,4,5];losses=[]
        for h in BANDWIDTHS:
            pred,_=kernel_predict(d,y,h)
            losses.append(float(np.mean(((pred[:,cols]-labels[:,cols])/sd[cols])**2)))
        h=BANDWIDTHS[int(np.argmin(losses))]
        models[key]={'z':z,'mu':mu,'W':W,'R':R,'B':q['B'],'bandwidth':h,'disconnected':g.get('disconnected',False)}
        results.append({'id':key,'bandwidth':h,'noisy_LOO_normalized_MSE':dict(zip(map(str,BANDWIDTHS),losses)),
            'disconnected':g.get('disconnected',False)})
    with (OUT/'TRAINED_MODELS_V2.npz').open('xb') as f:
        np.savez_compressed(f,y=y,**{key+'_'+k:v for key,m in models.items() for k,v in m.items()})
    write(OUT/'TRAINING_V2.json',{'rows':results,'trained_without_new64_outputs':True,'noise_seed':6126,
        'wall_s':time.perf_counter()-st,'CPU_s':time.process_time()-cpu,'review_state':'PENDING_INDEPENDENT_REVIEW'})
    print('bandwidths',[(r['id'],r['bandwidth']) for r in results],flush=True)

def evaluate():
    st=time.perf_counter();cpu=time.process_time()
    test=[json.loads((OUT/f'sample{i:03d}.json').read_text()) for i in range(64)]
    yt=np.array([r['outputs'] for r in test]);rng=np.random.default_rng(6125);rows=[];raw={}
    with np.load(OUT/'TRAINED_MODELS_V2.npz') as f:
        y=f['y'];baseline=np.mean((yt-y.mean(0))**2,axis=0)
        for key,g in GROUPS.items():
            z=np.array([[r['measurements'][k] for k in g['obs']] for r in test])
            R=f[key+'_R'];noise=rng.multivariate_normal(np.zeros(len(g['obs'])),R,size=(64,32))
            noisy=(z[:,None,:]+noise).reshape(-1,len(g['obs']))
            labels=np.repeat(yt,32,axis=0)
            mu,W=f[key+'_mu'],f[key+'_W'];x=(f[key+'_z']-mu)@W.T
            lin=y.mean(0)+(noisy-mu)@f[key+'_B'].T
            if bool(f[key+'_disconnected']):pred=np.repeat(y.mean(0)[None,:],len(noisy),axis=0)
            else:pred,_=kernel_predict(distances((noisy-mu)@W.T,x),y,float(f[key+'_bandwidth']))
            mse=np.mean((pred-labels)**2,axis=0);mse_lin=np.mean((lin-labels)**2,axis=0)
            red=np.divide(baseline-mse,baseline,out=np.zeros_like(baseline),where=baseline>1e-20)
            red_lin=np.divide(baseline-mse_lin,baseline,out=np.zeros_like(baseline),where=baseline>1e-20)
            q={'id':key,'baseline_RMS_prediction_error':np.sqrt(baseline).tolist(),
                'kernel_RMS_prediction_error':np.sqrt(mse).tolist(),'linear_RMS_prediction_error':np.sqrt(mse_lin).tolist(),
                'kernel_MSE_reduction_fraction':red.tolist(),'linear_MSE_reduction_fraction':red_lin.tolist(),
                'workload_units':g['cost'],'kernel_value_per_cost':float(np.mean(red[[0,1,2,4]]))/g['cost'],
                'kernel_vs_prior_hypoxia_gate':'PASS' if mse[5]<=baseline[5]+1e-10 else 'FAIL',
                'strength_injury_vs_linear_gate':'PASS' if np.mean(mse[[0,1,4]]/baseline[[0,1,4]])<=1.1*np.mean(mse_lin[[0,1,4]]/baseline[[0,1,4]]) else 'FAIL'}
            rows.append(q);raw[key+'_kernel_predictions']=pred;raw[key+'_linear_predictions']=lin;raw[key+'_noisy_measurements']=noisy
    gate='PASS' if all(q['kernel_vs_prior_hypoxia_gate']=='PASS' and q['strength_injury_vs_linear_gate']=='PASS' for q in rows if q['id'] in ['M1','M2']) else 'FAIL'
    with (OUT/'VALIDATION_RAW.npz').open('xb') as f:np.savez_compressed(f,test_truth=yt,**raw)
    d={'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'gate':gate,'new_validation_samples':64,
       'training_samples':256,'noise_draws_per_sample':32,'rows':rows,
       'rank':sorted(rows,key=lambda q:q['kernel_value_per_cost'],reverse=True),
       'native_uncertainty':None,'uncertainty_interpretation':'Out-of-sample RMSpredictionrisk with noisy observations, not a native posterior SD',
       'strongest_control_gate':'TIE: conventional identical Gaussian positive-kernel regression receives same inputs',
       'new64_solver_wall_s':sum(r['solver_costs']['full_wall_s'] for r in test),
       'new64_solver_CPU_s':sum(r['solver_costs']['CPU_s'] for r in test),
       'evaluation_wall_s':time.perf_counter()-st,'evaluation_CPU_s':time.process_time()-cpu}
    write(OUT/'VALIDATION.json',d)
    for q in rows:print(q['id'],q['kernel_MSE_reduction_fraction'],q['kernel_vs_prior_hypoxia_gate'],q['strength_injury_vs_linear_gate'],flush=True)
    print('pivot gate',gate,'rank',[r['id'] for r in d['rank']],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prior','build','fit','evaluate'])
    p.add_argument('--start',type=int,default=0);p.add_argument('--stop',type=int,default=32);a=p.parse_args()
    if a.action=='prior':prior()
    elif a.action=='build':build(a.start,a.stop)
    elif a.action=='fit':fit()
    else:evaluate()
