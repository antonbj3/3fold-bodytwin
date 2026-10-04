"""Frozen MSK coagulation, radius-resolved hydraulic/shear ports.

Pressure-driven independent tubes, uniform concentrations, synthetic plug-to-
lumen closure. Growing platelet coverage does not repair transected perfusion.
"""
from pathlib import Path
import ast
import numpy as np
from scipy.integrate import solve_ivp
V=Path(__file__).resolve().parent/'vendor'

def definitions(path,stop):
 tree=ast.parse(path.read_text());nodes=[]
 for n in tree.body:
  if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id==stop for t in n.targets):break
  if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef)) or isinstance(n,ast.Assign) and not any(isinstance(t,ast.Name) and t.id=='OUT_DIR' for t in n.targets):nodes.append(n)
 ns={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),ns)
 return ns
CO=definitions(V/'coagulation_hemostasis.py.source','sol_cat')
PL=definitions(V/'platelet_hemostasis.py.source','_baseline_probe')

def hydraulic(radius,density,area,pressure,viscosity,length):
 r,n=np.asarray(radius,float),np.asarray(density,float)
 if r.shape!=n.shape or r.ndim!=1 or not np.isfinite(r).all() or not np.isfinite(n).all():raise ValueError('Invalid vessel histogram')
 if np.any(r<=0) or np.any(n<0) or min(area,pressure)<0 or min(viscosity,length)<=0:raise ValueError('Invalid hydraulic ports')
 q=np.pi*pressure*r**4/(8*viscosity*length)
 gamma=pressure*r/(2*viscosity*length)
 return dict(vessel_count=n*area,Q_per_bin_m3_s=n*area*q,gamma_per_bin_s_inverse=gamma,
             initial_flow_m3_s=float(n@q*area))

def shear_on_rate(gamma,mode='two'):
 # Exact closure algebra from frozen swarm source; sigmoid probabilities are
 # chosen proxies, not measurements of catalytic adhesion rates.
 x=np.log10(np.maximum(np.asarray(gamma),np.finfo(float).tiny))
 z=np.log(19.);w=(np.log10(20000.)-np.log10(6000.))/(2*z)
 center=(np.log10(6000.)+np.log10(20000.))/2
 high=1/(1+np.exp(np.clip(-(x-center)/w,-700,700)))
 wint=(np.log10(900.)-np.log10(600.))/(2*z)
 low=1/(1+np.exp(np.clip((x-np.log10(np.sqrt(600.*900.)))/wint,-700,700)))
 if mode=='constant':return np.full_like(x,PL['K_ON0'])
 if mode=='two':return PL['K_ON0']*(low+(1-low)*high)
 raise ValueError('Unknown shear closure')

def bleeding(hist,*,minutes=120.,shear=True,aggregation_factor=1.,rtol=1e-9,seal_rate=20.):
 bins=len(hist['gamma_per_bin_s_inverse']);par=dict(CO['LOCKED'],TF=1000.,XIIa=0.,S_amp=5.)
 y0=np.r_[CO['y0_vec'](),np.zeros(3*bins+1)]
 q0=hist['Q_per_bin_m3_s'];gam=hist['gamma_per_bin_s_inverse'];threshold=PL['THETA_C']
 def rhs(t,y):
  a,p=y[15:15+bins],y[15+bins:15+2*bins]
  cov=PL['W_A']*a+PL['W_P']*p
  sealed=y[15+2*bins:15+3*bins]
  lumen=np.maximum(.05,1-cov) if shear else np.ones(bins)
  ka=shear_on_rate(gam*lumen,'two' if shear else 'constant')
  clot=y[CO['IDX']['IIa']]>=CO['CLOT_THRESHOLD_NM']
  leak=np.maximum(1-cov/threshold,0)**2
  leak=np.where(clot & (cov>=threshold),0.,leak)
  flow=q0*lumen**4*leak*(1-sealed)
  dseal=seal_rate*(1-sealed)*(clot & (cov>=threshold))
  return np.r_[CO['rhs'](t,y[:15],par),ka*(1-a)-PL['KOFF_A']*a,
                PL['K_AGG0']*aggregation_factor*a*(1-p)-PL['KOFF_P']*p,dseal,60*flow.sum()]
 t=np.unique(np.r_[np.linspace(0,min(10,minutes),1001),np.linspace(min(10,minutes),minutes,221)])
 sol=solve_ivp(rhs,(0,minutes),y0,method='LSODA',rtol=rtol,atol=1e-11,t_eval=t,max_step=.1)
 if not sol.success:raise RuntimeError(sol.message)
 cov=PL['W_A']*sol.y[15:15+bins]+PL['W_P']*sol.y[15+bins:15+2*bins]
 clot=sol.y[CO['IDX']['IIa']]>=CO['CLOT_THRESHOLD_NM']
 shut=(cov>=threshold)&clot
 times=[float(t[np.flatnonzero(a)[0]]) if a.any() else None for a in shut]
 active=q0>0
 hits=np.flatnonzero(np.all(shut[active],axis=0)) if active.any() else np.array([0])
 weights=q0/q0.sum() if q0.sum()>0 else np.zeros(bins)
 seal_bins=sol.y[15+2*bins:15+3*bins]
 seal=np.sum(weights[:,None]*seal_bins,axis=0)
 fgn=sol.y[CO['IDX']['Fg']]+sol.y[CO['IDX']['Fbn']]
 return dict(t_min=t,coverage=cov,seal=seal,seal_bins=seal_bins,blood_m3=sol.y[-1],thrombin_nM=sol.y[CO['IDX']['IIa']],
  summary=dict(initial_flow_m3_s=float(q0.sum()),blood_volume_uL=float(sol.y[-1,-1]*1e9),
  hemostasis_min=float(t[hits[0]]) if len(hits) else None,per_bin_hemostasis_min=times,
  censor_horizon_min=minutes,seal_final=float(seal[-1]),
  max_fibrinogen_inventory_error_nM=float(abs(fgn-fgn[0]).max()),rhs_evaluations=sol.nfev,
  interpretation='Synthetic pressure-driven tubes; censoring is not a prediction of permanent nonclosure; seal does not restore cut vessels'))
