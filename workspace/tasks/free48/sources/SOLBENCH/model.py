"""Shared synthetic nonlinear region family anchored to frozen native DallaMan RHS.
State: organ12=native13 without Gt, regions reshape(4,M)=[free,bound,gate,exposure],ledger5.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
import ast,copy,hashlib,json,resource,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.sparse import csr_matrix,lil_matrix
R=Path(__file__).resolve().parent
scope=dict(np=np,solve_ivp=solve_ivp,brentq=brentq)
allowed={'P','IDX','NSTATE','k_empt','rhs','basal_state'};nodes=[]
for n in ast.parse((R/'inputs/native_glucose.py').read_text()).body:
 if isinstance(n,ast.FunctionDef) and n.name in allowed:nodes.append(n)
 elif isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in allowed for t in n.targets):nodes.append(n)
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(R/'inputs/native_glucose.py'),'exec'),scope)
P=scope['P'];P['HE_b']=-P['m5']*P['S_b']+P['m6'];P['m3_b_diag_only']=P['HE_b']*P['m1']/(1-P['HE_b']);MOL_MG=180156.;BW=P['BW'];VPT=P['V_G']*.1*1e-3*BW;VTT=VPT*P['k1']/P['k2'];CONV=MOL_MG*60/BW
KEEP=np.array([0,2,3,4,5,6,7,8,9,10,11,12]);YB=scope['basal_state']()[0]
DEFAULT=dict(mu=1.,diff=.08,A=.2,jstar=1e-4,K=5.,Bcap=20.,kon=.0008,koff=.02,tau_open=30.,tau_close=45.,KI=P['I_b'],KF=float(YB[1]),meal_mg=60000.)

def row_hash(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()
class Coupled:
 def __init__(self,n,params=None,partition=None,edits=None,edge_edits=None,extra_edges=None):
  self.n=int(n);self.params={**DEFAULT,**(params or {})};self.partition=np.arange(n) if partition is None else np.array(partition,int)
  if len(self.partition)!=n or min(self.partition)<0:raise ValueError('partition')
  self.m=int(max(self.partition)+1);self.count=np.bincount(self.partition,minlength=self.m);self.weights=self.count/n
  if np.any(self.count==0):raise ValueError('empty cluster')
  self.edits=copy.deepcopy(edits or {});self.edge_edits=copy.deepcopy(edge_edits or {});self.extra_edges=list(extra_edges or []);physical_m=np.ones(n);physical_cap=np.ones(n)
  for i,x in self.edits.items():physical_m[int(i)]=x.get('m',1.);physical_cap[int(i)]=x.get('cap_scale',1.)
  self.perfusion=np.zeros(self.m);self.cap_scale=np.zeros(self.m)
  for i in range(self.m):
   selected=self.partition==i
   if np.ptp(physical_m[selected]) or np.ptp(physical_cap[selected]):raise ValueError('heterogeneous constitutive law merged into one state')
   self.perfusion[i]=physical_m[selected][0];self.cap_scale[i]=physical_cap[selected][0]
  self.cap=self.params['Bcap']*self.params['mu']*self.cap_scale;self.kon=self.params['kon']/self.params['mu'];self.koff=self.params['koff']
  self.fac=CONV*self.params['A']*self.params['jstar']*self.params['mu'];self.K=self.params['K'];self.tauo=self.params['tau_open']*self.params['mu'];self.tauc=self.params['tau_close']/self.params['mu']
  rates=lil_matrix((self.m,self.m));self.edges=[]
  for i in range(n):
   j=(i+1)%n;factor=self.edge_edits.get(i,self.edge_edits.get(str(i),1.));self.edges.append((i,j,self.params['diff']*factor))
  self.edges.extend(self.extra_edges)
  for i,j,d in self.edges:
   a,b=self.partition[i],self.partition[j]
   if a!=b:
    rates[a,b]+=d/self.count[a];rates[a,a]-=d/self.count[a];rates[b,a]+=d/self.count[b];rates[b,b]-=d/self.count[b]
  self.L=rates.tocsr();self.dim=17+4*self.m
  self.source_row=dict(n=self.n,params=self.params,edits=self.edits,edge_edits=self.edge_edits,extra_edges=self.extra_edges,scope='SYNTHETIC_COUPLED_REGION_CLOSURE',empirical='UNKNOWN');self.source_id=row_hash(self.source_row)
 def regions(self,y):return np.asarray(y)[12:12+4*self.m].reshape(4,self.m)
 def initial_state(self):
  y=np.zeros(self.dim);y[:12]=YB[KEEP];y[3]=self.params['meal_mg'];z=self.regions(y);z[0]=YB[1];z[1]=.2*self.cap;z[2]=.3;return y
 def native_state(self,y):
  x=np.zeros(13);x[KEEP]=y[:12];x[1]=self.weights@self.regions(y)[0];return x
 def quantities(self,y):
  f,b,r,E=self.regions(y);o=y[:12];cb=o[0]*BW/(MOL_MG*VPT);ct=f*BW/(MOL_MG*VTT)
  sat=cb/(self.K+cb)-ct/(self.K+ct);q=self.perfusion*(P['k1']*o[0]-P['k2']*f)+self.fac*r*sat
  vm=(1-P['part'])*(P['Vm0']+P['VmX']*o[6]);uid=vm*f/(P['Km0']+f)
  bind=self.kon*f*(self.cap-b)-self.koff*b;I=o[1]/P['V_I'];ar=I/(self.params['KI']+I);fr=f/(self.params['KF']+f)
  gate=ar*(1-r)/self.tauo-fr*r/self.tauc
  return q,uid,bind,gate
 def organ(self,t,o,fmean,qmean,uidmean):
  x=np.zeros(13);x[KEEP]=o;x[1]=fmean;dx=scope['rhs'](t,x,self.params['meal_mg'],0,1,0,0,True)
  q0=P['k1']*o[0]-P['k2']*fmean;dx[0]+=q0-qmean;dx[10]+=P['K']*(q0-qmean)/P['V_G']
  egp=max(P['kp1']-P['kp2']*o[0]-P['kp3']*o[8]-P['kp4']*o[9],0);ra=P['f']*P['k_abs']*o[5]/BW if self.params['meal_mg']>0 else 0;ex=P['ke1']*max(o[0]-P['ke2'],0)
  return dx[KEEP],np.array([egp,ra,P['Uii'],ex,uidmean])
 def rhs(self,t,y):
  q,uid,bind,gate=self.quantities(y);z=self.regions(y);o,ledger=self.organ(t,y[:12],self.weights@z[0],self.weights@q,self.weights@uid)
  return np.r_[o,np.array([q-uid-bind+self.L@z[0],bind,gate,z[0]]).ravel(),ledger]
 def jac(self,t,y):
  z=self.regions(y);f,b,r,E=z;o=y[:12];q,uid,bind,gate=self.quantities(y);qmean=self.weights@q;umean=self.weights@uid;fmean=self.weights@f
  j=lil_matrix((self.dim,self.dim));base=np.r_[*self.organ(t,o,fmean,qmean,umean)]
  rows=np.r_[np.arange(12),np.arange(self.dim-5,self.dim)]
  for k in range(12):
   step=1e-6*(abs(o[k])+1);oo=o.copy();oo[k]+=step;dif=(np.r_[*self.organ(t,oo,fmean,qmean,umean)]-base)/step
   for a,v in zip(rows,dif):
    if v:j[a,k]=v
  cb=o[0]*BW/(MOL_MG*VPT);ct=f*BW/(MOL_MG*VTT);csb=BW/(MOL_MG*VPT);cst=BW/(MOL_MG*VTT);sat=cb/(self.K+cb)-ct/(self.K+ct)
  qgp=self.perfusion*P['k1']+self.fac*r*self.K/(self.K+cb)**2*csb;qf=-self.perfusion*P['k2']-self.fac*r*self.K/(self.K+ct)**2*cst;qr=self.fac*sat
  vm=(1-P['part'])*(P['Vm0']+P['VmX']*o[6]);uf=vm*P['Km0']/(P['Km0']+f)**2;ux=(1-P['part'])*P['VmX']*f/(P['Km0']+f)
  bf=self.kon*(self.cap-b);bb=-self.kon*f-self.koff;I=o[1]/P['V_I'];ar=I/(self.params['KI']+I);fr=f/(self.params['KF']+f)
  gf=-self.params['KF']/(self.params['KF']+f)**2*r/self.tauc;gr=-ar/self.tauo-fr/self.tauc;gip=self.params['KI']/(self.params['KI']+I)**2*(1-r)/self.tauo/P['V_I']
  ids=np.arange(self.m);fi=12+ids;bi=12+self.m+ids;ri=12+2*self.m+ids;ei=12+3*self.m+ids
  for k in range(self.m):
   a,c,h,ee=int(fi[k]),int(bi[k]),int(ri[k]),int(ei[k]);j[a,a]=qf[k]-uf[k]-bf[k];j[a,c]=-bb[k];j[a,h]=qr[k];j[a,0]=qgp[k];j[a,6]=-ux[k];j[c,a]=bf[k];j[c,c]=bb[k];j[h,a]=gf[k];j[h,h]=gr[k];j[h,1]=gip[k];j[ee,a]=1
   for orgrow,factor in [(0,1),(9,P['K']/P['V_G'])]:j[orgrow,a]=-factor*self.weights[k]*qf[k];j[orgrow,h]=-factor*self.weights[k]*qr[k]
   j[self.dim-1,a]=self.weights[k]*uf[k]
  ll=self.L.tocoo()
  for a,c,v in zip(ll.row,ll.col,ll.data):j[12+int(a),12+int(c)]+=v
  qdg=-self.weights@qgp;j[0,0]+=qdg;j[9,0]+=P['K']/P['V_G']*qdg;j[self.dim-1,6]=self.weights@ux
  return j.tocsr()
 def inventory(self,y):z=self.regions(y);return y[0]+self.weights@(z[0]+z[1])
 def exposure_concentration(self,y):return self.regions(y)[3]*BW/(MOL_MG*VTT)*60

def homogeneous_checkpoint(params=None,t=120.,rtol=1e-9):
 model=Coupled(1,params,partition=[0]);y=model.initial_state();sol=solve_ivp(model.rhs,(0,t),y,method='BDF',jac=model.jac,rtol=rtol,atol=1e-10,dense_output=True)
 if not sol.success:raise RuntimeError(sol.message)
 return dict(t_min=float(t),state=sol.y[:,-1].tolist(),params=model.params,source_id=model.source_id,history='fully coupled homogeneous0:120min; no imposedfuturefeaturehistory',history_sha256=row_hash({'params':model.params,'t':t,'state':sol.y[:,-1].tolist()}),cost_nfev=sol.nfev)
def lift_checkpoint(checkpoint,model):
 if checkpoint['params']!=model.params:raise ValueError('checkpoint source/history mismatch')
 src=np.array(checkpoint['state']);y=np.zeros(model.dim);y[:12]=src[:12];model.regions(y)[:]=src[12:16,None];y[-5:]=src[-5:];return y

def active_partition(n,nodes,radius):
 active=sorted({(int(i)+k)%n for i in nodes for k in range(-radius,radius+1)});p=np.full(n,len(active),int)
 if len(active)==n:return np.arange(n)
 for j,i in enumerate(active):p[i]=j
 return p

def integrate(model,y0,tspan=(120.,150.),sc_dose_U=0.,event_delay_min=0.,method='BDF',rtol=2e-8,atol=1e-9):
 wall=time.perf_counter();cpu=time.process_time();y0=np.array(y0,float).copy();before=model.inventory(y0);ledger0=y0[-5:].copy();pieces=[];event=tspan[0]+event_delay_min
 if not tspan[0]<=event<=tspan[1]:raise ValueError('event outside horizon')
 start=tspan[0];nfev=nlu=njev=0
 for stop,injection in [(event,True),(tspan[1],False)]:
  if stop>start:
   sol=solve_ivp(model.rhs,(start,stop),y0,method=method,jac=model.jac,rtol=rtol,atol=atol,dense_output=True)
   if not sol.success:raise RuntimeError(sol.message)
   pieces.append(sol);nfev+=sol.nfev;nlu+=sol.nlu;njev+=sol.njev;y0=sol.y[:,-1].copy();start=stop
  if injection:y0[11]+=sc_dose_U*P['pmol_per_unit']/BW
 def evaluate(t):
  aa=np.atleast_1d(t);vals=[]
  for ti in aa:
   pick=next((s for s in reversed(pieces) if s.t[0]-1e-12<=ti<=s.t[-1]+1e-12),None)
   if pick is None:raise ValueError('query outside solvedclock')
   vals.append(pick.sol(ti))
  result=np.array(vals).T;return result[:,0] if np.ndim(t)==0 else result
 end=y0;led=end[-5:]-ledger0;defect=model.inventory(end)-before-led[0]-led[1]+led[2]+led[3]+led[4]
 return SimpleNamespace(model=model,end=end,sol=evaluate,nfev=nfev,nlu=nlu,njev=njev,cost=dict(cpu_s=time.process_time()-cpu,wall_s=time.perf_counter()-wall,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,state_dim=model.dim,regions=model.m),ledger_error_mg_kg=float(defect))
