"""Graph-context histories and sparse future defect operator; synthetic laws.

The exact unedited prefix is periodic, not a supplied future input. Eight
contexts retain distinct neighborhoods even when there are only two materials.
"""
import copy,time
import numpy as np
from scipy.sparse import lil_matrix
import model as m
from phase_operator import PHASE_LAWS
WORD=(0,0,0,0,1,1,1,1)

class ContextCoupled(m.Coupled):
 def __init__(self,n,params=None,word=WORD,laws=None,edits=None,edge_edits=None,extra_edges=None):
  if n%len(word):raise ValueError('period must divide graph')
  self.word=tuple(word);self.context_laws=copy.deepcopy(laws or PHASE_LAWS);self.local_edits=copy.deepcopy(edits or {})
  physical={i:{**self.context_laws[word[i%len(word)]],**self.local_edits.get(i,self.local_edits.get(str(i),{}))} for i in range(n)}
  super().__init__(n,params,edits=physical,edge_edits=edge_edits,extra_edges=extra_edges)
  self.context_of_cluster=np.arange(n)%len(word);self.source_row.update(word=list(word),material_laws=self.context_laws,local_edits=self.local_edits,scope='SYNTHETIC_CONTEXT_NATIVE_REGIONS');self.source_id=m.row_hash(self.source_row)

class ContextPartition:
 def __init__(self,n,active,period):
  self.n=n;self.period=period;self.active=tuple(sorted(set(active)));self.lookup={i:k for k,i in enumerate(self.active)};self.bulk={};k=len(self.active)
  for p in range(period):
   if n//period>sum(i%period==p for i in self.active):self.bulk[p]=k;k+=1
  self.m=k
 def __getitem__(self,i):return self.lookup.get(int(i),self.bulk.get(int(i)%self.period))
 def __len__(self):return self.n
 def dense(self):return np.array([self[i] for i in range(self.n)])

def compile_context(n,centres,radius,params=None,word=WORD,laws=None,edits=None,edge_edits=None,extra_edges=None,retain_nodes=()):
 wall=time.perf_counter();cpu=time.process_time();period=len(word)
 if n%period:raise ValueError('period must divide graph')
 laws=copy.deepcopy(laws or PHASE_LAWS);edits=copy.deepcopy(edits or {});edge_edits=copy.deepcopy(edge_edits or {});extra_edges=list(extra_edges or [])
 active=set(retain_nodes)|{(int(c)+j)%n for c in centres for j in range(-radius,radius+1)}
 if not {int(i) for i in edits}<=active:raise ValueError('material edit must be retained')
 if not {int(i) for e in extra_edges for i in e[:2]}<=active:raise ValueError('chord endpoints must be retained')
 out=m.Coupled(1,params);out.n=n;out.word=tuple(word);out.context_laws=laws;out.local_edits=edits;out.partition=ContextPartition(n,active,period);out.m=out.partition.m;out.count=np.ones(out.m)
 out.context_of_cluster=list(i%period for i in out.partition.active)
 for p,k in out.partition.bulk.items():out.count[k]=n//period-sum(i%period==p for i in active);out.context_of_cluster.append(p)
 out.context_of_cluster=np.array(out.context_of_cluster);out.weights=out.count/n;out.perfusion=np.array([laws[word[p]]['m'] for p in out.context_of_cluster]);out.cap_scale=np.array([laws[word[p]]['cap_scale'] for p in out.context_of_cluster])
 for i,changes in edits.items():
  k=out.partition[int(i)];out.perfusion[k]=changes.get('m',out.perfusion[k]);out.cap_scale[k]=changes.get('cap_scale',out.cap_scale[k])
 out.cap=out.params['Bcap']*out.params['mu']*out.cap_scale;out.dim=17+4*out.m;rates=lil_matrix((out.m,out.m));edge_ids={i for node in active for i in [node,(node-1)%n]};visited=[]
 def add(a,b,d):
  if a!=b:
   rates[a,b]+=d/out.count[a];rates[a,a]-=d/out.count[a];rates[b,a]+=d/out.count[b];rates[b,b]-=d/out.count[b]
 for i in sorted(edge_ids):
  d=out.params['diff']*edge_edits.get(i,edge_edits.get(str(i),1.));visited.append((i,(i+1)%n,d));add(out.partition[i],out.partition[(i+1)%n],d)
 interface_counts=[]
 for p in range(period):
  remaining=n//period-sum(i%period==p for i in edge_ids);interface_counts.append(remaining)
  if remaining:
   d=remaining*out.params['diff']
   for key,factor in edge_edits.items():
    if int(key) not in edge_ids and int(key)%period==p:d+=out.params['diff']*(factor-1)
   add(out.partition.bulk[p],out.partition.bulk[(p+1)%period],d)
 for i,j,d in extra_edges:add(out.partition[i],out.partition[j],d)
 out.L=rates.tocsr();out.edits=edits;out.edge_edits=edge_edits;out.extra_edges=extra_edges;out.edges=visited+extra_edges;out.interface_counts=interface_counts
 out.source_row=dict(n=n,params=out.params,word=list(word),material_laws=laws,local_edits=edits,edge_edits=edge_edits,extra_edges=extra_edges,scope='SYNTHETIC_CONTEXT_NATIVE_REGIONS',empirical='UNKNOWN');out.source_id=m.row_hash(out.source_row)
 out.construction_cost=dict(cpu_s=time.process_time()-cpu,wall_s=time.perf_counter()-wall,physical_nodes_visited=len(active),physical_edges_visited=len(visited)+len(extra_edges),bulk_interface_counts=interface_counts,materialized_N_partition=False)
 return out

def context_checkpoint(params=None,t=120,word=WORD,laws=None):
 model=ContextCoupled(len(word),params,word,laws);res=m.integrate(model,model.initial_state(),tspan=(0,t),rtol=1e-10,atol=1e-11)
 return dict(state=res.end.tolist(),t_min=t,params=model.params,word=list(word),material_laws=model.context_laws,history_sha256=m.row_hash(res.end.tolist()),native_feedback_evolved=True,cost=res.cost)

def lift_context(cp,model):
 if cp['params']!=model.params or cp['word']!=list(model.word) or cp['material_laws']!=model.context_laws:raise ValueError('context history/source mismatch')
 z=np.array(cp['state']);P=len(model.word);regions=z[12:12+4*P].reshape(4,P);y=np.zeros(model.dim);y[:12]=z[:12];y[-5:]=z[-5:];model.regions(y)[:]=regions[:,model.context_of_cluster];return y

def transfer_context_history(old,state,new):
 if old.n!=new.n or old.params!=new.params or old.word!=new.word or old.context_laws!=new.context_laws:raise ValueError('physical source/history event map needed')
 z=old.regions(state);y=np.zeros(new.dim);y[:12]=state[:12];y[-5:]=state[-5:];out=new.regions(y)
 for node in new.partition.active:out[:,new.partition[node]]=z[:,old.partition[node]]
 for p,k in new.partition.bulk.items():
  mask=old.context_of_cluster==p;total=z[:,mask]@old.count[mask];taken=out[:,[new.partition[i] for i in new.partition.active if i%len(new.word)==p]].sum(axis=1);out[:,k]=(total-taken)/new.count[k]
 error=float(new.inventory(y)-old.inventory(state))
 if abs(error)>1e-8 or np.any(out[1]>new.cap+1e-7):raise ValueError('inventory/capacity physical event required')
 return y,dict(inventory_error=error,history_replayed=False,source_ids=[old.source_id,new.source_id],history_hashes=[m.row_hash(state.tolist()),m.row_hash(y.tolist())])
