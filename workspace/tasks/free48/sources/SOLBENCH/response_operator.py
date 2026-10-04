"""Implicit local partition and conservative retained-history continuation.
Does not materialize N-state histories or an N-length partition during local updates.
Error qualification remains external; ordinary local ROM alone is not a discovery.
"""
import copy,time
import numpy as np
from scipy.sparse import lil_matrix
import model as m
class SparsePartition:
 def __init__(self,n,active):
  self.n=n;self.active=tuple(sorted(set(active)));self.lookup={node:k for k,node in enumerate(self.active)};self.bulk=len(self.active) if len(self.active)<n else None
 def __len__(self):return self.n
 def __getitem__(self,node):
  if np.ndim(node):return np.array([self[int(x)] for x in node])
  return self.lookup.get(int(node),self.bulk)
 def dense(self):return np.array([self[i] for i in range(self.n)])

def compile_local(n,centres,radius,params=None,edits=None,edge_edits=None,extra_edges=None,retain_nodes=()):
 wall=time.perf_counter();cpu=time.process_time();edits=copy.deepcopy(edits or {});edge_edits=copy.deepcopy(edge_edits or {});extra_edges=list(extra_edges or [])
 active={int(i)%n for i in retain_nodes}|{(int(c)+j)%n for c in centres for j in range(-radius,radius+1)}
 if not {int(i) for i in edits}<=active:raise ValueError('an edited constitutive law cannot be represented by a homogeneous bulk')
 # Scalar initialization computes native/source coefficients, not an N-field scan.
 out=m.Coupled(1,params);out.n=n;out.partition=SparsePartition(n,active);out.m=len(active)+(len(active)<n);out.count=np.ones(out.m)
 if out.partition.bulk is not None:out.count[-1]=n-len(active)
 out.weights=out.count/n;out.edits=edits;out.edge_edits=edge_edits;out.extra_edges=extra_edges;out.perfusion=np.ones(out.m);out.cap_scale=np.ones(out.m)
 for node,law in edits.items():idx=out.partition[int(node)];out.perfusion[idx]=law.get('m',1.);out.cap_scale[idx]=law.get('cap_scale',1.)
 out.cap=out.params['Bcap']*out.params['mu']*out.cap_scale;out.dim=17+4*out.m;rates=lil_matrix((out.m,out.m));out.edges=[]
 edge_ids={i for node in active for i in [node,(node-1)%n]}
 for i in sorted(edge_ids):out.edges.append((i,(i+1)%n,out.params['diff']*edge_edits.get(i,edge_edits.get(str(i),1.))))
 out.edges.extend(extra_edges)
 for i,j,d in out.edges:
  a,b=out.partition[i],out.partition[j]
  if a!=b:rates[a,b]+=d/out.count[a];rates[a,a]-=d/out.count[a];rates[b,a]+=d/out.count[b];rates[b,b]-=d/out.count[b]
 out.L=rates.tocsr();out.source_row=dict(n=n,params=out.params,edits=edits,edge_edits=edge_edits,extra_edges=extra_edges,scope='SYNTHETIC_COUPLED_REGION_CLOSURE',empirical='UNKNOWN');out.source_id=m.row_hash(out.source_row);out.construction_cost=dict(cpu_s=time.process_time()-cpu,wall_s=time.perf_counter()-wall,physical_nodes_visited=len(active),physical_edges_visited=len(out.edges),materialized_N_partition=False)
 return out

def transfer_history(old,state,new):
 if old.n!=new.n or old.params!=new.params:raise ValueError('state history/source changes need an explicit event map')
 start=time.process_time();z=old.regions(state);y=np.zeros(new.dim);y[:12]=state[:12];y[-5:]=state[-5:];out=new.regions(y);active=new.partition.active
 for node in active:out[:,new.partition[node]]=z[:,old.partition[node]]
 if new.partition.bulk is not None:
  out[:,new.partition.bulk]=(z@old.count-out[:,:len(active)].sum(axis=1))/new.count[new.partition.bulk]
 physical=np.max(np.abs(new.inventory(y)-old.inventory(state)))
 if physical>1e-8:raise RuntimeError('conservative state transfer failed')
 if np.any(out[1]>new.cap+1e-7):raise ValueError('site-capacity shrink requires a declared bound release event')
 provenance=dict(parent_source_id=old.source_id,new_source_id=new.source_id,parent_state_hash=m.row_hash(np.asarray(state).tolist()),mapped_state_hash=m.row_hash(y.tolist()),inventory_transfer_error_mg_kg=float(physical),old_region_states=old.m,new_region_states=new.m,CPU_s=time.process_time()-start,old_history_replayed=False)
 return y,provenance
