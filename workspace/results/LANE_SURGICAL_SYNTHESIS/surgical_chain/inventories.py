"""U/I/M joint chemistry-direction-bridge state driven by FV healing history.

Linear direction-independent kinetics plus an affine stress observation. No
identification from HP, tracer contrast or Γ. Exact leaf scope is synthetic.
"""
import numpy as np
from scipy.linalg import expm

def query_strength(qm,angle_deg,bridge_fraction=1.):
 if not 0<=bridge_fraction<=1:raise ValueError('bridge probability outside [0,1]')
 e=np.array([np.cos(np.deg2rad(angle_deg)),np.sin(np.deg2rad(angle_deg))])
 q=np.asarray(qm)
 if q.shape[:2]!=(2,2) or not np.isfinite(q).all():raise ValueError('Invalid tensor')
 if not np.allclose(q,np.swapaxes(q,0,1),atol=1e-12):raise ValueError('Tensor not symmetric')
 if np.linalg.eigvalsh(np.moveaxis(q,-1,0) if q.ndim==3 else q).min()<-1e-10:raise ValueError('Tensor not PSD')
 return 150*bridge_fraction*np.einsum('i,ij...,j->...',e,q,e)

def reservoir(d0,h0,a0,elapsed,k,eta,law='first_order'):
 if min(d0,h0,a0,k)<0 or np.any(np.asarray(elapsed)<0) or not 0<=eta<=1:raise ValueError('Invalid chemical inventory')
 if law=='first_order':d=d0*np.exp(-k*np.asarray(elapsed))
 elif law=='second_order':d=d0/(1+k*d0*np.asarray(elapsed))
 else:raise ValueError('Unknown chemical law')
 return np.array([d,h0+eta*(d0-d)/2,a0+(1-eta)*(d0-d)])

class JointInventory:
 def __init__(self,n,k1=.1,k2=.05,death=.03,initial=None):
  self.q=np.zeros((3,2,2,n)) if initial is None else np.array(initial,float,copy=True)
  self.k1,self.k2,self.death=k1,k2,death
  self.born=np.zeros(n);self.removed=np.zeros(n)
  self.initial=np.trace(self.q.sum(axis=0),axis1=0,axis2=1).copy()
  self.maximum_balance=0.
 def step(self,dt,birth,oxygen_gate):
  k1=self.k1*np.asarray(oxygen_gate);k2=self.k2*np.asarray(oxygen_gate);de=self.death
  def rhs(q):
   u,i,m=q
   source=np.zeros_like(u);source[0,0]=.5*birth;source[1,1]=.5*birth
   return np.array([source-(k1+de)*u,k1*u-(k2+de)*i,k2*i])
  old=self.q.copy();mid=old+.5*dt*rhs(old)
  self.q=old+dt*rhs(mid)
  self.born+=dt*birth
  self.removed+=dt*de*np.trace(mid[:2].sum(axis=0),axis1=0,axis2=1)
  mass=np.trace(self.q.sum(axis=0),axis1=0,axis2=1)
  self.maximum_balance=max(self.maximum_balance,float(abs(mass+self.removed-self.initial-self.born).max()))
  if np.linalg.eigvalsh(np.moveaxis(self.q,-1,0)).min()<-1e-10:raise RuntimeError('Joint inventory positivity broken; refine dt')
