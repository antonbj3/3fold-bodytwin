"""Dimensioned values with explicit empirical scope; null is never zero."""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math
from typing import Any
import numpy as np

STATUSES={"MEASURED",'SYNTETISKT',"UNKNOWN"}
class UnknownPort(ValueError):pass
@dataclass(frozen=True)
class Port:
 value: Any
 unit: str
 status: str
 uncertainty: Any
 source: str
 scope: str
 joint_id: str | None=None
 def __post_init__(self):
  if self.status not in STATUSES:raise ValueError('Invalid port status')
  if not self.unit or not self.source or not self.scope:raise ValueError('unit, source and scope required')
  if (self.value is None)!=(self.status=="UNKNOWN"):raise ValueError("null iff UNKNOWN; no zero-filled unknowns")
  if self.value is not None:
   a=np.asarray(self.value)
   if np.issubdtype(a.dtype,np.number) and not np.isfinite(a).all():raise ValueError('Nonfinite port')
 def require(self,unit,*,nonnegative=False,positive=False):
  if self.unit!=unit:raise ValueError(f'Unit mismatch: {self.unit} != {unit}')
  if self.status=="UNKNOWN":raise UnknownPort(self.scope)
  a=np.asarray(self.value)
  if nonnegative and np.any(a<0) or positive and np.any(a<=0):raise ValueError('Out-of-range physical port')
  return self.value
 def dict(self):return asdict(self)

def synth(value,unit,source,scope,uncertainty='Specified scenario, no CI',joint_id='scenario_nominal'):
 return Port(value,unit,'SYNTETISKT',uncertainty,source,scope,joint_id)
def unknown(unit,source,scope):return Port(None,unit,"UNKNOWN",None,source,scope)
def from_dict(d):return Port(**d)
def closure(value,unit,source,scope,parents):
 if any(p.status=="UNKNOWN" for p in parents):return unknown(unit,source,scope)
 return synth(value,unit,source,scope,{'kind':'propagated conditional scenario','parents':[p.source for p in parents],'posterior':None},parents[0].joint_id if parents else None)
