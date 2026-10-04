"""Per-observable conventional calibration with no shared physical state."""
import numpy as np
from .ports import synth,unknown

class SeparateCurve:
 def __init__(self,days,values,unit,source):
  self.days=np.asarray(days,float);self.values=np.asarray(values,float);self.unit=unit;self.source=source
  if self.days.shape!=self.values.shape or np.any(np.diff(self.days)<=0):raise ValueError('Separate cohort curve must have strictly increasing times')
 def query(self,day,*,changed_inputs=False):
  if changed_inputs:return unknown(self.unit,self.source,'Independent curve has no intervention/history map; requires separately specified calibration/closure')
  if day<self.days[0] or day>self.days[-1]:return unknown(self.unit,self.source,'Extrapolation law not identified')
  return synth(float(np.interp(day,self.days,self.values)),self.unit,self.source,'Calibrated observable interpolation only, no heldout crossprediction; fitted points not independent evidence')
