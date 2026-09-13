"""Bracketed inversion of an explicit two-wavelength pulse path model.

Owned path arrays, no runtime initialization or implicit artifact writes.
This solves a supplied optical model; it is not a clinical calibration.
"""
import numpy as np


class PulsePathInverse:
    def __init__(self, paths, background, fraction, extinction, *, hb_g_dl=15., pulse_fraction=.01):
        self.paths=np.array(paths,dtype=float,copy=True)
        self.background=np.array(background,dtype=float,copy=True)
        self.fraction=np.array(fraction,dtype=float,copy=True)
        self.extinction=np.array(extinction,dtype=float,copy=True)
        if self.paths.ndim!=2 or len(self.paths)==0 or self.paths.shape[1]==0:
            raise ValueError('Nonempty region paths required')
        n=self.paths.shape[1]
        if self.background.shape!=(n,) or self.fraction.shape!=(n,) or self.extinction.shape!=(2,2):
            raise ValueError('Region coefficients and two extinction rows required')
        for x in (self.paths,self.background,self.fraction,self.extinction):
            if not np.isfinite(x).all() or (x<0).any():raise ValueError('Finite nonnegative model arrays required')
            x.flags.writeable=False
        if (self.fraction>1).any() or not np.any(self.fraction>0):raise ValueError('Nonzero blood fractions within[0,1] required')
        if not np.isfinite([hb_g_dl,pulse_fraction]).all() or hb_g_dl<=0 or pulse_fraction<=0:
            raise ValueError('Positive finite hemoglobin and pulse fraction required')
        if (self.fraction*(1+pulse_fraction)>1).any():raise ValueError('Pulse exceeds unit blood fraction')
        self.scale=np.log(10.)*hb_g_dl/64500.;self.pulse_fraction=float(pulse_fraction)
        self.weighted_length=self.paths@self.fraction
        self.background_depth=self.paths@self.background
        self.weighted_length.flags.writeable=False;self.background_depth.flags.writeable=False

    def ratio(self, saturation):
        if not np.isfinite(saturation) or not 0<=saturation<=1:raise ValueError('Saturation outside[0,1]')
        blood=self.scale*(self.extinction[:,0]*saturation+self.extinction[:,1]*(1-saturation))
        base=self.background_depth[:,None]+self.weighted_length[:,None]*blood
        pulse=base+self.pulse_fraction*self.weighted_length[:,None]*blood
        # Shift both states by one common exponent to avoid total underflow.
        shift=base.min(axis=0)
        log0=np.log(np.exp(-(base-shift)).sum(axis=0))
        log1=np.log(np.exp(-(pulse-shift)).sum(axis=0))
        od=log0-log1
        if not np.isfinite(od).all() or np.any(od<=0):raise ValueError('Positive finite pulse signal required')
        return float(od[0]/od[1])

    def solve(self, observed_ratio, bracket=(.5,1.)):
        lo,hi=map(float,bracket)
        if not np.isfinite([observed_ratio,lo,hi]).all() or observed_ratio<=0 or not 0<=lo<hi<=1:
            raise ValueError('Finite positive observation and ordered physical bracket required')
        a=self.ratio(lo)-observed_ratio;b=self.ratio(hi)-observed_ratio
        if a==0:return lo
        if b==0:return hi
        if np.signbit(a)==np.signbit(b):raise ValueError('Observed ratio is not bracketed')
        for _ in range(48):
            mid=(lo+hi)/2;c=self.ratio(mid)-observed_ratio
            if c==0:return mid
            if np.signbit(c)==np.signbit(a):lo=mid;a=c
            else:hi=mid
        return (lo+hi)/2


def selftest():
    paths=np.array([[3.],[3.]])
    extinction=np.array([[586.,1548.52],[1058.,691.32]])
    model=PulsePathInverse(paths,[.004],[.02],extinction)
    saturation=.73
    expected=(extinction[:,0]*saturation+extinction[:,1]*(1-saturation))
    forward_error=abs(model.ratio(saturation)-expected[0]/expected[1])
    inverse_error=abs(model.solve(expected[0]/expected[1])-saturation)
    before=model.ratio(saturation);paths[:]=99
    rejected=0
    for action in (lambda:PulsePathInverse([[-1.]],[.004],[.02],extinction),
                   lambda:PulsePathInverse([[3.]],[.004],[0.],extinction),
                   lambda:model.solve(100.),lambda:model.solve(.8,(1.,.5))):
        try:action()
        except ValueError:rejected+=1
    gates=dict(homogeneous_identity=forward_error<=1e-10,inverse_identity=inverse_error<=1e-10,
               owned_paths=model.ratio(saturation)==before,invalid_refusals=rejected==4)
    return dict(gates={k:bool(v) for k,v in gates.items()},forward_error=float(forward_error),inverse_error=float(inverse_error),rejected=rejected)


if __name__=='__main__':
    import json
    result=selftest();print(json.dumps(result))
    raise SystemExit(0 if all(result['gates'].values()) else 2)
