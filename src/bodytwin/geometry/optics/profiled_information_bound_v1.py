"""Evaluate scalar nuisance-projection bounds from componentwise perturbation radii.

The formulas use floating-point arithmetic, not outward-rounded intervals.
"""
from dataclasses import dataclass
import numpy as np
from bodytwin.geometry.optics.profiled_information_v1 import profile_information


@dataclass(frozen=True)
class ProfileBound:
    nominal: float
    lower: float
    upper: float
    sensitivity_error: float
    nuisance_error: float
    nuisance_sigma_min: float


def profile_bound(interest,nuisance,interest_radius,nuisance_radius):
    s,b,rs,rb=[np.array(x,dtype=float,copy=True) for x in (interest,nuisance,interest_radius,nuisance_radius)]
    if s.ndim!=2 or s.shape[1]!=1 or b.ndim!=2 or b.shape[0]!=s.shape[0] or not b.shape[1] or rs.shape!=s.shape or rb.shape!=b.shape:
        raise ValueError('One interest column and nonempty aligned nuisance/radius matrices required')
    if not all(np.isfinite(x).all() for x in (s,b,rs,rb)) or np.any(rs<0) or np.any(rb<0):
        raise ValueError('Finite arrays and nonnegative radii required')
    baseline=profile_information(s,b)
    if baseline.nuisance_rank!=b.shape[1]:raise ValueError('Full-column-rank nuisance required')
    sigma=float(np.linalg.svd(b,compute_uv=False)[-1]);es=float(np.linalg.norm(rs));eb=float(np.linalg.norm(rb))
    if eb>=sigma:raise ValueError('Perturbation may destroy nuisance rank')
    difference=es+min(1.,eb/(sigma-eb))*float(np.linalg.norm(s))
    nominal=float(np.linalg.norm(baseline.projected))
    return ProfileBound(nominal,max(0.,nominal-difference),nominal+difference,difference,eb,sigma)


def selftest():
    import itertools
    s=np.array([[1.],[2.],[3.]]);b=np.array([[1.],[0.],[0.]])
    rs=np.full_like(s,.01);rb=np.full_like(b,.02)
    copies=[x.copy() for x in (s,b,rs,rb)]
    first=profile_bound(s,b,rs,rb);second=profile_bound(s,b,rs,rb)
    contained=[]
    for signs in itertools.product((-1.,1.),repeat=6):
        a=s+rs*np.array(signs[:3])[:,None];n=b+rb*np.array(signs[3:])[:,None]
        # Independent rank-one formula, without the production SVD projection.
        score=float(np.linalg.norm(a-n*float((n.T@a).item())/float((n.T@n).item())))
        contained.append(first.lower<=score<=first.upper)
    zero=profile_bound(s,b,np.zeros_like(s),np.zeros_like(b))
    refused=0
    for args in ((s,b,-rs,rb),(s,b,rs,np.ones_like(rb)),(s,np.column_stack((b,b)),rs,np.zeros((3,2))),
                 (s.ravel(),b,rs,rb),(s,b,rs[:2],rb),(s*np.nan,b,rs,rb)):
        try:profile_bound(*args)
        except ValueError:refused+=1
    return dict(zero_radius=zero.lower==zero.upper==float(np.sqrt(13.)),corner_enclosure=all(contained) and len(contained)==64,
                repeat=first==second,nonmutation=all(np.array_equal(a,b) for a,b in zip(copies,(s,b,rs,rb))),refusals=refused==6)


if __name__=='__main__':
    import json
    gates=selftest();print(json.dumps(dict(gates=gates)))
    raise SystemExit(0 if all(gates.values()) else 2)
