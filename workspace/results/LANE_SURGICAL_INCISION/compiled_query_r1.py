#!/usr/bin/env python3
"""NumPy-only consumer of the centered residual capsule.

Rejected queries require the executed full fallback in centered_residual_r1.py.
They are not returned as accepted predictions.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np


def evaluate(op,length_mm,kappa_multiplier=1.):
    t0=time.process_time()
    PB,RC=op["PB"],op["RC"]
    x=op["x"]
    cut=(x>=.008)&(x<.008+length_mm*1e-3-1e-12)
    diag=op["k0"]*np.exp(-op["kappa"]*kappa_multiplier/op["delta0"])*(~cut)
    K=PB.T @ (diag[:,None]*PB)
    offset=op["offset_PB"]
    force=PB.T @ (diag*offset)
    combos=op["combos"]
    def mono(z):
        vals=np.r_[1.,z][combos]
        v=np.prod(vals,axis=1)
        dv=np.zeros((len(combos),len(z)))
        for j in range(3):
            other=[k for k in range(3) if k!=j]
            valid=combos[:,j]>0
            np.add.at(dv,(np.flatnonzero(valid),combos[valid,j]-1),vals[valid,other[0]]*vals[valid,other[1]])
        return v,dv
    z=op["zref"].copy()
    for _ in range(50):
        v,dv=mono(z)
        r=RC @ v+K @ z+force
        if np.linalg.norm(r)/op["reduced_scale"]<1e-10:
            break
        step=np.linalg.solve(RC @ dv+K,-r)
        alpha=1.
        for _ in range(20):
            nv,_=mono(z+alpha*step)
            nr=RC @ nv+K @ (z+alpha*step)+force
            if np.linalg.norm(nr)<np.linalg.norm(r):
                z+=alpha*step
                break
            alpha*=.5
        else:
            raise RuntimeError("Newton line search failed")
    else:
        raise RuntimeError("Newton failed")
    s=diag*(PB @ z+offset)
    norm2=float(v @ op["H"] @ v+2*s @ op["D"] @ v+s @ op["G"] @ s)
    absv,abss=np.abs(v),np.abs(s)
    allowance=128*np.finfo(float).eps*float(absv @ op["H_abs"] @ absv+2*abss @ op["D_abs"] @ absv+abss @ op["G_abs"] @ abss)
    bound=float(op["compliance"]*np.sqrt(max(0,norm2)+allowance))
    return {"status":"CONDITIONAL_DISCRETE_ACCEPT" if bound<=1e-7 else "REQUIRES_FULL_FALLBACK",
        "gap_m":(2*(PB @ z+offset)).tolist() if bound<=1e-7 else None,
        "gap_bound_m":bound,"cpu_s":time.process_time()-t0,
        "global_solves":0,"global_state_reconstructions":0,
        "validity":"supplied frozen cohesive history; synthetic scalar discrete convex model; engineering arithmetic allowance, pending review",
        "review":"PENDING_INDEPENDENT_REVIEW"}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--length-mm",type=float,default=9.)
    parser.add_argument("--kappa-multiplier",type=float,default=1.)
    args=parser.parse_args()
    path=Path(__file__).resolve().parent/"r1/centered_residual/operator_tier2.npz"
    with np.load(path) as f:
        op={k:f[k] for k in f.files}
    result=evaluate(op,args.length_mm,args.kappa_multiplier)
    # Keep stdout compact; full gap is present in Python API, no field reconstruction.
    if result["gap_m"] is not None:
        result["gap_max_m"]=max(result.pop("gap_m"))
    print(json.dumps(result,allow_nan=False))


if __name__=="__main__":
    main()
