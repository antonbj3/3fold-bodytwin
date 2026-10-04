#!/usr/bin/env python3
"""Retain raw normalized-state overshoot and refine the inherited HX RHS."""
import json
import time
import numpy as np
from scipy.integrate import solve_ivp
from incision_r1 import LANE,ROOT,parent_module,write_json


def main():
    raw=json.loads((LANE/"r1/physical_patch/HX_chain.json").read_text())
    q36=parent_module("q36_refinement",ROOT/"results/BT-HX-Q036/model.py")
    p=q36.Parameters()
    t0=time.process_time()
    rows=[]
    bounded=np.array([0,1,2,3,8,9,10,11])
    for seed in (0.,raw["mechanical_initial_state_port"]["value"]):
        sc=q36.Scenario(name="edge_refined_synthetic",pre_ischemia_min=0,ischemia_min=10,reperfusion_min=10,q_ischemia=max(.01,1-seed))
        state=q36.healthy_state(p)
        state[8]=seed*.2;state[9]=seed*.1
        snapshots=[];max_overshoot=0.
        for lo,hi in ((0,10),(10,20)):
            solve=solve_ivp(lambda t,y:q36.rhs(t,y,p,sc),(lo,hi),state,method="DOP853",rtol=1e-10,atol=1e-12,max_step=.05)
            assert solve.success
            state=solve.y[:,-1]
            max_overshoot=max(max_overshoot,float(max(0,np.max(solve.y[bounded]-1),np.max(-solve.y[bounded]))))
            raw_state=state.copy()
            # Same normalized snapshot projection as inherited simulate().
            state[bounded]=np.clip(state[bounded],0,1)
            state[6]=np.clip(state[6],0,5)
            state[[4,5,7]]=np.maximum(state[[4,5,7]],0)
            snapshots.append({"time_min":hi,"raw_state":q36.unpack(raw_state),"reported_projected_state":q36.unpack(state)})
        rows.append({"mechanical_seed":seed,"max_raw_fraction_bound_violation":max_overshoot,"snapshots":snapshots})
    write_json(LANE/"r1/HX_CHAIN_BOUND_REFINEMENT.json",{"original":"r1/physical_patch/HX_chain.json","reason":"Raw RHS integration had ATP overshoot2.4989e-7 at20min; retain it and refine, then apply inherited snapshot projection explicitly","rtol":1e-10,"atol":1e-12,"max_step_min":.05,"cases":rows,"cpu_s":time.process_time()-t0,"review":"PENDING_INDEPENDENT_REVIEW"})


if __name__=="__main__":
    main()
