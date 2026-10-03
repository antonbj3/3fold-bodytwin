"""Dimensioned conservative oxygen inventory, testing the quasi-steady reduction."""
from dataclasses import replace,asdict
import hashlib
import json
import time
import numpy as np
from . import response_r1 as r

def early_storage(g,p,alpha,face_history=None):
    start=time.perf_counter()
    solver=r.Oxygen(g,p)
    Y=np.repeat(r.q036.healthy_state(r.QP)[:,None],g.n,axis=1)
    P=np.full(g.n,p.healthy_Torr)
    H=np.zeros(g.n);cf=r.cut_flow(g,p)
    r.QNS['external_flow']=lambda t,sc:cf
    cumulative_net=0.;cumulative_supply=0.;cumulative_use=0.
    initial_inventory=float(alpha*np.sum(g.dx*P))
    series=[];maxproject=0.
    def kinetics(Y,P):
        Ys=Y.copy();Ys[1]=np.clip(P/r.O_SCALE,0,1)
        dy=r.QRHS(0,Ys,r.QP,r.SCENARIO);dy[1]=0.
        h=r.QP.k_injury*Ys[8]*(1-Ys[2])*(1+np.maximum(0,Ys[6]-1))
        return dy,h
    steps=int(np.ceil(p.early_minutes/p.early_dt_min));dt=p.early_minutes/steps
    for i in range(steps):
        dy,h=kinetics(Y,P)
        mid=r.project_q(Y+.5*dt*dy)
        demand=p.consumption*(.35+.65*mid[2]*(1-mid[8]))
        if face_history is not None:solver.p=replace(p,face_conductance=float(face_history((i+.5)*dt)))
        Pmid=solver.solve(mid[0],demand,storage_alpha=alpha,dt_min=.5*dt,previous_pressure=P)
        dm,hm=kinetics(mid,Pmid)
        raw=Y+dt*dm
        Yn=r.project_q(raw)
        maxproject=max(maxproject,float(np.max(np.abs(Yn-raw))))
        demand_end=p.consumption*(.35+.65*Yn[2]*(1-Yn[8]))
        if face_history is not None:solver.p=replace(p,face_conductance=float(face_history((i+1)*dt)))
        Pn=solver.solve(Yn[0],demand_end,storage_alpha=alpha,dt_min=dt,previous_pressure=P)
        cumulative_net+=dt*(solver.last_supply-solver.last_consumption)
        cumulative_supply+=dt*solver.last_supply;cumulative_use+=dt*solver.last_consumption
        Y=Yn;P=Pn;Y[1]=np.clip(P/r.O_SCALE,0,1);H+=dt*np.maximum(hm,0)
        if i % max(1,round(1/dt))==0 or i==steps-1:
            series.append([(i+1)*dt,float(P[0]),float(Y[2,0]),float(Y[6,0]),float(Y[8,0]),
                           float(H[0]),float(alpha*np.sum(g.dx*P)),cumulative_net])
    final_inventory=float(alpha*np.sum(g.dx*P))
    balance=abs(final_inventory-initial_inventory-cumulative_net)/max(initial_inventory,abs(cumulative_supply),abs(cumulative_use),1e-12)
    return {'state':Y,'hazard':H,'pressure':P,'snapshots':[],
            'wall_s':time.perf_counter()-start,'oxygen':solver.stats(),'projection':maxproject,
            'inventory':{'initial_ml_m2':initial_inventory,'final_ml_m2':final_inventory,
                         'integrated_net_supply_ml_m2':cumulative_net,'cumulative_inflow_ml_m2':cumulative_supply,
                         'cumulative_consumption_ml_m2':cumulative_use,'relative_conservation_error':float(balance),
                         'terminal_face_flux_ml_m2_min':solver.last_left_flux,
                         'alpha_ml_m3_Torr':alpha,'alpha_status':'SYNTHETIC_SENSITIVITY'},'series':np.array(series)}

def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--out',default=str(r.ROOT/'r1'/'storage'))
    out=r.Path(parser.parse_args().out).resolve();out.mkdir(parents=True,exist_ok=False)
    pre=r.ROOT/'PREREG_R1_STORAGE.json'
    if hashlib.sha256(pre.read_bytes()).hexdigest()!=(r.ROOT/'PREREG_R1_STORAGE.sha256').read_text().split()[0]:raise RuntimeError('PREREG changed')
    rates=json.loads((r.ROOT/'r1'/'coupled_v2'/'calibration.json').read_text())['rates_day_inverse']
    p=r.Params(width_m=.00075,k_deposit_day=rates[0],k_mature_day=rates[1]);g=r.Grid(p)
    start=time.perf_counter();cases=[];stored={}
    for alpha in [10.,30.,100.]:
        e=early_storage(g,p,alpha);l=r.long(g,p,e)
        name=f'alpha{int(alpha)}_dt0p1';r.save_case(out,name,g,p,e,l)
        np.savez_compressed(out/(name+'_inventory.npz'),series=e['series'])
        cases.append({'alpha':alpha,'dt_min':.1,'inventory':e['inventory'],'observables':r.observables(g,p,e,l)})
        stored[alpha]=(e,l)
        print(name,e['inventory'],flush=True)
    refinements=[]
    prev=stored[30.]
    for dt in [.05,.025]:
        pp=replace(p,early_dt_min=dt);e=early_storage(g,pp,30.);l=r.long(g,pp,e)
        name='alpha30_dt'+str(dt).replace('.','p');r.save_case(out,name,g,pp,e,l)
        refinements.append({'dt_min':dt,'previous_dt':.1 if dt==.05 else .05,
                            'difference':r.compare(g,e,l,g,*prev),'inventory':e['inventory']})
        prev=(e,l)
    # Same coupled transient equations on uniform fine spatial mesh.
    gf=r.Grid(p,'uniform');ef=early_storage(gf,p,30.);lf=r.long(gf,p,ef)
    r.save_case(out,'alpha30_uniform',gf,p,ef,lf)
    spatial=r.compare(gf,ef,lf,g,*stored[30.])
    raw=np.load(r.ROOT/'r1'/'coupled_v2'/'w750_sealed_graded.npz')
    quasi={'pressure':raw['early_pressure'],'hazard':raw['hazard'],'state':raw['early_state']}
    ql={'strength_pct':raw['strength_pct'],'days':raw['days']}
    ablations=[]
    for alpha,(e,l) in stored.items():
        ablations.append({'alpha':alpha,'quasi_steady_difference':r.compare(g,e,l,g,quasi,ql),
                          'injury_proxy_width_difference_m':float(abs(np.sum(g.dx*(1-np.exp(-e['hazard'])>.2))-np.sum(g.dx*(1-np.exp(-quasi['hazard'])>.2))))})
    r.write_json(out/'summary.json',{'parameters':asdict(p),'cases':cases,'time_refinements':refinements,
        'spatial_control':spatial,'quasi_steady_ablations':ablations,
        'strongest_control':{'method':'Conventional implicit FV storage + identical Q036 state','gate':'TIE','gain':1.},
        'full_storage_experiment_s':time.perf_counter()-start,'review_state':'PENDING_INDEPENDENT_REVIEW',
        'alpha_empirical_calibration':'UNKNOWN','empirical_joint_prediction':'UNKNOWN'})

if __name__=='__main__':main()
