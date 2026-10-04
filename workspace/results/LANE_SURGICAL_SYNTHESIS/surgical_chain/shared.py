"""Single collagen reference and synchronized FV/mark fluxes (R3).

C is dimensionless collagen mass per fixed registered reference volume. Mark
traces partition exactly that C. Proportional removal and isotropic precursor
birth are synthetic closures, not a measured species-specific removal law.
"""
import time
import numpy as np
from .vendor import response_r1 as r
from .inventories import query_strength

VERSION = 'shared_reference_r3'

def gate(P, p):
    if not p.oxygen_feedback:
        return np.ones_like(P)
    return np.minimum(P/(10+P)/(p.healthy_Torr/(10+p.healthy_Torr)), 1.)

def collagen_ports(Y, h, p, remove=True):
    return p.k_deposit_day*Y[2]*h*(1-Y[3]), (.008*Y[1] if remove else np.zeros_like(h))

def mark_rhs(Q, birth, loss_rate, h, k1, k2):
    u, i, m = Q
    source = np.zeros_like(u)
    source[0,0] = source[1,1] = .5*birth
    j1, j2 = k1*h*u, k2*h*i
    return np.array([source-j1-loss_rate*u, j1-j2-loss_rate*i, j2-loss_rate*m])

def cohort_rhs(z, birth, loss_rate, h, k1, k2):
    # Independent conventional six diagonal direction cohorts, same information.
    u, i, m = z
    j1, j2 = k1*h*u, k2*h*i
    return np.array([.5*birth-j1-loss_rate*u, j1-j2-loss_rate*i, j2-loss_rate*m])

def chemical_rhs(z, birth, loss_rate, h, params, sites):
    # D -> 2H with yield eta or A with yield 1-eta. Site-equivalent ledger
    # is D+2H+A = sites*C; sites/mol and collagen-reference mass remain distinct.
    d, hp, a = z
    conversion = params['chemistry_rate']*d
    eta = params['chemistry_yield']
    return np.array([sites*birth-conversion-loss_rate*d,
                     .5*eta*conversion-loss_rate*hp,
                     (1-eta)*conversion-loss_rate*a])

def late_response(g, p, e, params, until, dt, *, checkpoint=None,
                  removal=True, compare_cohorts=False):
    if not np.isclose(dt,p.dt_day):
        raise ValueError('FV and inventory dt must match')
    if checkpoint is not None and str(checkpoint.get('reference_version','')) != VERSION:
        raise ValueError('Legacy collagen reference requires explicit conversion or day0 rebuild')
    st = time.perf_counter()
    oxygen = r.Oxygen(g,p)
    H = e['hazard'] if p.memory_transfer else np.zeros(g.n)
    if checkpoint is None:
        Y = np.zeros((6,g.n)); Y[5] = np.clip(e['state'][0],.01,1)
        Q = np.zeros((3,2,2,g.n)); chem = np.zeros((3,g.n))
        born = np.zeros(g.n); removed = np.zeros(g.n); initial = np.zeros(g.n)
        start = 0.
    else:
        Y = checkpoint['long_state'].copy(); Q = checkpoint['Q'].copy()
        chem = checkpoint['chemical_marks'].copy()
        born = checkpoint['born'].copy(); removed = checkpoint['removed'].copy()
        initial = checkpoint['initial'].copy(); start = float(checkpoint['absolute_day'])
    z = np.stack([Q[:,0,0],Q[:,1,1]],axis=1)
    if compare_cohorts and np.max(abs(Q[:,0,1]))>1e-12:
        raise ValueError('Diagonal cohort control only for diagonal initial tensors')
    sites = params['chemistry_D10']+2*params['chemistry_H10']+params['chemistry_A10']
    if checkpoint is not None:
        mass=np.trace(Q.sum(0),axis1=0,axis2=1)
        if (np.max(abs(mass-Y[3]))>1e-10 or
            np.max(abs(mass+removed-initial-born))>1e-10 or
            np.max(abs(chem[0]+2*chem[1]+chem[2]-sites*Y[3]))>1e-10):
            raise ValueError('Checkpoint reference or birth/removal ledger mismatch')
    steps = round((until-start)/dt)
    if steps<=0 or not np.isclose(steps*dt,until-start):
        raise ValueError('Positive integer-step suffix required')
    monitor = np.maximum(0,np.minimum(g.edges[1:],p.width_m)-g.edges[:-1])/p.width_m
    days=[]; strengths=[]; masses=[]; Cseries=[]; chemseries=[]; ledger=[]
    fields=[]; pressures=[]; field_days=[]; legacy_strength=[]; saved={}; max_balance=0.; max_ref=0.
    max_sites=0.; min_mark=0.; max_control=0.; max_projection=0.
    def evaluate(t, state, q, ch, cz):
        N,M,F,C,X,V = state
        flow = V*p.outside_flow/(1+.20*M+.50*(1-np.exp(-H))*np.exp(-t/7))
        pressure = oxygen.solve(flow,p.consumption*(1+.12*N+.20*M+.08*F),C)
        h = gate(pressure,p)
        birth, loss = collagen_ports(state,h,p,removal)
        dy = r.long_rhs(t,state,pressure,g,p,H)
        dy[3] = birth-loss*C
        return dy,mark_rhs(q,birth,loss,h,params['k_U_I'],params['k_I_M']),chemical_rhs(ch,birth,loss,h,params,sites),cohort_rhs(cz,birth,loss,h,params['k_U_I'],params['k_I_M']),birth,loss,pressure
    for j in range(steps):
        t=start+j*dt
        Y[2] = r.diffuse_fibro(Y[2],g,p.fibro_diffusion_m2_day,.5*dt)
        dy,dq,dc,dz,_,_,_ = evaluate(t,Y,Q,chem,z)
        mid=Y+.5*dt*dy
        # Keep vendor state guards, but C is never projected separately from Q.
        mid[:2]=np.maximum(mid[:2],0); mid[[2,4,5]]=np.clip(mid[[2,4,5]],0,1)
        if np.any(mid[3]<-1e-12) or np.any(mid[3]>1+1e-12):
            raise RuntimeError('Collagen midpoint outside bounds; refine dt')
        dy,dq,dc,dz,b,loss,P = evaluate(t+.5*dt,mid,Q+.5*dt*dq,chem+.5*dt*dc,z+.5*dt*dz)
        raw=Y+dt*dy; Y=raw.copy(); Y[:2]=np.maximum(Y[:2],0)
        Y[[2,4,5]]=np.clip(Y[[2,4,5]],0,1)
        max_projection=max(max_projection,float(abs(raw-Y).max()))
        Q+=dt*dq; chem+=dt*dc; z+=dt*dz
        born+=dt*b; removed+=dt*loss*mid[3]
        Y[2]=r.diffuse_fibro(Y[2],g,p.fibro_diffusion_m2_day,.5*dt)
        mass=np.trace(Q.sum(0),axis1=0,axis2=1)
        ref=float(abs(mass-Y[3]).max()); bal=float(abs(mass+removed-initial-born).max())
        site_err=float(abs(chem[0]+2*chem[1]+chem[2]-sites*Y[3]).max())
        max_ref=max(max_ref,ref); max_balance=max(max_balance,bal); max_sites=max(max_sites,site_err)
        mineig=float(np.linalg.eigvalsh(np.moveaxis(Q,-1,0)).min())
        min_mark=min(min_mark,mineig,float(chem.min()),float(Y[3].min()))
        if min_mark<-1e-10 or np.any(Y[3]>1+1e-10):
            raise RuntimeError('Collagen marks out of bounds; refine dt')
        if compare_cohorts:
            max_control=max(max_control,float(abs(z[:,0]-Q[:,0,0]).max()),float(abs(z[:,1]-Q[:,1,1]).max()))
        tt=start+(j+1)*dt; days.append(tt)
        strengths.append(float(query_strength(Q[2],params['load_angle'],params['bridge_fraction'])@monitor))
        legacy_strength.append(float(75*(Y[3]*Y[4])@monitor))
        masses.append([float(np.trace(a,axis1=0,axis2=1)@monitor) for a in Q])
        Cseries.append(float(Y[3]@monitor)); chemseries.append(chem@monitor)
        ledger.append([float(born@monitor),float(removed@monitor)])
        if any(abs(tt-d)<dt/3 for d in [1,3,7,10,14,21,28,42,90]):
            fields.append(Y.copy()); pressures.append(P.copy()); field_days.append(tt)
        if abs(tt-21)<dt/3:
            saved=dict(absolute_day=tt,Q=Q.copy(),born=born.copy(),removed=removed.copy(),initial=initial.copy(),long_state=Y.copy(),chemical_marks=chem.copy(),reference_version=np.array(VERSION))
    ch=np.array(chemseries); cs=np.array(Cseries)
    # Explicit diagnostic transfer only; no native chemistry->rupture law.
    hp=np.divide(ch[:,1],cs,out=np.zeros_like(cs),where=cs>0)
    fraction=np.clip((hp-.003)/.040,0,1)
    chemical_strength=75*params['bridge_fraction']*cs*fraction
    legacy=dict(days=np.array(days),strength_pct=np.array(legacy_strength),state=Y,pressure=oxygen.P.copy(),saved_days=np.array(field_days),fields=np.array(fields),pressures=np.array(pressures),wall_s=time.perf_counter()-st,oxygen=oxygen.stats(),projection=max_projection)
    return dict(days=np.array(days),strength_pct=np.array(strengths),UIM=np.array(masses),C=np.array(Cseries),chemical_marks_series=ch,chemical_strength_pct=chemical_strength,ledger=np.array(ledger),Q_final=Q,born=born,removed=removed,initial=initial,checkpoint=saved,legacy=legacy,max_mass_balance=max_balance,
                shared_reference=dict(version=VERSION,reference='Fixed registered reference volume; C and tensor traces normalized by same intact collagen mass',formation='k_deposit F h(1-C)',removal='0.008 macrophage C; proportional across all marks',additional_UI_death=0.,unused_legacy_ports=['birth_scale','turnover_U_I'],max_pointwise_trace_error=max_ref,max_collagen_ledger_error=max_balance,max_site_ledger_error=max_sites,min_mark_eigenvalue=min_mark,max_cohort_tensor_error=max_control if compare_cohorts else None,site_equivalents_per_collagen_reference=sites,chemical_maturity_upper_bound=max(0,(.5*params['chemistry_yield']*sites-.003)/.040),synthetic=True,native_rupture_law=None))

def run(config, **kwargs):
    from .chain import run as legacy_run
    return legacy_run(config, inventory='shared', **kwargs)
