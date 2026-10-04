"""Data-constrained knee load, feasible contact set, and measurement gap.

The law band is a literature-parameter sensitivity envelope, not a probability
interval. LP bounds concern affine compressive forces, not vector norms.
"""
from __future__ import annotations

from collections.abc import Mapping
import numpy as np
from .solver import solve_bounds

RQ_GRID = np.array([.030, .040, .050, .060])
RA_GRID = np.array([.035, .045, .060])
BETA_GRID = np.array([.20, .29, .40])
SENSOR_SET = ("patellar tendon force", "semimembranosus tendon force",
              "biceps femoris short-head tendon force")


def _array(trial, key, n=None):
    x = np.asarray(trial[key], float)
    if n is not None and len(x) != n:
        raise ValueError(f"{key}: frame count mismatch")
    return x


def knee_load(trial: Mapping, *, system: Mapping | None = None,
              bounds: Mapping | None = None):
    """Return aligned per-frame load arrays and explicit evidence status.

    Required trial keys: grf (T,3), knee_moment_Nm, ankle_moment_Nm,
    knee_flex_deg. Optional measured_N, n1g_N, bw_N and any selected_N.
    ``system`` uses the L1 contact map (A,b,F0,cj,C0,Mkx,mx,d), with
    optional per-frame cap. ``bounds`` is a previously solved, frame-aligned
    result from the *same* system and constraints; it avoids repeat LP work.
    """
    grf = _array(trial, 'grf')
    if grf.ndim != 2 or grf.shape[1] != 3:
        raise ValueError('grf must have shape (frames,3)')
    n = len(grf)
    mk = np.abs(_array(trial, 'knee_moment_Nm', n))
    ma = np.abs(_array(trial, 'ankle_moment_Nm', n))
    angle = _array(trial, 'knee_flex_deg', n)
    g = np.linalg.norm(grf, axis=1)
    rq = .040 + .010 * np.maximum(0, 1 - np.abs(angle - 45)/45)
    nominal = g + mk/rq + .29*ma/.045
    low = g + mk/RQ_GRID.max() + BETA_GRID.min()*ma/RA_GRID.max()
    high = g + mk/RQ_GRID.min() + BETA_GRID.max()*ma/RA_GRID.min()
    if not (np.isfinite(nominal).all() and np.isfinite(low).all() and np.isfinite(high).all()):
        raise ValueError('nonfinite ID inputs')
    if system is not None and bounds is not None:
        raise ValueError('provide system or precomputed bounds')
    if system is not None:
        bounds = joint_load(system)
    if bounds is not None:
        intervals = {k: np.asarray(v, float) for k,v in bounds.items()}
        if any(v.shape != (n,2) for v in intervals.values()):
            raise ValueError('bounds must be frame-aligned (frames,2) arrays')
        state = 'feasible_set_available'
    else:
        intervals = None
        state = 'no_verified_equilibrium_contact_map'
    measured = np.asarray(trial['measured_N'],float) if 'measured_N' in trial else None
    if measured is not None and measured.shape != (n,):
        raise ValueError('measured_N: frame count mismatch')
    data = {'grf_N':g, 'knee_moment_Nm':_array(trial,'knee_moment_Nm',n),
            'ankle_moment_Nm':_array(trial,'ankle_moment_Nm',n),
            'law_N':nominal, 'law_band_N':np.stack([low,high],axis=1),
            'set_N':intervals, 'set_status':state,
            'measured_N':measured, 'n1g_N':np.asarray(trial['n1g_N'],float) if 'n1g_N' in trial else None,
            'selected_N':np.asarray(trial['selected_N'],float) if 'selected_N' in trial else None,
            'measurement':measurement_plan(),
            'interpretation':'ID determines external load and net moments; contact law remains lever-arm dependent, and muscle sharing requires the feasible set.'}
    if intervals is not None and 'total' in intervals:
        lo,hi=intervals['total'].T
        data['set_feasible']=np.isfinite(lo)&(hi>=lo)
    if measured is not None:
        data['law_contains_measured'] = (low-1e-6 <= measured)&(measured <= high+1e-6)
        if intervals is not None and 'total' in intervals:
            lo,hi=intervals['total'].T
            data['set_contains_measured']=(lo-1e-5<=measured)&(measured<=hi+1e-5)
            data['measured_set_position']=np.divide(measured-lo,hi-lo,
                out=np.full(n,np.nan),where=np.isfinite(hi-lo)&(hi>lo))
    if intervals is not None and data['selected_N'] is not None and 'total' in intervals:
        data['selected_inside_set']=(intervals['total'][:,0]-1e-5<=data['selected_N'])&(data['selected_N']<=intervals['total'][:,1]+1e-5)
    if intervals is not None:
        for compartment in ('medial','lateral'):
            key='measured_'+compartment+'_N'
            if key in trial and compartment in intervals:
                value=_array(trial,key,n)
                lo,hi=intervals[compartment].T
                data[key]=value
                data[compartment+'_contains_measured']=(lo-1e-5<=value)&(value<=hi+1e-5)
    return data


def joint_load(system: Mapping, *, on_infeasible='mark'):
    """Exact LP limits for affine L1 total, medial and lateral compression.

    Zero moment reserves, f >= 0, f <= Fmax (and optional FLFV cap), and
    nonnegative medial/lateral compression are imposed on every frame.
    Infeasible frames receive NaN bounds by default; ``on_infeasible='raise'``
    preserves solve_bounds' strict behavior. Other solver errors still raise.
    """
    if on_infeasible not in ('mark','raise'):
        raise ValueError("on_infeasible must be 'mark' or 'raise'")
    cj=np.asarray(system['cj'],float)
    if cj.ndim != 2: raise ValueError('cj must be (frames,muscles)')
    t,n=cj.shape
    mkx=np.asarray(system['Mkx'],float)
    c0=np.asarray(system['C0'],float)
    mx=np.asarray(system['mx'],float)
    d=float(system['d'])
    if d<=0: raise ValueError('contact spacing d must be positive')
    maps={'total':cj[:,None,:],
          'medial':(.5*cj+mkx/d)[:,None,:],
          'lateral':(.5*cj-mkx/d)[:,None,:]}
    offsets={'total':c0[:,None], 'medial':(.5*c0-mx/d)[:,None],
             'lateral':(.5*c0+mx/d)[:,None]}
    g=np.stack([-maps['medial'][:,0,:],-maps['lateral'][:,0,:]],axis=1)
    h=np.stack([offsets['medial'][:,0],offsets['lateral'][:,0]],axis=1)
    limit=np.broadcast_to(np.asarray(system['F0'],float),(t,n)).copy()
    if 'cap' in system: limit=np.minimum(limit,np.asarray(system['cap'],float))
    a=np.asarray(system['A'],float);b=np.asarray(system['b'],float)
    if a.shape!=(t,b.shape[1],n) or b.shape[0]!=t:
        raise ValueError('A/b shape mismatch')
    result={key:np.full((t,2),np.nan) for key in maps}
    for i in range(t):
        try:
            raw=solve_bounds(a[i],b[i],{key:value[i] for key,value in maps.items()},
                {key:value[i] for key,value in offsets.items()},upper=limit[i],
                A_ub=g[i],b_ub=h[i])
        except ValueError as exc:
            if on_infeasible=='raise' or 'infeasible frame' not in str(exc):raise
            continue
        for key,value in raw.items():result[key][i]=value[0,0]
    return result


def measurement_plan():
    return {'minimal_local_sensor_set':SENSOR_SET,
            'conditional_base':'moment + six stiffness rows + eleven EMG ratios in gait2392',
            'expected_95pct_width_N':{
                'stance':{'cost_4_moment':[1959,1611], 'cost_20_base':[590,848], 'cost_46_triple':[452,403]},
                'swing':{'cost_20_base':[539,708], 'cost_46_triple':[408,365]}},
            'units':'[medial,lateral] N; relative cost units',
            'status':'local gait2392 operator and assumed prior/noise; no GC or the collaborator closure demonstrated',
            'current_data_enough':False,
            'plain_text':'Current GRF, moments and surface EMG do not determine medial/lateral contact; separate calibrated tendon-force measurements and a verified person-specific contact map are needed.'}
