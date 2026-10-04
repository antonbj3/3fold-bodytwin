"""A369 fixed-support implicit KKT derivative, copied unchanged to vendor/whatif/implicit_kkt.py."""
from .vendor.whatif.implicit_kkt import frame
from .config import paths
import numpy as np

def _read(prefix, tag, t):
    z=np.load(prefix/f'sys_{tag}.npz')
    l=np.load(prefix/f'leg_sys_{tag}.npz')
    return ({k:z[k][t:t+1] if k in ('A','b','m') else z[k][:,t:t+1] if k=='N' else z[k] for k in ('A','b','m','N','names')},
            {k:l[k][t:t+1] for k in ('Whip','dhip')})

def derivative(tag, parameter, t, baseline_force, config=None):
    p=paths(config); base=p['solver_data']; mat=p['whatif_matrices']
    z0,l0=_read(base,tag,t)
    pairs=[]
    for delta in (-1,1):
        stem=f'{tag}_{parameter}_{"m" if delta<0 else "p"}1'
        z,l=_read(mat/'data',stem,t)
        from .solver import solve_step
        m=int(z['m'][0]);f1,receipt=solve_step(z['A'][0,:m],z['b'][0,:m],z['N'][:,0],z['names'])
        if receipt['kkt']>1e-10:raise RuntimeError(f'perturbed {stem} KKT {receipt["kkt"]}')
        f=f1[None]
        pairs.append((z,l,f))
    zm,lm,fm=pairs[0];zp,lp,fp=pairs[1]
    return frame(0,z0,l0,baseline_force[None],zm,lm,fm,zp,lp,fp)

def reported_curve(tag, parameter, config=None):
    import json
    p=paths(config)
    return json.loads((p['whatif_results']/'results.json').read_text())['curves'][tag][parameter]

def reported_derivative(tag, parameter, preferred_frame=None, config=None):
    """Select a stable frame from A369's archived 141-frame implicit calculation."""
    import json
    p=paths(config)
    row=json.loads((p['whatif_results']/'implicit_kkt_results.json').read_text())[tag][parameter]
    if 'rows' not in row:
        return {'status':row.get('status','unavailable'), 'source':'A369 archived implicit KKT'}
    good=[(i,r) for i,r in enumerate(row['rows']) if 'implicit_N_per_deg' in r]
    if not good:return {'status':'no stable active-set frame','source':'A369 archived implicit KKT'}
    chosen=min(good,key=lambda pair:abs(pair[0]-(preferred_frame if preferred_frame is not None else row['peak_frame'])))
    return {'frame':chosen[0],**chosen[1],'source':'A369 archived implicit KKT; fixed active set'}

def field_u380_shape_massprop(mesh, u, source_dir=None):
    """Call Field U380 at its source path; see CX-D1PARITY/U380_SOURCE.sha256.

    mesh uses V in mm, F triangle indices and an edit frame; u is
    (twist_rad, varus_rad, lengthening_mm). No Field source is vendored.
    """
    import importlib.util, sys
    from pathlib import Path
    src=Path(source_dir or '../3fold-motion-engine/_private/romi_collab/build/U380/code')
    file=src/'shape_massprop.py'
    if not file.is_file():raise FileNotFoundError(file)
    old=list(sys.path)
    try:
        sys.path.insert(0,str(src))
        spec=importlib.util.spec_from_file_location('bodytwin_field_u380',file)
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        return mod.shape_massprop(mesh,u)
    finally:
        sys.path[:]=old
