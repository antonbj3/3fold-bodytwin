"""Ordered, data-defined synthetic NIRS graphs over a fixed pure-function registry."""
import copy
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from bodytwin.geometry.optics.nirs_forward_v1 import signal,sao2_hill,BLOOD_FRACTION


def _record(po2,hb,fraction,sat,optical):
    return [po2,hb,fraction,sat,*np.array(optical['mua_per_mm']).ravel(),*optical['intensities'],optical['ratio_760_850']]


REGISTRY={
    'hill': (('po2',), lambda po2: float(sao2_hill(po2))),
    'blood_fraction': (('scale',), lambda scale: BLOOD_FRACTION*scale),
    'signal': (('paths','sat','hb','fraction'), lambda paths,sat,hb,fraction: signal(paths,sat,hb,fraction)[0]),
    'record': (('po2','hb','fraction','sat','optical'), _record),
}


INPUT_TYPES={'hill':{'po2':'scalar'},'blood_fraction':{'scale':'scalar'},
             'signal':{'paths':'paths','sat':'scalar','hb':'scalar','fraction':'fractions'},
             'record':{'po2':'scalar','hb':'scalar','fraction':'scalar','sat':'scalar','optical':'optical'}}
OUTPUT_TYPES={'hill':'scalar','blood_fraction':'fractions','signal':'optical','record':'record'}


class Chain:
    """Validate a JSON-compatible graph; run with explicitly supplied path arrays."""
    def __init__(self,spec):
        spec=copy.deepcopy(spec)
        required={'version','fixture','nominal','bands','seed','draws','nodes','output'}
        if not isinstance(spec,dict) or set(spec)!=required or spec['version']!=1 or spec['fixture'] not in ('box','curved'):
            raise ValueError('Unsupported chain schema')
        nominal=np.asarray(spec['nominal'],dtype=float);bands=np.asarray(spec['bands'],dtype=float)
        if nominal.shape!=(3,) or bands.shape!=(3,) or not np.isfinite([nominal,bands]).all() or np.any(bands<0) or np.any(nominal-bands<=0) or (nominal[2]+bands[2])*max(BLOOD_FRACTION)>1:
            raise ValueError('Invalid positive physical input bands')
        if type(spec['seed']) is not int or spec['seed']<0 or type(spec['draws']) is not int or not 1<=spec['draws']<=100000:
            raise ValueError('Invalid deterministic sampling parameters')
        available={'input.po2':'scalar','input.hb':'scalar','input.fraction':'scalar','input.paths':'paths'}
        if not isinstance(spec['nodes'],list) or not spec['nodes']:
            raise ValueError('Expected ordered nodes')
        for node in spec['nodes']:
            if not isinstance(node,dict) or set(node)!={'id','module','inputs'}:
                raise ValueError('Invalid node schema')
            name=node['id'];module=node['module'];inputs=node['inputs']
            if not isinstance(name,str) or not name.isidentifier() or name in available or not isinstance(module,str) or module not in REGISTRY:
                raise ValueError('Unknown module or duplicate node')
            if not isinstance(inputs,dict) or set(inputs)!=set(REGISTRY[module][0]) or any(not isinstance(v,str) or v not in available for v in inputs.values()):
                raise ValueError('Missing or forward seam reference')
            if any(available[v]!=INPUT_TYPES[module][k] for k,v in inputs.items()):
                raise ValueError('Incompatible seam type')
            available[name]=OUTPUT_TYPES[module]
        if spec['output'] not in {n['id'] for n in spec['nodes'] if n['module']=='record'}:
            raise ValueError('Output must reference a record node')
        self.spec=spec;self.nominal=nominal;self.bands=bands

    def propagate(self,paths,inputs):
        records=[]
        for po2,hb,fraction in inputs:
            values={'input.po2':po2,'input.hb':hb,'input.fraction':fraction,'input.paths':paths}
            for node in self.spec['nodes']:
                values[node['id']]=REGISTRY[node['module']][1](**{k:values[v] for k,v in node['inputs'].items()})
            records.append(values[self.spec['output']])
        return np.asarray(records)

    def observe(self,paths):
        return observe(self,paths)


def observe(chain,paths):
    NOMINAL=chain.nominal;BANDS=chain.bands;SEED=chain.spec['seed'];DRAWS=chain.spec['draws']
    propagate=chain.propagate
    random=np.random.default_rng(SEED).uniform(-1.,1.,(DRAWS,3))
    joint=propagate(paths,NOMINAL+BANDS*random);single=[];hashes={}
    for axis in range(3):
        offsets=np.zeros((DRAWS,3));offsets[:,axis]=BANDS[axis]*random[:,axis]
        values=propagate(paths,NOMINAL+offsets);single.append(float(np.ptp(values[:,-1])))
        hashes[f'isolated_{axis}']=hashlib.sha256(values.tobytes()).hexdigest()
    nominal=propagate(paths,NOMINAL[None,:]);zero=propagate(paths,NOMINAL+0.*random)
    grid=propagate(paths,np.array([NOMINAL+BANDS*np.array(x) for x in itertools.product((-1,0,1),repeat=3)]))
    grid_single=[]
    for axis in range(3):
        values=[]
        for side in (-1,0,1):
            x=NOMINAL.copy();x[axis]+=side*BANDS[axis];values.append(x)
        grid_single.append(float(np.ptp(propagate(paths,np.array(values))[:,-1])))
    blood_null=propagate(paths,np.array([[60.,15.,0.]]))[0,-1]
    oxygen_low=propagate(paths,np.array([[58.,15.,1.]]))[0,-1]
    oxygen_high=propagate(paths,np.array([[62.,15.,1.]]))[0,-1]
    hashes.update(joint=hashlib.sha256(joint.tobytes()).hexdigest(),grid=hashlib.sha256(grid.tobytes()).hexdigest(),
                  zero=hashlib.sha256(zero.tobytes()).hexdigest())
    return {'hashes':hashes,'nominal_ratio':float(nominal[0,-1]),'joint_span':float(np.ptp(joint[:,-1])),
        'isolated_spans':single,'sum_isolated_spans':sum(single),'sample_interval_95':np.quantile(joint[:,-1],[.025,.975]).tolist(),
        'grid_span':float(np.ptp(grid[:,-1])),'grid_isolated_span_sum':sum(grid_single),
        'random_within_grid':bool(joint[:,-1].min()>=grid[:,-1].min() and joint[:,-1].max()<=grid[:,-1].max()),
        'finite_physical':bool(np.isfinite(joint).all() and np.all(joint[:,3]>0) and np.all(joint[:,3]<1) and np.all(joint[:,8:10]>0) and np.all(joint[:,8:10]<1)),
        'zero_band_exact':bool(np.all(zero==nominal)),'zero_blood_ratio':float(blood_null),
        'oxygen_direction':bool(oxygen_low<nominal[0,-1]<oxygen_high)}


def generate(spec):
    """Return deterministic standalone source; caller owns writing and execution."""
    validated=Chain(spec).spec
    encoded=json.dumps(validated,sort_keys=True,separators=(',',':'))
    return ('"""Generated conditional NIRS graph; no clinical calibration."""\n'
            'import json\nimport sys\nimport numpy as np\n'
            'from bodytwin.chains.nirs_spec_v1 import Chain\n'
            'SPEC=json.loads('+repr(encoded)+')\n'
            'if __name__ == "__main__":\n'
            '    with np.load(sys.argv[1]) as archive: paths=archive["paths_0"]\n'
            '    print(json.dumps(Chain(SPEC).observe(paths),sort_keys=True))\n')


if __name__ == '__main__':
    import runpy
    runpy.run_path(str(Path(__file__).resolve().parents[3]/'probes/optics/nirs_spec_probe.py'),run_name='__main__')
