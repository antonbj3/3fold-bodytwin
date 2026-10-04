"""Conditional branch-preference certificates; no empirical skin admission.

Exact box envelopes and angle enclosures reuse the reviewed
PROOF_LANE_AVLANKNING_20261001/decision.py and mode3.py algebra. The built-in bank
is ONLY the declared finite-extension antiplane fixture (square [-1,1]^2,
extension 1/4, Dirichlet trace y). It is not a generic skin or mode-I provider.
Supplied interval and front certificates are conditional on producer premises.
"""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math

LABELS = ('REVIEW', 'DEFLECTION', 'UNCERTAIN')
ENERGIES = ('gp', 'gd', 'gamma_layer', 'gamma_interface')
BANK_SHA256 = '3c1ce365002d55dddbb4badec5ac02908280c2e968baf984ffe736a68b9783c7'
FIXTURE_CONTEXT = 'modeIII_square_-1_1_trace_y_extension_1/4_n64'


def _q(x):
    if type(x) not in (int, str, Q):
        raise TypeError('exact integer or rational string required; no floats/bools')
    return Q(x)


def _band(x, *, positive=False):
    if x is None:
        return None
    if not isinstance(x, (list, tuple)) or len(x) != 2:
        raise ValueError('two ordered endpoints required')
    lo, hi = map(_q, x)
    if not 0 <= lo <= hi or (positive and lo <= 0):
        raise ValueError('nonnegative ordered interval required')
    return lo, hi


def _hash(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True, ensure_ascii=False,
        separators=(',', ':'), allow_nan=False, default=str).encode()).hexdigest()


def _label(lo, hi):
    return 'DEFLECTION' if lo > 0 else 'REVIEW' if hi < 0 else 'UNCERTAIN'


def _output(decision, margin, *, source, context, missing=(), unit='(J/m²)^2'):
    def port(value, units, scope):
        return dict(value=value, unit=units,
            status='UNKNOWN' if value is None else 'SYNTETISKT',
            uncertainty=None if value is None else {'kind': 'conditional model; no probability'},
            source=source, scope=scope, joint_id=None)
    return dict(decision=decision, physical_admission=False, model_context=context,
        margin=None if margin is None else list(map(str, margin)), margin_unit=unit,
        missing=list(missing),
        path_preference=port(None if decision == 'UNCERTAIN' else decision, '1',
            'Conditional preference among declared branches; growth/onset unknown'),
        predicted_depth=port(None, 'm', 'Independently supplied cut depth required'),
        biological_damage_width=port(None, 'm', 'Independent injury observation required'))


def _certificate(result, inputs, *, rule, premises):
    result['certificate'] = dict(schema='avlankning-certificate-v1',
        grade='UNKNOWN' if result['decision'] == 'UNCERTAIN' else 'CONDITIONAL_MODEL',
        input_sha256=_hash(inputs), rule=rule, premises=premises,
        numerical_arithmetic='exact rational; certified angular enclosures for fixture',
        empirical_error_bound=None, review_state='PENDING_INDEPENDENT_REVIEW',
        end_to_end_lean_proof=False)
    result['certificate']['result_sha256'] = _hash(result)
    return result


def path_preference(gp, gd, gamma_layer, gamma_interface, *, unit='J/m²',
                    mode='evidence', source='unbound', model_context=None):
    """Retains the old four-interval API; output uses the lane's canonical labels."""
    if mode not in ('evidence', 'scenario') or unit != 'J/m²':
        raise ValueError('evidence/scenario mode and common J/m² required')
    if not isinstance(source, str) or not source.strip():
        raise ValueError('source required')
    inputs=dict(zip(ENERGIES, (gp, gd, gamma_layer, gamma_interface)))
    bands=tuple(_band(x) for x in inputs.values())
    missing=[];margin=None;decision='UNCERTAIN'
    if mode == 'evidence':
        missing.append('reviewed_empirical_mode_geometry_material_and_onset')
    elif not isinstance(model_context, str) or not model_context.strip():
        missing.append('common_load_geometry_extension_context')
    elif any(x is None for x in bands):
        missing.append('energy_or_toughness_enclosure')
    elif bands[0][0] <= 0 or bands[2][0] <= 0:
        missing.append('strict_positive_gp_and_gamma_layer')
    else:
        p,d,l,i=bands
        margin=(d[0]*l[0]-p[1]*i[1],d[1]*l[1]-p[0]*i[0])
        decision=_label(*margin)
        if decision == 'UNCERTAIN':missing.append('strict_preference_margin')
    result=_output(decision,margin,source=source,context=model_context,missing=missing)
    return _certificate(result,dict(inputs,mode=mode,unit=unit,source=source,context=model_context),
        rule='sign enclosure of Gd*Gamma_layer-Gp*Gamma_interface',
        premises=['supplied enclosures valid over whole common state box',
                  'two declared branches and positive normalization; conditional on growth'])


@lru_cache(None)
def _pi_bounds():
    def atan(x):
        s=sum((-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(40))
        return s,s+x**81/81
    a,b=atan(Q(1,5));c,d=atan(Q(1,239))
    return 16*a-4*d,16*b-4*c


@lru_cache(None)
def _sine(deg):
    if deg == 0:return Q(0),Q(0)
    if deg == 90:return Q(1),Q(1)
    pl,ph=_pi_bounds();xl,xh=deg*pl/180,deg*ph/180
    def bounds(x):
        s=sum((-1)**k*x**(2*k+1)/math.factorial(2*k+1) for k in range(12))
        return s,s+x**25/math.factorial(25)
    lo,_=bounds(xl);_,hi=bounds(xh);scale=10**18
    return Q((lo*scale).__floor__(),scale),Q((hi*scale).__ceil__(),scale)


def _angular(lo, hi, amplitude):
    low=Q(0) if lo <= 0 <= hi else min(abs(lo),abs(hi))
    high=max(abs(lo),abs(hi))
    return 1+amplitude*_sine(low)[0]**2,1+amplitude*_sine(high)[1]**2


def _box(b):
    if not isinstance(b, dict) or set(b) != {'r','gamma','theta'}:
        raise ValueError('r, gamma, theta coordinates required')
    out={k:tuple(map(_q,v)) for k,v in b.items()}
    if any(len(v)!=2 or v[0]>v[1] for v in out.values()):
        raise ValueError('ordered coordinate pairs required')
    if out['r'][0]<=0 or out['gamma'][0]<0 or out['theta'][0]<-90 or out['theta'][1]>90:
        raise ValueError('positive modulus, nonnegative gamma, angles within [-90,90]')
    return out


def _dump(b):
    return {k:list(map(str,b[k])) for k in ('r','gamma','theta')}


def _bank():
    data=Path(__file__).with_name('avlankning_fixture_bank.json').read_bytes()
    if hashlib.sha256(data).hexdigest()!=BANK_SHA256:
        raise ValueError('fixture bank hash mismatch')
    return json.loads(data)


def _release(base, branch, r):
    lower=[];upper=[]
    for x in r:
        lower.append(base['work']-base['Vleft']/x-base['Vright']
                     -branch['Uleft']*x-branch['Uright'])
        upper.append(base['Uleft']*x+base['Uright']-branch['work']
                     +branch['Vleft']/x+branch['Vright'])
    return max(Q(0),min(lower)),max(upper)


def _leaf(b, trials, amplitude):
    p=_release(trials[0],trials[1],b['r'])
    d=_release(trials[0],trials[2],b['r'])
    if p[0]<=0:return 'UNCERTAIN',None
    w=_angular(*b['theta'],amplitude)
    margin=(d[0]/p[1]*w[0]-b['gamma'][1],
            d[1]/p[0]*w[1]-b['gamma'][0])
    return _label(*margin),margin


def fixture_preference(box, *, anchor, closure, max_depth=6, mode='evidence'):
    b=_box(box)
    if mode not in ('evidence','scenario'):raise ValueError('invalid mode')
    if closure not in ('ISOTROPIC','SYNTHETIC_SIN2'):raise ValueError('explicit supported closure required')
    if type(max_depth) is not int or not 0<=max_depth<=10:raise ValueError('depth 0..10')
    data=_bank();anchor=str(_q(anchor))
    if anchor not in data:raise ValueError('verified fixture anchor required')
    trials=[]
    for rec,branch in zip(data[anchor],('base','penetrate','deflect')):
        if rec['branch']!=branch:raise ValueError('branch order mismatch')
        t={k:_q(rec[k]) for k in ('Uleft','Uright','Vleft','Vright','work')}
        if any(t[k]<0 for k in ('Uleft','Uright','Vleft','Vright')):raise ValueError('invalid energy coefficient')
        trials.append(t)
    inputs=dict(box=_dump(b),anchor=anchor,closure=closure,max_depth=max_depth,mode=mode)
    if mode=='evidence':
        return _certificate(_output('UNCERTAIN',None,source=BANK_SHA256,context=FIXTURE_CONTEXT,
            missing=['fixture_not_empirical_skin_provider'],unit='1'),inputs,
            rule='evidence admission remains unknown',premises=['fixture scope only'])
    counts={k:Q(0) for k in LABELS};leaves=[];evaluations=0
    amplitude=Q(0 if closure=='ISOTROPIC' else 1)
    def visit(cell,weight,depth):
        nonlocal evaluations
        evaluations+=1;status,margin=_leaf(cell,trials,amplitude)
        dims=('r','gamma','theta')
        dim=next((dims[(depth+i)%3] for i in range(3)
                  if cell[dims[(depth+i)%3]][0]<cell[dims[(depth+i)%3]][1]),None)
        if status!='UNCERTAIN' or depth==max_depth or dim is None:
            counts[status]+=weight
            leaves.append(dict(box=_dump(cell),weight=str(weight),status=status,
                margin=None if margin is None else list(map(str,margin))))
        else:
            low,high=cell[dim];mid=(low+high)/2
            for pair in ((low,mid),(mid,high)):
                child=dict(cell);child[dim]=pair;visit(child,weight/2,depth+1)
    visit(b,Q(1),0)
    decision=next((k for k in LABELS[:2] if counts[k]==1),'UNCERTAIN')
    result=_output(decision,_leaf(b,trials,amplitude)[1],source=BANK_SHA256,
        context=FIXTURE_CONTEXT,unit='1',missing=['whole_box_strict_preference'] if decision=='UNCERTAIN' else [])
    result.update(fractions={k:str(v) for k,v in counts.items()},leaves=leaves,
        fraction_kind='product coordinate volume, not probability',bound_evaluations=evaluations,
        routing_license='concave lower/convex upper release envelopes; angular split at zero',
        instance_statistics=dict(positive_r=True,angle_crosses_zero=b['theta'][0]<=0<=b['theta'][1],
                                 positive_release_denominator=_release(trials[0],trials[1],b['r'])[0]>0))
    return _certificate(result,inputs,rule='reviewed primal/dual envelopes + strict signs + full partition',
        premises=['verified upstream fixture trial fields, hash-bound bank',
                  'declared fixed mode/load/geometry/extension',
                  'synthetic directional resistance or isotropic closure',
                  'no onset, depth or biological injury inference'])


def front_preference(cells, *, front_complete=False, independent_strips=False,
                     coupling_margin_bound=None, source='unbound', model_context=None,
                     mode='evidence', parent_enclosures=None):
    """Local-front enclosure. Coupling bound is absolute per-cell margin correction.

    This is an intersection over spatial cells, not a probability of deflection.
    Valid full coverage and the coupling bound are producer premises. For real
    interacting fronts a producer must certify them; this routine cannot do so.
    """
    if type(front_complete) is not bool or type(independent_strips) is not bool:
        raise TypeError('explicit boolean coverage and independence required')
    if not isinstance(cells,list) or not cells:raise ValueError('nonempty front cells required')
    ids=[c['cell_id'] for c in cells]
    if any(not isinstance(i,str) or not i for i in ids) or len(set(ids))!=len(ids):
        raise ValueError('unique nonempty cell IDs required')
    areas=[_q(c['area_m2']) for c in cells]
    if min(areas)<=0:raise ValueError('positive area per cell required')
    if coupling_margin_bound is not None:
        delta=_q(coupling_margin_bound)
        if delta<0:raise ValueError('nonnegative absolute coupling bound required')
        if independent_strips and delta!=0:raise ValueError('independent strips have zero coupling')
    else:delta=Q(0) if independent_strips else None
    rows=[];counts={k:Q(0) for k in LABELS};total=sum(areas)
    for c,area in zip(cells,areas):
        if parent_enclosures is not None:
            for k in ENERGIES:
                parent=_band(parent_enclosures[k]);child=_band(c[k])
                if parent is None or child is None or not parent[0]<=child[0]<=child[1]<=parent[1]:
                    raise ValueError('child enclosure outside declared parent')
        r=path_preference(**{k:c[k] for k in ENERGIES},mode=mode,
            source=source,model_context=model_context)
        if delta is None or not front_complete or r['margin'] is None:
            status='UNCERTAIN';margin=None
        else:
            lo,hi=map(Q,r['margin']);margin=(lo-delta,hi+delta);status=_label(*margin)
        counts[status]+=area/total
        rows.append(dict(cell_id=c['cell_id'],area_fraction=str(area/total),decision=status,
            margin=None if margin is None else list(map(str,margin)),local_certificate=r['certificate']))
    decision=next((k for k in LABELS[:2] if counts[k]==1),'UNCERTAIN')
    missing=[]
    if not front_complete:missing.append('certified_full_front_coverage')
    if delta is None:missing.append('front_interaction_remainder')
    if decision=='UNCERTAIN':missing.append('uniform_strict_local_preference')
    margin=None if any(c['margin'] is None for c in rows) else (
        min(Q(c['margin'][0]) for c in rows),max(Q(c['margin'][1]) for c in rows))
    result=_output(decision,margin,source=source,context=model_context,missing=missing)
    result.update(cells=rows,area_fractions={k:str(v) for k,v in counts.items()},
        fraction_kind='declared spatial area, not probability',
        coupling_margin_bound=None if delta is None else str(delta),
        mixed_front=counts['REVIEW']>0 and counts['DEFLECTION']>0)
    return _certificate(result,dict(cells=cells,front_complete=front_complete,
        independent_strips=independent_strips,coupling_margin_bound=coupling_margin_bound,
        source=source,context=model_context,mode=mode,parent_enclosures=parent_enclosures),
        rule='intersection of local margin enclosures plus bounded coupling',
        premises=['full registered area cover', 'local G and Gamma bounds over each entire cell',
                  'independent strips OR certified absolute interaction correction',
                  'local path observable under this model; unmodeled front physics excluded'])


def evaluate_request(request, *, mode='evidence'):
    req=dict(request);kind=req.pop('kind');req.pop('id',None)
    if 'mode' in req:raise ValueError('mode is selected by consumer, not by request')
    if kind=='finite_extension_fixture':out=fixture_preference(mode=mode,**req)
    elif kind=='supplied_intervals':out=path_preference(mode=mode,**req)
    elif kind=='front_cells':out=front_preference(mode=mode,**req)
    else:raise ValueError('unsupported preference provider')
    if 'id' in request:out['id']=request['id']
    out['certificate']['request_sha256']=_hash(dict(request=request,consumer_mode=mode))
    out['certificate'].pop('result_sha256')
    out['certificate']['result_sha256']=_hash(out)
    return out


def incision_path_gate(config, *, mode='evidence'):
    """Optional current-chain hook. Unknown/deflection blocks straight continuation."""
    if config is None:return None
    if not isinstance(config,dict) or set(config)!={'requests'}:
        raise ValueError('avlankning config must contain requests only')
    requests=config['requests']
    if not isinstance(requests,list) or not requests:raise ValueError('nonempty requests required')
    ids=[r.get('id') for r in requests]
    if any(not isinstance(i,str) or not i for i in ids) or len(set(ids))!=len(ids):
        raise ValueError('unique nonempty request IDs required')
    rows=[evaluate_request(r,mode=mode) for r in requests]
    return dict(schema='bodytwin-avlankning-gate-v1',mode=mode,boxes=rows,
        blocks_straight_chain=any(r['decision']!='REVIEW' for r in rows),
        physical_admission=False,scope='opt-in conditional path check; supplied depth and injury remain independent')


def verify_integrity(result):
    """Integrity only; this does not validate producer premises or scientific scope."""
    copy=json.loads(json.dumps(result));expected=copy['certificate'].pop('result_sha256')
    return _hash(copy)==expected
