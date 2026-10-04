"""Run declared dataset-free kernel profiles and reject changed predictions."""
import argparse
import copy
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np
from .load import ROOT, kernel, functions

CONTRACTS=ROOT/'fixtures/profile_contracts.json'

def require(condition,message):
    if not condition:raise AssertionError(message)

def synthetic_response():
    return dict(schema='conditional-edit-response-v1',case='synthetic-four-support',geometry_sha256='synthetic-geometry',
        spectral_error_N_per_mm=0.,baseline_error_N=0.,A=[[1,-1,-1],[1,1,-1],[1,-1,1],[1,1,1]],
        basis_N=[[.5],[-.5],[-.5],[.5]],H_N_per_mm=[[100.]],baseline_force_N=[25.]*4,
        w_N=[100.,0.,0.],baseline_height_mm=[0.]*4,predicted_fdi=[16,17,26,27],
        complete_active_calibration=True,radius_model='full_response')

def x14():
    with kernel('X14','seating.py') as m:
        x=np.array([-12.,-5.,7.,18.]);a=np.array([.02,-.02,.02,-.02]);b=np.array([.02,.02,-.02,-.02])
        qa=m.seat_vertices(x,a);qb=m.seat_vertices(x,b)
        lp=m.seat_lp(x,a)
        summary=lambda z:np.array([z.mean(),np.mean(z*z),z.min(),z.max()])
        bound=m.measured_seating_bound(x,a,.001)
        require(bound['max_gap_bound_mm'][0]<=qa['max_gap_mm']<=bound['max_gap_bound_mm'][1],'seating bound')
        return dict(lp_error_mm=abs(qa['max_gap_mm']-lp['max_gap_mm']),summary_identity_error=float(np.max(abs(summary(a)-summary(b)))),
            downstream_difference_mm=abs(qa['max_gap_mm']-qb['max_gap_mm']),minimum_extension='two chord-relative slopes',max_gap_mm=qa['max_gap_mm'])

def x51():
    with kernel('X51','code/acceptance.py') as m:
        p=m.plan();z=m.n_zero();g=m.n_zero(groups=4)
        # Scalar binomial arithmetic independently verifies the consumer risk.
        from math import comb
        scalar=sum(comb(p['n'],k)*.1**k*.9**(p['n']-k) for k in range(p['failures_allowed']+1))
        require(abs(scalar-p['false_acceptance_at_p_bad'])<1e-13,'independent binomial control')
        return dict(n=p['n'],failures_allowed=p['failures_allowed'],zero_failure_n=z['n'],four_group_n=g['n'],
                    three_survivor_risk_upper=m.upper(3,0),consumer_risk=scalar,power=p['power_at_p_good'])

def x56():
    with kernel('X56','code/motion_port.py') as m:
        z=np.array([-1.,1.]);C=np.array([[.001,.0001],[.0001,.0002]]);H=np.column_stack([np.ones(2),z]);W=np.eye(2)
        result=m.identify(z,W,H@C, np.zeros((2,2)))
        boxes=m.project_intervals([.01,-.001],[.02,.001],[-10.,0.,10.])
        for u in [.01,.02]:
            for theta in [-.001,.001]:
                values=u+np.array([-10,0,10])*theta
                require(np.all(values>=boxes[:,0]) and np.all(values<=boxes[:,1]),'outward region bound')
        return dict(matrix_error=float(np.max(abs(np.array(result['C'])-C))),missing_reference_status=m.identify(z,W,H@C)['status'],
                    missing_biology_status=m.local_threshold([1,2],True,True),region_bounds_mm=boxes.tolist())

def x59():
    funcs=functions('X59','code/analyze_r2.py',['contrast','witness','ecological_check'],{'np':np,'D':Decimal})
    rows=json.loads((ROOT/'fixtures/cunali2017_table1.json').read_text())['rows'];mapped=[]
    for r in rows:
        mapped.append(dict(system=r['system'],replica_mean_um=r['mean_1_um'],ct_mean_um=r['mean_2_um'],
            replica_sd_um=r['sd_1_um'],ct_sd_um=r['sd_2_um'],source_locator=r['locator']))
    contrasts=[float(funcs['contrast'](r)['bias_replica_minus_ct_um']) for r in mapped]
    w=funcs['witness']();coeff=funcs['ecological_check'](mapped)
    return dict(contrasts_um=contrasts,summary_identity_error=w['identity_max_abs_error_um'],
        classification_difference=w['downstream_classification_error_difference'],paired_sd_difference_um=w['downstream_paired_sd_difference_um'],
        minimum_extension=w['minimum_extension_for_nominal_normal_loa'],regression_gate_Amann=coeff[0]['external_coefficient_gate'],
        regression_gate_Dentsply_Sirona=coeff[1]['external_coefficient_gate'],source_locator=rows[0]['locator'],
        external_regression_errors=coeff[1]['errors'],individual_prediction='UNKNOWN')

def x63():
    with kernel('X63','code/assembly.py') as m:
        K=m.beam_matrix([10],exact=True);cantilever=m.rational_solve([r[2:] for r in K[2:]],[Fraction(1),Fraction(0)])
        lengths=[10,10,10];base=m.beam_matrix(lengths);zero=np.zeros(4);alt=np.array([.001,-.001,.001,-.001]);load=np.zeros(4)
        a=m.state(base,zero,300,load);b=m.state(base,alt,300,load)
        ea=m.exact_first_event(lengths,zero,300,[0,0,0,1]);eb=m.exact_first_event(lengths,alt,300,[0,0,0,1])
        require(ea['status']==eb['status']=='EXACT_RATIONAL_FIXED_MODEL','initial branch missing')
        require(ea['identity_residual_exact']=='0' and eb['identity_residual_exact']=='0','rational equilibrium')
        independent=m.independent_state(base,alt,300,load)
        require(np.max(abs(independent['bolt']-b['bolt']))<1e-5,'independent convex control')
        public=json.loads((ROOT/'fixtures/preload_publication.json').read_text())
        return dict(cantilever_compliance_fraction=str(cantilever[0]),closed_gap_identity_error=float(np.max(abs(a['gap']-b['gap']))),
            event_difference_N=abs(ea['first_event_N']-eb['first_event_N']),first_event_N=ea['first_event_N'],
            minimum_extension='regional preload and distortion state',preload_comparator_relative_error_percent=100*(public['first_mean_N']/public['tenth_mean_N']-1),
            physical_enclosure='MISSING')

def x68():
    with kernel('X68','code/enclosures.py') as m:
        poly=[.1,2.,0.];D=4.;lo=2.;hi=3.;r=m.certified_derivative(poly,D,lo,hi)
        lower,upper=map(Fraction,r['exact_rational_bounds']);c,b,_=map(lambda x:Fraction(float(x)),poly)
        A=b+2*c*Fraction(D)**2
        for i in range(101):
            d=Fraction(2)+Fraction(i,100);g=-2*A*d+4*c*d**3
            require(lower<=g<=upper,'exact derivative control')
        return dict(tested_rational_points=101,physical_remainder=r['physical_model_remainder_bound'],derivative_interval=r['dT_dDf_Ncm_per_mm'],proof=r['proof'])

def x74():
    with kernel('X74','code/optics.py') as m:
        a=[50.,2.,10.];b=[55.,4.,12.];value=m.de00(a,b)
        bounds=m.interval_de00([[49.9,50.1],[1.9,2.1],[9.9,10.1]],b)
        require(bounds[0]<=value<=bounds[1],'colour interval')
        return dict(identity_deltaE=m.de00(a,a),symmetry_error=abs(value-m.de00(b,a)),deltaE=value,interval=bounds,physical_validation='UNKNOWN')

def x82():
    with kernel('X82','code/height_force.py') as m:
        s=synthetic_response();a=m.predict(s,m.make_edit(s,0,.02));b=m.predict(s,m.make_edit(s,0,0))
        require(a['status']=='CONDITIONAL_PASS','conditional contact gate')
        A=np.array(s['A']);summary=lambda f:np.r_[A.T@f,np.count_nonzero(np.array(f)>0)]
        unknown=m.predict(dict(s,complete_active_calibration=False),m.make_edit(s,0,.02))
        return dict(force_N=a['force_N'],summary_identity_error=float(np.max(abs(summary(a['force_N'])-summary(b['force_N'])))),
            downstream_difference_N=float(np.max(abs(np.array(a['force_N'])-b['force_N']))),missing_calibration_status=unknown['status'],
            joint_l2_error_radius_N=a['joint_l2_error_radius_N'],minimum_extension='local force response on the balanced nullspace',
            mathematical_scope='Conditional synthetic response; numerical roundoff is not formally enclosed')

def x85():
    with kernel('X85','code/lab_compare.py') as m:
        f=dict(case='synthetic',frame='synthetic-frame',geometry_sha256='synthetic-geometry',force_interval_N=[[7.,9.],[7.,9.]])
        y=dict(case=f['case'],frame=f['frame'],geometry_sha256=f['geometry_sha256'],unit='N',matched_load_rate_history=True,reaction_N=[8.,8.],reaction_error_N=[.1,.1])
        normal=m.compare(f,y);bad=copy.deepcopy(y);bad['reaction_N'][0]+=3;fault=m.compare(f,bad)
        identity='ACCEPT'
        try:m.compare(f,dict(y,geometry_sha256='changed'))
        except ValueError:identity='REJECT'
    with kernel('X85','code/height_uncertainty.py') as m:
        probes=dict(schema='regional-actuation-probes-v1',units={'height':'mm','force':'N'},observation_operator='calibrated_regional_axial_reactions',
            linear_support_closure=True,matched_load_rate_history=True,height_error_mm=.001,step_mm=.01,
            baseline_force_N=[10.,10.],plus_force_N=[[10.09,9.91]],minus_force_N=[[9.91,10.09]],raw_force_channel_bound_N=0.,
            geometry_sha256='synthetic',basis_sha256='synthetic',case='synthetic',frame='synthetic',kind='SIMULATION')
        bounded=m.fit_bounded_heights(probes);denominator='ACCEPT'
        try:m.fit_bounded_heights(dict(probes,height_error_mm=.01))
        except ValueError:denominator='REJECT'
        require(bounded['beta_N_per_mm'][0]>0,'missing height error enclosure')
    return dict(heldout_status=normal['status'],injected_status=fault['status'],invalid_identity=identity,denominator_zero=denominator,
                measurement='SYNTHETIC; no lab measurement performed',height_column_error_N_per_mm=bounded['beta_N_per_mm'][0])

RUNNERS={'X14':x14,'X51':x51,'X56':x56,'X59':x59,'X63':x63,'X68':x68,'X74':x74,'X82':x82,'X85':x85}

def check(demo,result,contract):
    c=contract
    if demo=='X14':
        require(result['lp_error_mm']<=c['lp_parity_tolerance_mm'],'LP parity');require(result['summary_identity_error']==0,'summary identity');require(result['downstream_difference_mm']>=c['minimum_downstream_difference_mm'],'seating sufficiency')
    elif demo=='X51':
        for key in ['n','failures_allowed','zero_failure_n','four_group_n']:require(result[key]==c[key],key)
        require(abs(result['three_survivor_risk_upper']-c['three_survivor_risk_upper'])<=c['tolerance'],'risk upper')
    elif demo=='X56':
        require(result['matrix_error']<=c['matrix_error_tolerance'],'metrology recovery')
        for key in ['missing_reference_status','missing_biology_status']:require(result[key]==c[key],key)
    elif demo=='X59':
        require(np.max(abs(np.array(result['contrasts_um'])-c['contrasts_um']))<=c['tolerance_um'],'published contrasts')
        for key in ['summary_identity_error','classification_difference','regression_gate_Amann','regression_gate_Dentsply_Sirona']:require(result[key]==c[key],key)
    elif demo=='X63':
        require(result['cantilever_compliance_fraction']==c['cantilever_compliance_fraction'],'beam control');require(result['closed_gap_identity_error']==0,'closed gaps');require(result['event_difference_N']>=c['minimum_event_difference_N'],'event sufficiency')
        require(abs(result['preload_comparator_relative_error_percent']-c['preload_comparator_relative_error_percent'])<=c['relative_error_tolerance_percent'],'published preload')
    elif demo=='X68':
        require(result['tested_rational_points']==c['tested_rational_points'],'derivative coverage');require(result['physical_remainder']==c['physical_remainder'],'missing physical remainder')
        lo,hi=result['derivative_interval']
        cc=Fraction(float(.1));bb=Fraction(2);AA=bb+2*cc*16
        for i in range(101):
            d=Fraction(2)+Fraction(i,100);v=-2*AA*d+4*cc*d**3
            require(Fraction(float(lo))<=v<=Fraction(float(hi)),'reported derivative interval')
    elif demo=='X74':
        require(result['identity_deltaE']==c['identity_deltaE'],'colour identity');require(result['symmetry_error']<=c['symmetry_tolerance'],'colour symmetry')
        interval=np.asarray(result['interval'],float);value=float(result['deltaE'])
        require(interval.shape==(2,) and np.isfinite(interval).all() and np.isfinite(value),'finite colour result')
        require(0<=interval[0]<=value<=interval[1],'reported colour enclosure')
    elif demo=='X82':
        require(np.max(abs(np.array(result['force_N'])-c['force_N']))<=c['force_tolerance_N'],'force oracle');require(result['summary_identity_error']==0,'wrench identity');require(result['downstream_difference_N']>=c['minimum_downstream_difference_N'],'force sufficiency');require(result['missing_calibration_status']==c['missing_calibration_status'],'missing calibration')
    elif demo=='X85':
        for key in ['heldout_status','injected_status','invalid_identity','denominator_zero']:require(result[key]==c[key],key)

MUTATIONS={'X14':('lp_error_mm',1.),'X51':('n',47),'X56':('matrix_error',1.),'X59':('contrasts_um',[1000.]*8),
           'X63':('cantilever_compliance_fraction','1/7000'),'X68':('derivative_interval',[0.,1.]),
           'X74':('identity_deltaE',1.),'X82':('force_N',[1000.]*4),'X85':('heldout_status','UNKNOWN')}

def verify_freeze():
    f=json.loads((ROOT/'FROZEN_PREDICTIONS.json').read_text())
    for name,expected in f['files'].items():
        require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,'frozen input or kernel changed: '+name)

def figure(rows,path):
    contrasts=rows['X59']['contrasts_um']
    labels=['Amann margin','Amann axial','Amann transition','Amann occlusal','Dentsply margin','Dentsply axial','Dentsply transition','Dentsply occlusal']
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="760" height="390" viewBox="0 0 760 390">',
           '<rect width="760" height="390" fill="white"/>','<g font-family="sans-serif" font-size="14" fill="#17252c">',
           '<text x="22" y="27">Published bench group means: replica minus micro-CT</text>',
           '<text x="22" y="49">Cunali 2017, Table 1; 10 copings/system; per surface region</text>']
    for i,(label,value) in enumerate(zip(labels,contrasts)):
        y=78+32*i;x=230+min(0,value)*6
        parts += [f'<text x="22" y="{y+14}">{label}</text>',f'<rect x="{x}" y="{y}" width="{abs(value)*6}" height="18" fill="#287f8e"/>',f'<text x="{240+max(0,value)*6}" y="{y+14}">{value:.2f}</text>']
    parts += ['<path d="M230 72V333" stroke="#17252c"/>','<text x="230" y="356">Regional difference (micrometres); individual agreement UNKNOWN</text>','</g></svg>']
    path.write_text('\n'.join(parts)+'\n')

def run(selected,output):
    verify_freeze();contracts=json.loads(CONTRACTS.read_text())['profiles'];rows={}
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    for ident in selected:
        if ident not in RUNNERS:raise ValueError(ident+' has source-only code; supply external inputs and port its full pipeline as described in docs/DATA.md')
        result=RUNNERS[ident]();check(ident,result,contracts[ident]);mutated=copy.deepcopy(result);key,value=MUTATIONS[ident];mutated[key]=value
        rejected=False
        try:check(ident,mutated,contracts[ident])
        except AssertionError:rejected=True
        require(rejected,'injected wrong result passed '+ident)
        rows[ident]=dict(result,status='PASS',profile=contracts[ident]['profile'],claim_type='capability',resolution=contracts[ident]['resolution'],injected_wrong_result_rejected=True)
    (output/'results.json').write_text(json.dumps({'physical_validation':'UNKNOWN','full_original_pipelines_replayed':0,'kernel_profiles':rows},indent=2,allow_nan=False)+'\n')
    if 'X59' in rows:figure(rows,output/'published_method_contrasts.svg')
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('ids',nargs='*');p.add_argument('--output',type=Path);p.add_argument('--list',action='store_true');args=p.parse_args()
    if args.list:
        for ident in RUNNERS:print(ident, json.loads(CONTRACTS.read_text())['profiles'][ident]['profile'])
        return
    output=args.output or Path(tempfile.mkdtemp(prefix='dental-demo-'));selected=args.ids or list(RUNNERS)
    rows=run(selected,output);print(json.dumps({'status':'PASS','kernel_profiles':len(rows),'output':str(output),'physical_validation':'UNKNOWN'}))

if __name__=='__main__':main()
