"""Portable, scoped source operators with synthetic or public aggregate inputs."""
import argparse
from collections import defaultdict
from contextlib import contextmanager
import copy
from fractions import Fraction as Q
import hashlib
import importlib
import itertools
import json
import math
from pathlib import Path
import sys
import tempfile
import time
import warnings
import xml.etree.ElementTree as ET
import zipfile
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, linprog as scipy_linprog
from .load import ROOT, functions, kernel

CONTRACTS=ROOT/'fixtures/capability_contracts.json'

def fixture(name): return json.loads((ROOT/'fixtures'/name).read_text())
def require(value, message):
    if not value: raise AssertionError(message)
def linprog(*args,**kw):
    kw['options']=dict(kw.get('options',{}),threads=1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return scipy_linprog(*args,**kw)
def defs(demo,path,names,extra=None):
    return functions(demo,path,names,dict(np=np,Q=Q,F=Q,math=math,json=json,Path=Path,defaultdict=defaultdict,itertools=itertools,warnings=warnings,linprog=linprog,**(extra or {})))

@contextmanager
def package(demo,name):
    oldpath=list(sys.path);saved={k:v for k,v in sys.modules.items() if k==name or k.startswith(name+'.')}
    for k in saved:del sys.modules[k]
    sys.path.insert(0,str(ROOT/'implementations'/demo))
    try:yield importlib.import_module(name)
    finally:
        sys.path[:]=oldpath
        for k in list(sys.modules):
            if k==name or k.startswith(name+'.'):del sys.modules[k]
        sys.modules.update(saved)

def cad():
    ns=defs('CAD','code/three_mf.py',['tag','write','readback'],dict(ET=ET,zipfile=zipfile,NS='http://schemas.microsoft.com/3dmanufacturing/core/2015/02',clean=lambda x:x))
    vertices=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
    faces=np.array([[0,1,2],[0,2,3],[0,3,1]]);roles=np.array([0,1,2])
    with tempfile.TemporaryDirectory(prefix='dental-cad-') as tmp:
        p=Path(tmp)/'synthetic.3mf';ns['write'](p,vertices,faces,roles,{'kind':'SYNTHETIC_SERIALIZATION_FIXTURE'})
        v,f,r=ns['readback'](p)
        with zipfile.ZipFile(p) as z:raw={k:z.read(k) for k in z.namelist()}
        raw['3D/3dmodel.model']=raw['3D/3dmodel.model'].replace(b'unit="millimeter"',b'unit="meter"')
        bad=Path(tmp)/'wrong.3mf'
        with zipfile.ZipFile(bad,'w') as z:
            for k,value in raw.items():z.writestr(k,value)
        rejected=False
        try:ns['readback'](bad)
        except ValueError:rejected=True
    cycles=defs('CAD','code/margin.py',['cycles'])['cycles']
    theta=np.arange(8)*math.pi/4
    cycle_vertices=np.c_[np.cos(theta),np.sin(theta),np.zeros(8)]
    edges=np.c_[np.arange(8),np.roll(np.arange(8),-1)]
    require(len(cycles(edges,cycle_vertices)[0])==1,'closed source cycle')
    require(len(cycles(edges[:-1],cycle_vertices)[0])==0,'open cycle must fail')
    a=np.array([3.,2.5,2.,1.5]);b=np.array([3.,2.,2.5,1.5])
    summary=lambda r:np.array([r[0],r[-1],r.mean(),3.,(r[0]-r[-1])/3])
    return dict(coordinate_error_mm=float(abs(v-vertices).max()),face_identity=np.array_equal(f,faces),role_identity=np.array_equal(r,roles),invalid_unit_rejected=rejected,summary_identity_error=float(abs(summary(a)-summary(b)).max()),waist_rebound_difference_mm=float(max(0,np.diff(b).max())-max(0,np.diff(a).max())),evidence='Synthetic export/cycle replay; full CAD anatomy replay requires external inputs',resolution='PER_POINT')

def crown_prep():
    fn=defs('PROOF_LANE_CROWN_FIX_PREP','code/core.py',['normals_exact'])['normals_exact']
    v=np.array([[0.,0,0],[1,0,0],[0,1,0]]);f=np.array([[0,1,2]])
    good=fn(v,f,[0,0,1]);bad=fn(v,f[:,::-1],[0,0,1])
    top=defs('PROOF_LANE_CROWN_FIX_PREP_R2','code/lower_topology.py',['topology'])['topology']
    faces=[[0,1,2],[0,2,3]];one=top(faces);extra=top(faces+[[10,11,12]])
    return dict(valid_normals=good['pass_exact'],flipped_normal_rejected=not bad['pass_exact'],valid_lower_topology=one['all_pass'],isolated_component_rejected=not extra['all_pass'],complete_original_chain=0,original_denominator=18,evidence='Rational stored-coordinate normals and lower-section incidence; original reviewed complete chain remains 0/18',physical_fit='UNKNOWN',resolution='PER_POINT/PER_SURFACE_REGION')

def x87():
    ns=defs('X87','science.py',['guide','classify','controls'],dict(brentq=brentq))
    errors=[];budgets=[]
    for p in fixture('guide_publication.json')['profiles']:
        g=ns['guide'](p,.05);c=ns['controls'](p,.05,g)
        errors.append(c['error_mm']);budgets.append(dict(guide=p['id'],bound_mm=g['guide_budget_mm']))
        require(c['injected_plus1mm_rejected'],'guide numerical control injection')
    return dict(max_control_error_mm=max(errors),zero_residual_equality=ns['classify']({'lower_mm':0.,'upper_mm':0.},0,strict=True),positive_residual_equality=ns['classify']({'lower_mm':2.,'upper_mm':2.},2),budgets=budgets,physical_margin='UNKNOWN',resolution='PHENOMENOLOGICAL',debt='Replace sample-as-population moments and missing signed pose/anatomical closure with paired independent measurements')

def x88():
    public=fixture('pulpotomy_publication.json')
    with kernel('X88','code/model.py') as m:
        a=m.r1_result(public['primary'],public['external']);rows=a['rows']
        ps=[Q(r['k'],r['n']) for r in rows]
        pmf=m.poisson_binomial_exact(ps);control=m.enumeration_control(ps)
        w=a['sufficiency']
        return dict(counts=[[r['k'],r['n']] for r in rows],pooled=[a['pooled']['k'],a['pooled']['n']],max_interval_control_error=a['interval_control_max_error'],summary_identity_error=w['identity_error'],category_difference_exact=w['downstream_difference_exact'],panel_control_error_exact=str(max(abs(x-y) for x,y in zip(pmf,control))),target_site_risk='UNKNOWN',rows=rows,missing_category_count=1,source_locator=public['source_locator'],resolution='POPULATION',scope='Retrospective cohort transcription and binomial intervals; no treatment choice or individual calibration')

def x89():
    with kernel('X89','code/calibration.py') as m:out=m.calibration_controls()
    ns=defs('X89','code/codesign.py',['compliance'])
    tool=dict(diameter_mm=.6,neck_reach_mm=2.7,gauge_mm=22,shank_mm=3)
    E=572000.;c=ns['compliance'](tool,E);n=3.;L=22.3
    q=quad(lambda x:(L-x)**2*64/(math.pi*E*(3 if x<L-n else .6)**4),0,L-n,epsabs=1e-13)[0]+quad(lambda x:(L-x)**2*64/(math.pi*E*.6**4),L-n,L,epsabs=1e-13)[0]
    s=out['sufficiency']
    return dict(upper_rational=out['baseline']['upper_rational'],missing_curvature='UNKNOWN_NO_NONLINEAR_ENCLOSURE',summary_identity_error=s['identity_error_mm'],midpoint_difference_rational=s['midpoint_difference_rational'],max_quadrature_error_mm_N=abs(c-q),published_RMS_um=[r['mean_RMS_um'] for r in fixture('milling_publication.json')['rows']],physical_qualification='UNKNOWN',synthetic_displacement_bound=out['baseline'],resolution='PER_POINT; published RMS PER_SURFACE_REGION',debt='Installed force/pose/profile/curvature measurements, CAM, runout and batch sintering')

def x96():
    fn=defs('X98','implant_safety/vendor/x8_geometry.py',['project_cylinder'])['project_cylinder']
    p=np.array([[12.,3.5,0]]);e=np.zeros(3);a=np.array([1.,0,0])
    nominal=float(np.linalg.norm(p-fn(p,e,a,10,2)))
    full=float(np.linalg.norm(p-fn(p,e,a,12,2)))
    axis=float(np.linalg.norm(p-fn(p,np.array([10.,0,0]),a,2,0)))
    # Both hypothetical tools share the original cylinder, maximum added depth,
    # and pose; their appended radial occupancy differs.
    summary_full=np.r_[nominal,e,a,2.,12.-10.]
    summary_axis=np.r_[nominal,e,a,2.,np.linalg.norm([2.,0,0])]
    return dict(nominal_gap_mm=nominal,full_radius_gap_mm=full,minimal_tool_gap_mm=min(nominal,axis),summary_identity_error=float(abs(summary_full-summary_axis).max()),shape_gap_difference_mm=abs(min(nominal,axis)-full),nominal_label_is_actual=False,published_depth_examples_mm=[p['extra_mm'] for p in fixture('drill_protocols.json')['protocols']],physical_safety='UNKNOWN',resolution='PER_POINT; protocol depth PHENOMENOLOGICAL',minimum_extension='Retain radial occupancy along the appended tool segment and the manufacturer depth datum')

def synthetic_bite():
    meta=dict(case_id='synthetic',frame_id='synthetic-mm',pose_id='synthetic-static',unit='mm')
    xy=np.array([[0.,0],[1,0],[1,1],[0,1]]);ff=np.array([[0,1,2],[0,2,3]])
    v=np.r_[np.c_[xy,np.ones(4)],np.c_[xy,np.zeros(4)]]
    f=np.r_[ff,[[4,6,5],[4,7,6]]];roles=np.array([0,0,1,1])
    b=dict(schema='registered-bite-v1',**meta,xy_mm=xy,faces=ff,antagonist_z_mm=np.full(4,1.05),reference_gap_mm=np.full(4,.05),source_locator='synthetic mathematical fixture',evidence_kind='SIMULATION',teeth_fdi=[16,26,17],pair_model='AXIAL_POSITIVE_SUPPORT_V1',tooth_pairs=[dict(upper_fdi=16,lower_fdi=46,gap_interval_mm=['0','0']),dict(upper_fdi=26,lower_fdi=36,gap_interval_mm=['0.1','0.1'])])
    c=dict(schema='research-crown-v1',**meta,representation='full_mesh',tooth_fdi=16,vertices_mm=v,faces=f,face_roles=roles,margin_z_mm=0.,wall_lower_mm=1.,relief_cap_mm=.1)
    return b,c

def x97():
    with package('X97','occlusion_module') as m:
        api=importlib.import_module('occlusion_module.api');reach=importlib.import_module('occlusion_module._vendor.reachability');contact=importlib.import_module('occlusion_module._vendor.contact')
        b,c=synthetic_bite();a=m.analyze(b,c,{'propose_adjustment':False})
        w=reach.sufficiency(contact)
        classes={r['fdi']:r['classification'] for r in a['force_intervals']['per_tooth']}
        return dict(contact_area_mm2=a['contact_map']['metrics']['predicted']['area_mm2'],force_classes=[classes[i] for i in [16,26,17]],missing_force_upper=a['force_intervals']['upper_bound_missing'],summary_identity_error=w['summary_identity_error'],spatial_difference_mm2=w['downstream_difference_mm2'],physical_validation=a['physical_validation'],result=a,resolution='PER_POINT/PER_SURFACE_REGION/PER_TOOTH',scope='Synthetic open surface carrier tests API composition; not a manufactured crown or original ten-crown replay')

def x98():
    with package('X98','implant_safety') as m, tempfile.TemporaryDirectory(prefix='dental-safety-') as tmp:
        root=Path(tmp);(root/'data').mkdir()
        for name,value in [('sites.json',[]),('guide_profiles.json',fixture('guide_publication.json')),('protocols.json',fixture('drill_protocols.json'))]:
            (root/'data'/name).write_text(json.dumps(value))
        (root/'SOURCE_MANIFEST.json').write_text('[]')
        s=m.SafetyModule(root);frame=dict(order='zyx',units='mm',spacing_mm=[1.,1.,1.],origin_mm=[0.,0.,0.],direction=np.eye(3).tolist())
        p=root/'synthetic.npz';np.savez(p,nominal_voxels_zyx=np.array([[4,0,6]],int),revision_voxels_zyx=np.array([[2,0,4]],int))
        s.register_site(site_id='synthetic',points_path=p,points_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),frame=frame,nominal_key='nominal_voxels_zyx',source_locator='synthetic occupancy; no anatomical measurement')
        pose=dict(entry_zyx_mm=[0.,0.,0.],axis_zyx=[0.,0.,1.],frame=frame)
        query=dict(site_id='synthetic',cbct_pose=pose,guide_type='fully_guided',implant_system='Straumann_BLX_2024_VeloDrill',length_mm=4.,radius_mm=1.,confidence=.95)
        a=s.query(**query);unknown=s.query(**dict(query,implant_system='missing'))
        bad=False
        try:s.query(**dict(query,cbct_pose=dict(pose,axis_zyx=[0.,0.,2.])))
        except m.ContractError:bad=True
        geo=a['digital_geometry'];nom=geo['nominal']['interval_mm'];union=geo['revision_union']['interval_mm']
        require(nom[0]<=math.sqrt(8.5)<=nom[1]+1e-12,'independent box-to-cylinder nominal distance')
        require(union[0]<=.5<=union[1]+1e-12,'independent union distance')
        return dict(nominal_distance_squared_mm2=8.5,nominal_interval_mm=nom,revision_union_gap_mm=union[1],max_bracket_width_mm=max(x['width_mm'] for x in geo.values() if isinstance(x,dict) and 'width_mm' in x),physical_safety=a['physical_safety']['status'],injury_probability=a['physical_safety']['injury_probability'],bad_pose_rejected=bad,missing_protocol=unknown['combined_conditional']['class_vs_2mm'],result=a,resolution='PER_TOOTH; guide assumptions PHENOMENOLOGICAL')

def point_triangle_distance(point,tri):
    # Independent finite-fixture control: projected interior plus three edges.
    a,b,c=tri;u=b-a;v=c-a;w=point-a
    gram=np.array([[u@u,u@v],[u@v,v@v]])
    st=np.linalg.solve(gram,[u@w,v@w]);candidates=[]
    if st.min()>=0 and st.sum()<=1:candidates.append(a+st[0]*u+st[1]*v)
    for start,end in [(a,b),(b,c),(c,a)]:
        edge=end-start;t=np.clip((point-start)@edge/(edge@edge),0,1)
        candidates.append(start+t*edge)
    return min(float(np.linalg.norm(point-q)) for q in candidates)

def insertion():
    names=['rational','nextdown','nextup','dot_interval','separator_batch','rref_solve','primal_matrix','validate_witness','lp_witness','independent_lp']
    ns=defs('B31_INSERTION','code/sweep.py',names,dict(DOWN=-np.inf,UP=np.inf))
    A,B,C,D=defs('B31_INSERTION','code/fixtures.py',['fixture'])['fixture']()
    sep,_=ns['separator_batch'](np.repeat(A[None],len(B),axis=0),B,D)
    witness=None;wtri=None
    for tri in C:
        w,_=ns['lp_witness'](A,tri,D)
        if w is not None:witness=w;wtri=tri;break
    require(witness is not None,'source continuous collision witness missing')
    wrong=copy.deepcopy(witness);wrong['weights'][0]='-1'
    safe=all(not ns['independent_lp'](A,t,D)['feasible'] for t in B)
    bad=any(ns['independent_lp'](A,t,D)['feasible'] for t in C)
    summaries=[[min(point_triangle_distance(p,t) for p in A+shift for t in triangles) for shift in [D*0,D]] for triangles in [B,C]]
    identity=float(abs(np.array(summaries[0])-summaries[1]).max())
    return dict(safe_separators=int(sep.sum()),collision_witness_valid=ns['validate_witness'](A,wtri,D,witness),wrong_weight_rejected=not ns['validate_witness'](A,wtri,D,wrong),control_matches=safe and bad,summary_identity_error=identity,endpoint_gap_summaries_mm=summaries,path_decisions_differ=bool(sep.all() and witness),witness=witness,resolution='PER_POINT',minimum_extension='Keep whole translation and spatial obstacle triangles, beyond identical endpoint gap',endpoint_scope='All fixture source vertices; identical planar endpoint region. No general continuous mesh clearance certificate.')

def support():
    ns=defs('X95_NATIVE_SUPPORT_SHAPE','code/port_certificate.py',['down','up','box_lower','vertex_upper'])
    lo=np.array([[3.,0,0]]);hi=np.array([[4.,1,1]]);point=np.zeros(3)
    lower=float(ns['box_lower'](point,lo,hi).min());upper=float(ns['vertex_upper'](point,lo).min())
    require(lower<=3<=upper,'outward box/vertex interval')
    p2=defs('X95_NATIVE_SUPPORT_SHAPE','code/native_port.py',['p2'])['p2']
    weights=p2(np.array([[.25,.25,.5]]))
    wronglo=lo.copy();wronglo[0,0]=3.1
    return dict(distance_enclosed_mm=3.,interval_mm=[lower,upper],shape_function_partition_error=float(abs(weights.sum()-1)),expanded_box_fault_rejected=not np.all(lo>=wronglo),rigorous_bound='OUTWARD_ROUNDED_STRAIGHT_BOX_VERTEX',resolution='PER_POINT',physical_support='UNKNOWN')

def pulp_port():
    with kernel('X95_CT_PULP_PREP_PORT','code/x12_consumer_ports.py') as m:
        a=m.uniform_shell_decision(2.,.5,1.,.5)
        bad=False
        try:m.uniform_shell_decision(-1.,.5,1.,.5)
        except ValueError:bad=True
        return dict(required_low_mm=a['distance_lower_mm'],decision=a['decision'],bad_distance_rejected=bad,clinical_feasibility=a['clinical_feasibility'],result=a,resolution='PER_POINT',debt='Annotation-derived total hard tissue is not a calibrated dentin measurement')

def seating():
    from mpmath import iv
    ns=defs('X95_SEATING_IDENTIFIABLE_ASSAY','assay.py',['interval','ik','integrate_nominal','segment_coverage'],dict(iv=iv,Unknown=ValueError))
    z=ns['ik'](0,1.5);nom=ns['integrate_nominal']([[0,1,1],[1,2,.5]],0,2,0,0)
    require(float(z.a)<=1.5<=float(z.b),'constant-age integral')
    missing=False
    try:ns['segment_coverage']([[0,1,1],[1.5,2,.5]],0,2)
    except ValueError:missing=True
    return dict(integral_exact='3/2',nominal_control_tolerance=abs(nom-1.5),gap_without_history='UNKNOWN' if missing else 'ACCEPTED',physical_status='UNKNOWN',resolution='PHENOMENOLOGICAL',scope='Conditional constant-age hydraulic dose primitive; full identifiable laboratory CSV port is included as source')

def contact_repair():
    with package('X97','occlusion_module'):
        reach=importlib.import_module('occlusion_module._vendor.reachability');contact=importlib.import_module('occlusion_module._vendor.contact');w=reach.sufficiency(contact)
    ns=defs('B31_CONTACT','code/reachability.py',['floor_compatible'])
    return dict(floor_fault_rejected=not ns['floor_compatible'](1.,{'attained':.5}),summary_identity_error=w['summary_identity_error'],spatial_difference_mm2=w['downstream_difference_mm2'],physical_force='UNKNOWN',resolution='PER_SURFACE_REGION',minimum_extension='Reference overlap for fixed score; spatial masks for changed pose or repair')

def layer_nesting():
    ns=defs('X95_LAYER_THICKNESS_NESTING','code/height_core.py',['largest_rectangle','base_footprint','cap_prefix','vertex_caps','formula_volume','sequential_control'])
    shape=(2,3,1);idx=np.arange(6);allowed=np.ones(6,bool);occ=np.ones(shape,bool);base=ns['base_footprint'](occ)
    caps=ns['vertex_caps'](ns['cap_prefix'](allowed,shape,idx,base));control=ns['sequential_control'](allowed,shape,idx,base)
    a=np.array([[1,1,1],[1,1,1],[0,0,0]],bool);b=np.array([[1,1,0],[1,0,1],[0,1,1]],bool)
    areas=[ns['largest_rectangle'](x)[0] for x in [a,b]]
    return dict(caps_match=np.array_equal(caps,control),volume_mm3=ns['formula_volume'](caps,.5),max_volume_error_mm3=abs(ns['formula_volume'](caps,.5)-.75),rectangle_area=areas[0],summary_identity_error=float(abs(int(a.sum())-int(b.sum()))),changed_area=areas[0]!=areas[1],physical_manufacture='UNKNOWN',resolution='PER_POINT/PER_SURFACE_REGION',minimum_extension='Spatial occupancy and oriented layer/pose, beyond identical occupancy count')

RUNNERS={'CAD':cad,'CROWN_PREP':crown_prep,'X87':x87,'X88':x88,'X89':x89,'X96':x96,'X97':x97,'X98':x98,'INSERTION':insertion,'SUPPORT':support,'PULP_PORT':pulp_port,'SEATING':seating,'CONTACT_REPAIR':contact_repair,'LAYER_NESTING':layer_nesting}
MUTATIONS={'CAD':('coordinate_error_mm',1.),'CROWN_PREP':('valid_normals',False),'X87':('zero_residual_equality','ABOVE'),'X88':('counts',[[0,39],[10,27],[9,19]]),'X89':('upper_rational','1/1000'),'X96':('full_radius_gap_mm',2.5),'X97':('contact_area_mm2',0.),'X98':('physical_safety','PASS'),'INSERTION':('collision_witness_valid',False),'SUPPORT':('distance_enclosed_mm',1.),'PULP_PORT':('decision','PASS'),'SEATING':('integral_exact','1'),'CONTACT_REPAIR':('floor_fault_rejected',False),'LAYER_NESTING':('volume_mm3',1.)}
TOLERANCES={'coordinate_error_mm':0,'summary_identity_error':0,'waist_rebound_difference_mm':0,'max_control_error_mm':1e-10,'max_interval_control_error':1e-10,'max_quadrature_error_mm_N':1e-12,'max_bracket_width_mm':1e-6,'shape_function_partition_error':1e-15,'nominal_control_tolerance':1e-14,'max_volume_error_mm3':1e-12}

def check(name,result,contract):
    for key,value in contract.items():
        require(key in result,name+': missing '+key)
        actual=result[key]
        if key in TOLERANCES and key.startswith(('max_','shape_function','nominal_control')):
            require(isinstance(actual,(int,float)) and math.isfinite(actual) and 0<=actual<=value,name+': '+key)
        elif isinstance(value,float):require(isinstance(actual,(int,float)) and math.isfinite(actual) and abs(actual-value)<=1e-10,name+': '+key)
        else:require(actual==value,name+': '+key)

def verify_freeze():
    f=json.loads((ROOT/'FROZEN_CAPABILITIES.json').read_text())
    for name,h in f['files'].items():require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,'Capability hash drift: '+name)

def clean(value):
    if isinstance(value,np.ndarray):return clean(value.tolist())
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(x) for x in value]
    return value

def figure(rows,out):
    if 'X96' not in rows:return
    r=rows['X96'];parts=['<svg xmlns="http://www.w3.org/2000/svg" width="720" height="260" viewBox="0 0 720 260">','<rect width="720" height="260" fill="white"/>','<g font-family="sans-serif" font-size="15" fill="#192c3a">','<text x="24" y="30">Identical nominal gap and added depth; changed swept shape</text>']
    for i,(label,key) in enumerate([('Original cylinder','nominal_gap_mm'),('Full-radius extension','full_radius_gap_mm'),('Axis-only extension','minimal_tool_gap_mm')]):
        y=58+48*i;v=r[key];parts.extend([f'<text x="24" y="{y+19}">{label}</text>',f'<rect x="230" y="{y}" width="{v*120}" height="28" fill="#397c91"/>',f'<text x="{242+v*120}" y="{y+19}">{v:g} mm</text>'])
    parts.extend(['<text x="24" y="223">Synthetic Euclidean fixture, PER_POINT. Physical safety: UNKNOWN.</text>','</g></svg>'])
    (out/'tool_shape_witness.svg').write_text('\n'.join(parts)+'\n')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('profiles',nargs='*');p.add_argument('--list',action='store_true');p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.list:
        print('\n'.join(RUNNERS));return
    selected=a.profiles or list(RUNNERS)
    if any(x not in RUNNERS for x in selected):p.error('Unknown capability; use --list')
    verify_freeze();contracts=fixture('capability_contracts.json')['profiles'];rows={}
    for name in selected:
        tick=time.perf_counter();r=clean(RUNNERS[name]());check(name,r,contracts[name])
        bad=copy.deepcopy(r);key,value=MUTATIONS[name];bad[key]=value;rejected=False
        try:check(name,bad,contracts[name])
        except AssertionError:rejected=True
        require(rejected,'Wrong result not rejected: '+name)
        rows[name]=dict(r,status='PASS',injected_wrong_result_rejected=rejected,elapsed_s=time.perf_counter()-tick)
        print(name+' PASS; injected wrong result rejected',flush=True)
    out=a.output or Path(tempfile.mkdtemp(prefix='dental-capabilities-'));out.mkdir(parents=True,exist_ok=True)
    (out/'capability_results.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    figure(rows,out)
    print(str(out))

if __name__=='__main__':main()
