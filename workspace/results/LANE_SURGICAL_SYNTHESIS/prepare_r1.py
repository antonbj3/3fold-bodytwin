from pathlib import Path
import hashlib,json,shutil,datetime
ROOT=Path(__file__).resolve().parent
W=ROOT.parent.parent

def write(p,x):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x') as f:json.dump(x,f,indent=2,ensure_ascii=False,allow_nan=False)

def snap(src,dst):
 dst.parent.mkdir(parents=True,exist_ok=True)
 if dst.exists():raise FileExistsError(dst)
 shutil.copyfile(src,dst)
 return {'source':str(src),'snapshot':str(dst.relative_to(ROOT)),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'bytes':src.stat().st_size}

files=[]
for lane,n in [('LANE_SURGICAL_INCISION',4),('LANE_SURGICAL_BINDINGS',4),('LANE_SKIN_TOUGHNESS_GAP',2),('LANE_SURGICAL_RESPONSE',3)]:
 src=W/'results'/lane
 names=['RESULTS.md','WORK_STATUS.json','NEXT_ROUND.md',f'VERIFICATION_R{n}.json',f'MEASUREMENT_SPEC_R{n}.md']
 names+=['PORTS_R2_FINAL_V2.json','MECHANISM_TABLE_R2_FINAL.json'] if n==2 else ['RESPONSE_PORTS_R3.json'] if n==3 else ['PORTS_R4_FINAL.json']
 for name in names:
  if (src/name).exists():files.append(snap(src/name,ROOT/'sources'/lane/name))
 for p in sorted((src/'night_rounds').glob('*.json'))[-1:]:files.append(snap(p,ROOT/'sources'/lane/'night_rounds'/p.name))
# Full reference reports are preserved, never treated as a common posterior.
for lane,name in [('LANE_SURGICAL_INCISION','r4/fracture_transfer.json'),('LANE_SURGICAL_INCISION','r4/summary.json'),('LANE_SURGICAL_INCISION','r3/geometry/needle_geometry_results.json'),('LANE_SURGICAL_BINDINGS','INTERFACE_BUDGET_RESULTS_R4.json'),('LANE_SKIN_TOUGHNESS_GAP','r2/rotation_vector_v2/summary.json'),('LANE_SURGICAL_RESPONSE','r3/crosslink_port/summary.json'),('LANE_SURGICAL_RESPONSE','r3/crosslink_port/native.npz'),('LANE_SURGICAL_RESPONSE','r3/kinetics_uncertainty/summary.json'),('LANE_SURGICAL_RESPONSE','r3/joint_history/summary.json'),('LANE_SURGICAL_RESPONSE','r1/coupled_v2/calibration.json'),('LANE_SURGICAL_RESPONSE','r2/world/summary.json'),('LANE_SURGICAL_RESPONSE','PREREG_R3_CROSSLINK_PORT.json'),('LANE_SURGICAL_RESPONSE','WORLD_TARGETS_R3_v2.json')]:
 files.append(snap(W/'results'/lane/name,ROOT/'sources'/lane/name))
files.append(snap(W/'tasks/build_night/NIGHT_LOG.md',ROOT/'sources/NIGHT_LOG_initial.md'))
code=ROOT/'surgical_chain/vendor';code.mkdir(parents=True,exist_ok=False)
for name in ['response_r1.py','storage_r1.py']:
 files.append(snap(W/'results/LANE_SURGICAL_RESPONSE'/name,code/(name+'.original')))
 text=(code/(name+'.original')).read_text()
 if name=='response_r1.py':
  text=text.replace("QPATH = ROOT.parent / 'BT-HX-Q036' / 'model.py'","QPATH = ROOT / 'q036.py'")
 else:text=text.replace('import response_r1 as r','from . import response_r1 as r')
 (code/name).write_text(text)
files.append(snap(W/'results/BT-HX-Q036/model.py',code/'q036.py'))
# Never import executable MSK scripts: use their AST definitions before their runners.
for name in ['coagulation_hemostasis.py','platelet_hemostasis.py']:
 files.append(snap(Path('source_repository/scripts/msk')/name,code/(name+'.source')))
for name in ['__init__.py']:(code/name).write_text('')
for job,names in [('55a5902a27c3f6',['vessel_map_flow.py','vessel_map_flow_results.json','RESULTS.md']),('27a6dac22a7a23',['platelet_hemostasis_shear.py','run_out/flow_coupled.json','RESULTS.md'])]:
 base=Path('/mnt/games-240/research/bunny48_20260926/bodytwin')/('BT-FW48-AUTO-'+job)
 for name in names:files.append(snap(base/name,ROOT/'sources'/('BT-FW48-AUTO-'+job)/name))
write(ROOT/'SOURCE_MANIFEST_R1.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':files,'scope':'Frozen local research and synthetic closures; no data sent externally'})
pre={'round':1,'lane':'LANE_SURGICAL_SYNTHESIS','review_state':'PENDING_INDEPENDENT_REVIEW','authority':'Latest synthesis steer: runnable tested chain and unified measurements','capability':'Tool+path+layer stack -> injury/gap -> radius-resolved bleeding+seal -> oxygen -> U/I/M directional bridge inventories -> strength over time','operation':'Strict dimensioned provenance ports plus explicit conditional executable chain; reuse frozen FV/Q036/coagulation equations','controls':['Same-information conventional FV/cohort/Poiseuille/coagulation operators: expect TIE','Separately calibrated observable curves: in-sample interpolation only; no heldout extrapolation declared without data','Count-only vessel histogram as information-loss diagnostic, never strongest control'],'unknowns':'Never zero-fill native cut, gap, viability/perfusion, chemistry or rupture law; default evidence mode stops at missing ports','defaults':'Explicit synthetic scenario, 10mm by2mm cut, dermis1.5mm/subcutis0.5mm; width250um,gap200um,5um front grid; 120min early/90day late; threadlimit2','measurements':['units/provenance/uncertainty per consumed and output port','regress source negative/positive outcomes with independent arithmetic/raw curves','energy and oxygen inventories; nonnegative mass and PSD tensors','r^4 and linear wall-shear analytic limits; zero wound, zero vessels, no chemistry','causal seals must not restore cut perfusion; unknown propagation','midpoint inventory refinement and early/late timestep refinement','full-replay vs checkpoint continuation with all required state','shared scenario draws stay joint; do not fabricate probabilities'],'error_requirements':{'numeric_mass':1e-9,'tensor_eigenvalue':-1e-10,'oxygen_relative_balance':1e-6,'source_numeric_regression':1e-7,'checkpoint_strength_pp':1e-7,'refinement_strength_pp':1.,'refinement_pressure_Torr':.5},'falsifiers':['Missing native port yields finite empirical strength','Mechanical width silently becomes ischemic width','Gamma becomes tensile strength','Lost radius or maturity-angle correlation','Unaccounted inventory or duplicate interface work','Claimed advantage against equally informed identical operator'],'scientific_gate':'UNKNOWN unless matched native acquisition; engineering PASS distinct from innovation TIE','costs':'Report setup/vendor import, full-chain and refinement/replay solves, storage/test cost separately; unknown laboratory acquisition cost=null'}
write(ROOT/'PREREG_R1.json',pre)
(ROOT/'PREREG_R1.sha256').write_text(hashlib.sha256((ROOT/'PREREG_R1.json').read_bytes()).hexdigest()+'\n')
print('Frozen',len(files),'source artifacts; PREREG_R1 written before chain calculations')
