"""Tests mathematical contracts and preserved empirical failures; no biological admission."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
import unittest,json,copy,hashlib,time,argparse
from pathlib import Path
import numpy as np
from surgical_chain.ports import Port,synth,unknown,UnknownPort
from surgical_chain.hemostasis import hydraulic,bleeding,CO,PL
from surgical_chain.inventories import query_strength,JointInventory,reservoir
from surgical_chain.chain import run,resume,incision,check_config
from surgical_chain.scenario import scenario_config,evidence_config
from surgical_chain.vendor import response_r1 as r
R=Path(__file__).resolve().parent
S=R/'sources'
TEST_OUT=None

class PortTests(unittest.TestCase):
 def test_unknown_never_zero(self):
  with self.assertRaises(UnknownPort):unknown('m','src','Unknown injury').require('m')
  with self.assertRaises(ValueError):Port(0.,'m',"UNKNOWN",None,'src','scope')
 def test_null_needs_unknown(self):
  with self.assertRaises(ValueError):Port(None,'m',"MEASURED",None,'src','scope')
 def test_units(self):
  with self.assertRaises(ValueError):synth(1.,'mm','src','length').require('m')
 def test_finite(self):
  with self.assertRaises(ValueError):synth(float('nan'),'m','src','length')
 def test_provenance(self):
  with self.assertRaises(ValueError):Port(1.,'m',"MEASURED",None,'','')
 def test_evidence_stops(self):
  d=run(evidence_config());self.assertIsNone(d['native_strength']['value'])
  self.assertTrue(all(x['value'] is None for x in d['stages'].values()))
 def test_synthetic_cannot_be_evidence(self):
  self.assertIsNone(run(scenario_config())['native_strength']['value'])
 def test_partial_scenario_stops(self):
  d=run(evidence_config(),mode='scenario');self.assertIsNone(d['native_strength']['value'])
 def test_unit_error_not_hidden_by_missing_data(self):
  d=evidence_config();d['parameters']['gap']['unit']='J/m²'
  with self.assertRaises(ValueError):run(d)
 def test_path_work_and_area(self):
  d=scenario_config();m=incision(d,check_config(d));self.assertAlmostEqual(m['area']['value'],2e-5)
  self.assertAlmostEqual(m['work']['value'],.003);self.assertAlmostEqual(m['force']['value'],.3)
  self.assertIsNone(m['total_tool_work']['value'])
 def test_zero_length_work(self):
  d=scenario_config();d['path']['value']=[[0.,0.,0.],[0.,0.,0.]]
  m=incision(d,check_config(d));self.assertEqual(m['work']['value'],0);self.assertEqual(m['area']['value'],0)
 def test_stack_depth(self):
  d=scenario_config();d['cut_depth']['value']=.01
  with self.assertRaises(ValueError):incision(d,check_config(d))

class PhysicsTests(unittest.TestCase):
 def hist(self,radius=8e-6,density=1e6,pressure=4000.,area=1e-5):return hydraulic([radius],[density],area,pressure,.0035,.001)
 def test_poiseuille_radius4_shear_radius1(self):
  a=self.hist();b=self.hist(radius=16e-6)
  self.assertAlmostEqual(b['initial_flow_m3_s']/a['initial_flow_m3_s'],16.)
  self.assertAlmostEqual(b['gamma_per_bin_s_inverse'][0]/a['gamma_per_bin_s_inverse'][0],2.)
 def test_zero_vessels_and_pressure(self):
  self.assertEqual(self.hist(density=0)['initial_flow_m3_s'],0.)
  self.assertEqual(self.hist(pressure=0)['initial_flow_m3_s'],0.)
  self.assertEqual(self.hist(area=0)['initial_flow_m3_s'],0.)
 def test_negative_histogram(self):
  with self.assertRaises(ValueError):self.hist(density=-1)
 def test_hemostasis_without_bleeding(self):
  b=bleeding(self.hist(density=0),minutes=10)
  self.assertAlmostEqual(b['summary']['blood_volume_uL'],0.)
 def test_fibrin_mass(self):
  b=bleeding(self.hist(),minutes=10);self.assertLess(b['summary']['max_fibrinogen_inventory_error_nM'],1e-6)
 def test_platelet_aggregation_block(self):
  # Low-shear bin permits adhesion; no-aggregation cannot reach threshold(.5)
  # since adhesion weight(.3) alone is insufficient.
  b=bleeding(self.hist(radius=3e-6),minutes=10,aggregation_factor=0.)
  self.assertIsNone(b['summary']['hemostasis_min']);self.assertLess(np.max(b['coverage']),PL['THETA_C'])
  self.assertEqual(float(b['seal'].max()),0.)
 def test_inventory_zero_birth(self):
  a=JointInventory(3)
  for i in range(20):a.step(.1,np.zeros(3),np.ones(3))
  self.assertEqual(float(abs(a.q).max()),0.)
 def test_inventory_conservation_and_PSD(self):
  a=JointInventory(3)
  for i in range(200):a.step(.1,np.array([.1,.2,.3]),np.array([0,.5,1]))
  self.assertLess(a.maximum_balance,1e-9)
  self.assertGreaterEqual(float(np.linalg.eigvalsh(np.moveaxis(a.q,-1,0)).min()),-1e-10)
 def test_no_oxygen_no_mature(self):
  a=JointInventory(1)
  for i in range(100):a.step(.1,np.ones(1)*.1,np.zeros(1))
  self.assertEqual(float(query_strength(a.q[2],0)[0]),0.)
 def test_joint_marginals(self):
  self.assertAlmostEqual(query_strength(np.diag([.1,.4]),0),15.)
  self.assertAlmostEqual(query_strength(np.diag([.4,.1]),0),60.)
 def test_rotation_and_no_bridges(self):
  q=np.diag([.4,.1]);self.assertAlmostEqual(query_strength(q,90),15.)
  self.assertEqual(query_strength(q,0,0),0.)
 def test_invalid_tensor(self):
  with self.assertRaises(ValueError):query_strength(np.diag([-.1,.4]),0)
 def test_healthy_oxygen_limit(self):
  p=r.Params(face_conductance=0.);g=r.Grid(p);o=r.Oxygen(g,p)
  v=o.solve(np.ones(g.n),np.full(g.n,p.consumption))
  self.assertLess(float(abs(v-p.healthy_Torr).max()),1e-7)
  self.assertLess(o.max_balance,1e-6)
 def test_reservoir_stoichiometry_and_restart(self):
  for law in ['first_order','second_order']:
   a=reservoir(.066,.018,0.,np.arange(81),.0695,.283333,law)
   self.assertLess(float(abs(a[0]+2*a[1]+a[2]-(.066+2*.018)).max()),1e-12)
   mid=reservoir(.066,.018,0.,11.,.0695,.283333,law)
   end=reservoir(*mid,69.,.0695,.283333,law)
   np.testing.assert_allclose(end,a[:,-1],atol=1e-12)

class SourceRegressionTests(unittest.TestCase):
 def load(self,lane,file):return json.loads((S/lane/file).read_text())
 def test_frozen_sources(self):
  m=json.loads((R/'SOURCE_MANIFEST_R1.json').read_text())
  for f in m['files']:self.assertEqual(hashlib.sha256((R/f['snapshot']).read_bytes()).hexdigest(),f['sha256'])
 def test_prereg_frozen(self):
  for name in ['PREREG_R1','PREREG_R1_SEAL_MEMORY','PREREG_R1_RESPONSE_R4','PREREG_R1_REFINEMENT']:
   self.assertEqual(hashlib.sha256((R/(name+'.json')).read_bytes()).hexdigest(),(R/(name+'.sha256')).read_text().strip())
 def test_needle_new_tools_fail_speed_pass(self):
  d=self.load('LANE_SURGICAL_INCISION','r4/fracture_transfer.json');errors=[]
  for row in d['heldout_tools']:
   err=abs(row['predicted_effective_J_J_m2']/row['heldout_effective_J_J_m2']-1)
   self.assertAlmostEqual(err,row['J_relative_error'],12);self.assertGreater(err,.2);errors.append(err)
  np.testing.assert_allclose(errors,[.3595505617977528,.2214285714285715,.5464190981432361],atol=1e-12)
  err=max(abs(a['predicted_J_J_m2']/a['observed_J_J_m2']-1) for a in d['heldout_speed_same_tool'])
  self.assertAlmostEqual(err,.028409090909090884,12)
 def test_human_needle_heldout(self):
  d=self.load('LANE_SURGICAL_INCISION','r3/geometry/needle_geometry_results.json')
  err=abs(d['heldout_predicted_force_N']/d['heldout_observed_force_N']-1)
  self.assertAlmostEqual(err,.11428571428571425,12);self.assertAlmostEqual(d['extreme_digitization_relative_error'],.38333333333333347,12)
 def test_exact_prospective_quotient_bound(self):
  bound=(8+5+5+3500*.002)/(1-.002)+4
  self.assertAlmostEqual(bound,29.050100200400802,10)
 def test_micro_interface_ceiling_is_conditional(self):
  d=self.load('LANE_SURGICAL_BINDINGS','INTERFACE_BUDGET_RESULTS_R4.json');x=d['cylindrical_inventory_sensitivity'][0]
  self.assertAlmostEqual(x['Gamma_cut_upper_scenario_J_m2']+x['M_max']*x['analogue_gamma_upper_J_m2'],6965.365853658535,7)
  self.assertIsNone(d['actual_native_dermal_prediction_J_m2'])
 def test_gap_bridge_bounds_and_mode_failure(self):
  d=self.load('LANE_SKIN_TOUGHNESS_GAP','MECHANISM_TABLE_R2_FINAL.json');tab={x['mechanism']:x for x in d['table']}
  self.assertLess(tab['Single crimp straightening']['predicted_extra_kJ_m2'][1],.14)
  h=tab['Finite-volume hydraulic ramp plus complete drainage']['predicted_extra_kJ_m2']['generous_h10mm_range'][1]
  self.assertAlmostEqual(h,.11046636166000416,12)
  x=self.load('LANE_SKIN_TOUGHNESS_GAP','r2/rotation_vector_v2/summary.json')
  self.assertEqual(x['predicted_fixed_inventory_I_III_ratio'],1.)
  self.assertAlmostEqual(x['observed_all_I_III_ratio'],30380/20600,12)
  self.assertEqual(x['native_joint_gate'],'UNKNOWN')
 def test_R3_strength_from_raw_curve(self):
  d=self.load('LANE_SURGICAL_RESPONSE','r3/crosslink_port/summary.json');native=next(x for x in d['curves'] if x['label']=='native')['strength']
  raw=np.load(S/'LANE_SURGICAL_RESPONSE/r3/crosslink_port/native.npz')
  held=[x for x in native['rows'] if x['independent_day']]
  pred=[float(np.interp(x['day'],raw['days'],raw['strength_pct'])) for x in held]
  rmse=float(np.sqrt(np.mean([(p-x['target_pct'])**2 for p,x in zip(pred,held)])))
  self.assertAlmostEqual(rmse,8.306367094671637,10);self.assertEqual(native['heldout_gate'],'FAIL')
 def test_R3_HP_input_not_early_prediction(self):
  d=self.load('LANE_SURGICAL_RESPONSE','RESPONSE_PORTS_R3.json')
  self.assertIsNone(d['conditional_operator']['native_stress_law']);self.assertIsNone(d['observed_inputs']['joint_biological_parameter_posterior'])
 def test_swarm_count_direction_and_wrong_text(self):
  d=self.load('BT-FW48-AUTO-55a5902a27c3f6','vessel_map_flow_results.json')
  radii=[np.sqrt(3*10)*1e-6,np.sqrt(10*30)*1e-6,np.sqrt(30*100)*1e-6,np.sqrt(100*300)*1e-6]
  n=np.array([15000,15000,300,300]);q=np.array(radii)**4
  ratio=(np.full(4,n.sum()/4)@q)/(n@q)
  self.assertAlmostEqual(ratio,25.375671641791048,10);self.assertGreater(ratio,1.)
  self.assertIn('R^3',d['equations']['gamma']) # preserved text error, never used as physics


 def test_R4_inventory_denominator_falsifier(self):
  o=self.load('LANE_SURGICAL_RESPONSE','OBSERVATIONS_R4.json');d=self.load('LANE_SURGICAL_RESPONSE','r4/SCORES_R4_v2.json')
  ratio=o['HP42_mol_mol']*o['collagen42_ug']/(o['HP_early_mol_mol'][-1]*o['collagen_early_ug'][-1])
  self.assertAlmostEqual(ratio,2.220261437908497,12)
  self.assertEqual(d['all_profile_joint_pass_count'],0)
  self.assertIn('FAIL',d['inventory']['selective_turnover_only_gate'])
 def test_R4_external_rates_and_native_null(self):
  pre=self.load('LANE_SURGICAL_RESPONSE','PREREG_R4_GLYCO_RESERVOIR.json');p=self.load('LANE_SURGICAL_RESPONSE','RESPONSE_PORTS_R4.json')
  for row in p['independent_analog_rates']:
   d=pre['source']['reported_average_endpoints'][row['group']]
   eta=2*(d['H28']-d['H0'])/(d['D0']-d['D28'])
   k=np.log(d['D0']/d['D28'])/28 if row['law']=='first_order' else (1/d['D28']-1/d['D0'])/28
   np.testing.assert_allclose(reservoir(d['D0'],d['H0'],0,28,k,eta,row['law']),row['source_endpoint_reconstruction'],atol=1e-12)
  self.assertTrue(all(x is None for x in p['native'].values()))
 def test_R4_nominal_and_all_glyco_FAIL(self):
  d=self.load('LANE_SURGICAL_RESPONSE','r4/glyco_reservoir/summary.json');nom=next(x for x in d['rows'] if x['nominal'])
  self.assertEqual(len(d['rows']),16);self.assertTrue(all(x['nominal_joint_gate']=='FAIL' for x in d['rows']))
  self.assertAlmostEqual(nom['prediction_HP']['42'],.026338454833019345,12)
  self.assertAlmostEqual(nom['strength']['heldout_RMSE_pp'],23.17981084535914,9)
 def test_R4_source_freeze(self):
  m=json.loads((R/'SOURCE_MANIFEST_R1_RESPONSE_R4.json').read_text())
  for f in m['files']:self.assertEqual(hashlib.sha256((R/f['snapshot']).read_bytes()).hexdigest(),f['sha256'])

class ExecutedChainTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.out=R/'r1/nominal_complete_v4';cls.d=json.loads((cls.out/'summary.json').read_text());cls.a=np.load(cls.out/'fields.npz')
 def test_chain_unknown_native_not_promoted(self):
  self.assertIsNone(self.d['native_rupture_strength']['value']);self.assertEqual(self.d['empirical_joint_gate'],'UNKNOWN')
  for key in ['oxygen_edge_pressure','hypoxic_width_proxy','inventory_UIM_final','strength_time']:self.assertEqual(self.d[key]['status'],'SYNTETISKT')
 def test_chain_mass_and_oxygen(self):
  n=self.d['numerics'];self.assertLess(n['UIM_balance'],1e-9)
  self.assertLess(n['oxygen_inventory']['relative_conservation_error'],1e-6)
  self.assertLess(n['early_oxygen']['max_relative_balance'],1e-6);self.assertLess(n['late_oxygen']['max_relative_balance'],1e-6)
 def test_persist_seal_after_thrombin(self):
  s=self.a['seal'];self.assertGreater(s[-1],0);self.assertGreaterEqual(float(np.diff(s).min()),-1e-9)
 def test_joint_mass_nonnegative(self):
  self.assertGreaterEqual(float(self.a['UIM'].min()),-1e-10)
 def test_finite_every_raw_field(self):
  for key in self.a.files:self.assertTrue(np.isfinite(self.a[key]).all(),key)
 def test_checkpoint_full_suffix(self):
  dest=TEST_OUT/'restart';late=resume(self.out/'checkpoint21.npz',dest)
  mask=self.a['days']>21+1e-8
  np.testing.assert_allclose(late['days'],self.a['days'][mask],atol=1e-9)
  self.assertLess(float(abs(late['strength_pct']-self.a['strength_pct'][mask]).max()),1e-7)
  np.testing.assert_allclose(late['Q_final'],self.a['Q_final'],atol=1e-9)


 def test_chemistry_ledger_and_checkpoint(self):
  a=self.a['chemistry_D_H_A'];self.assertLess(float(abs(a[0]+2*a[1]+a[2]-(.066+2*.018)).max()),1e-10)
  h42=float(np.interp(42,self.a['chemistry_days'],a[1]));self.assertAlmostEqual(h42,.026338454833019345,12)
 def test_explicit_inherited_closures(self):
  self.assertEqual(len(self.d['inherited_closures']['Q036_parameters']),38)
  self.assertTrue(all(x['status']=='SYNTETISKT' for x in self.d['inherited_closures']['coagulation_rates'].values()))

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
 TEST_OUT=args.out.resolve()
 if not TEST_OUT.is_relative_to(R):raise ValueError('Tests may write only in the synthesis lane')
 TEST_OUT.mkdir(parents=True,exist_ok=False)
 suite=unittest.defaultTestLoader.loadTestsFromModule(__import__(__name__));st=time.perf_counter()
 result=unittest.TextTestRunner(verbosity=2).run(suite)
 with (TEST_OUT/'verification.json').open('x') as f:json.dump({'review_state':'PENDING_INDEPENDENT_REVIEW','passed':result.testsRun-len(result.failures)-len(result.errors),'count':result.testsRun,'failures':[{'test':str(t),'trace':s} for t,s in result.failures+result.errors],'gate':'PASS' if result.wasSuccessful() else 'FAIL','wall_s':time.perf_counter()-st,'scope':'Model/port/numeric/source regression, not empirical biological validation'},f,indent=2)
 raise SystemExit(0 if result.wasSuccessful() else 1)
