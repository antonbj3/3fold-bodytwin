import copy,json,unittest
from pathlib import Path
from fractions import Fraction as Q
from surgical_chain.path_preference import (
    path_preference,fixture_preference,front_preference,incision_path_gate,
    verify_integrity,evaluate_request)

class ScalarTests(unittest.TestCase):
    def call(self,**kw):
        args=dict(gp=['1','1'],gd=['1/2','3/5'],gamma_layer=['1','1'],
            gamma_interface=['1/10','1/5'],mode='scenario',source='fixture',
            model_context='same load, geometry, extension')
        args.update(kw);return path_preference(**args)
    def test_labels(self):
        self.assertEqual(self.call()['decision'],"DEFLECTION")
        self.assertEqual(self.call(gamma_interface=['1','2'])['decision'],"REVIEW")
        self.assertEqual(self.call(gamma_interface=['1/2','3/5'])['decision'],'UNCERTAIN')
    def test_strict_boundary(self):
        self.assertEqual(self.call(gd=['1/2']*2,gamma_interface=['1/2']*2)['decision'],'UNCERTAIN')
        eps=Q(1,2**100)
        for sign,label in [(-1,"DEFLECTION"),(1,"REVIEW")]:
            self.assertEqual(self.call(gd=['1/2']*2,gamma_interface=[str(Q(1,2)+sign*eps)]*2)['decision'],label)
    def test_null_zero_context(self):
        for kw in (dict(gd=None),dict(model_context=None),dict(gp=['0','1']),dict(gamma_layer=['0','1'])):
            self.assertEqual(self.call(**kw)['decision'],'UNCERTAIN')
    def test_evidence(self):
        r=self.call(mode='evidence');self.assertEqual(r['decision'],'UNCERTAIN')
        self.assertIsNone(r['path_preference']['value'])
    def test_invalid(self):
        for kw in (dict(unit='J'),dict(gp=[1.,2.]),dict(gp=[True,2]),dict(gd=['-1','1']),
                   dict(gd=['2','1']),dict(mode='clinical'),dict(gp=['NaN','1'])):
            with self.assertRaises((ValueError,TypeError,ZeroDivisionError)):self.call(**kw)
    def test_no_unlicensed_outputs(self):
        r=self.call();self.assertFalse(r['physical_admission'])
        self.assertIsNone(r['predicted_depth']['value']);self.assertIsNone(r['biological_damage_width']['value'])
    def test_integrity_and_tamper(self):
        r=self.call();self.assertTrue(verify_integrity(r))
        r['decision']="REVIEW";self.assertFalse(verify_integrity(r))

class CornerRegressionTests(unittest.TestCase):
    def test_all_13_false_certain_corners(self):
        rows=json.loads(Path(__file__).with_name('corner_regressions_13.json').read_text())
        self.assertEqual(len(rows),13)
        for row in rows:
            with self.subTest(id=row['id']):
                result=evaluate_request(row['request'],mode='scenario')
                self.assertEqual(result['decision'],'UNCERTAIN')
                self.assertTrue(verify_integrity(result))
                self.assertNotEqual(row['corner_status'],row['reference_status'])
                self.assertEqual(result['decision'],row['expected']['status'])
    def test_fixture_evidence(self):
        row=json.loads(Path(__file__).with_name('corner_regressions_13.json').read_text())[0]
        r=evaluate_request(row['request']);self.assertEqual(r['decision'],'UNCERTAIN')
        self.assertIn('fixture_not_empirical_skin_provider',r['missing'])
    def test_angular_interior(self):
        args=dict(anchor='1',closure='SYNTHETIC_SIN2',mode='scenario')
        b=dict(r=['99/100','101/100'],gamma=['14/25','57/100'],theta=['-90','90'])
        full=fixture_preference(b,**args)
        self.assertEqual(full['decision'],'UNCERTAIN')
        for theta in ('-90','90'):
            corner=fixture_preference(dict(b,theta=[theta,theta]),**args)
            self.assertEqual(corner['decision'],"DEFLECTION")
        mid=fixture_preference(dict(b,theta=['0','0']),**args)
        self.assertEqual(mid['decision'],"REVIEW")
    def test_external_infinitesimal_anchor_is_separate(self):
        # Published homogeneous static mode III: R²=1/3. Synthetic sin² law.
        self.assertLess(Q(19,25)**2,Q(4,3))  # corner W=2
        self.assertGreater(Q(3,4)**2,Q(1,3)) # interior W=1
    def test_wrong_anchor_bad_depth(self):
        b=dict(r=['1','1'],gamma=['1','1'],theta=['0','0'])
        for kw in (dict(anchor='7'),dict(max_depth=True),dict(max_depth=11),dict(closure='inferred')):
            args=dict(anchor='1',closure='ISOTROPIC');args.update(kw)
            with self.assertRaises((ValueError,TypeError)):fixture_preference(b,**args)

class FrontTests(unittest.TestCase):
    def cells(self,gd=('3/2','3/2')):
        return [dict(cell_id=str(i),area_m2='1/2',gp=['2','2'],gd=[d,d],
            gamma_layer=['1','1'],gamma_interface=['1/2','1/2']) for i,d in enumerate(gd)]
    def call(self,cells=None,**kw):
        args=dict(front_complete=True,independent_strips=True,mode='scenario',
            source='parallel elastic strips',model_context='two local branches')
        args.update(kw);return front_preference(self.cells() if cells is None else cells,**args)
    def test_shape_changes_local_decision(self):
        a=self.cells();b=self.cells(('1/2','5/2'))
        for key in ('gp','gd','gamma_layer','gamma_interface'):
            self.assertEqual(sum(Q(c[key][0]) for c in a),sum(Q(c[key][0]) for c in b))
        self.assertEqual(self.call(a)['decision'],"DEFLECTION")
        r=self.call(b);self.assertEqual(r['decision'],'UNCERTAIN');self.assertTrue(r['mixed_front'])
        self.assertEqual(r['area_fractions'],{"REVIEW":'1/2',"DEFLECTION":'1/2','UNCERTAIN':'0'})
    def test_unknown_coupling_and_coverage(self):
        for kw in (dict(independent_strips=False),dict(front_complete=False)):
            self.assertEqual(self.call(**kw)['decision'],'UNCERTAIN')
    def test_bounded_coupling_and_boundary(self):
        self.assertEqual(self.call(independent_strips=False,coupling_margin_bound='1/4')['decision'],"DEFLECTION")
        self.assertEqual(self.call(independent_strips=False,coupling_margin_bound='1/2')['decision'],'UNCERTAIN')
    def test_refinement_preserves_parent_enclosure(self):
        parent={k:self.cells()[0][k] for k in ('gp','gd','gamma_layer','gamma_interface')}
        parent['gd']=['1/2','5/2']
        self.assertEqual(self.call(parent_enclosures=parent)['decision'],"DEFLECTION")
        parent['gd']=['1/2','1'];
        with self.assertRaises(ValueError):self.call(parent_enclosures=parent)
    def test_bad_cells(self):
        for cells in ([],self.cells()+[self.cells()[0]], [dict(self.cells()[0],area_m2='0')]):
            with self.assertRaises(ValueError):self.call(cells)
        with self.assertRaises(ValueError):self.call(independent_strips=False,coupling_margin_bound='-1')
    def test_evidence_does_not_admit(self):self.assertEqual(self.call(mode='evidence')['decision'],'UNCERTAIN')

class HookTests(unittest.TestCase):
    def request(self,gd='1/4'):
        return dict(kind='supplied_intervals',id='supplied-model',gp=['1','1'],gd=[gd,gd],
            gamma_layer=['1','1'],gamma_interface=['1/2','1/2'],source='test_model',
            model_context='explicit supplied common state')
    def test_optional_default(self):self.assertIsNone(incision_path_gate(None))
    def test_straight_guard(self):
        good=incision_path_gate({'requests':[self.request()]},mode='scenario')
        bad=incision_path_gate({'requests':[self.request('3/4')]},mode='scenario')
        self.assertFalse(good['blocks_straight_chain']);self.assertTrue(bad['blocks_straight_chain'])
        self.assertTrue(verify_integrity(good['boxes'][0]))
    def test_evidence_guard(self):
        self.assertTrue(incision_path_gate({'requests':[self.request()]})['blocks_straight_chain'])
    def test_mode_ids_and_empty(self):
        for cfg in ({'requests':[]},{'requests':[self.request(),self.request()]},
                    {'requests':[dict(self.request(),mode='scenario')]}, {'requests':[dict(self.request(),id=None)]}):
            with self.assertRaises(ValueError):incision_path_gate(cfg,mode='scenario')

if __name__=='__main__':unittest.main(verbosity=2)
