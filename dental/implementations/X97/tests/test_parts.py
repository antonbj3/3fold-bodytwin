import unittest, json, hashlib, copy
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
from scipy.optimize import minimize
from occlusion_module.api import _force_bounds, query_height, regional_design, ContractError, analyze, load_demo_input, geometry_sha256
from occlusion_module._vendor import contact, quality, reachability, height_force, deformable, force_port
R = Path(__file__).resolve().parents[1]

class PartTests(unittest.TestCase):

    def setUp(self):
        self.row = json.loads((R / 'inputs/cohort.json').read_text())[3]
        (self.b, self.c) = load_demo_input(R, self.row)

    def graph(self, total=None):
        b = {'teeth_fdi': [16, 26, 36, 46, 17], 'pair_model': 'AXIAL_POSITIVE_SUPPORT_V1', 'tooth_pairs': [dict(upper_fdi=16, lower_fdi=46, gap_interval_mm=['0', '0']), dict(upper_fdi=26, lower_fdi=36, gap_interval_mm=['0.1', '0.1'])]}
        if total is not None:
            b.update(total_force_interval_N=total, total_force_evidence={'kind': 'SIMULATED', 'locator': 'tests/test_parts.py:graph'})
        return b

    def test_tandlast_must_can_never(self):
        a = {r['fdi']: r for r in _force_bounds(self.graph([100, 100]), {'tooth_fdi': 16})['per_tooth']}
        self.assertEqual([a[t]['classification'] for t in [16, 26, 17]], ['MUST', 'CAN', 'NEVER'])
        self.assertEqual(a[16]['force_interval_N'], [0.0, 100.0])
        self.assertFalse(a[16]['infimum_attained'])

    def test_no_inherited_population_force(self):
        a = _force_bounds(self.graph(), {'tooth_fdi': 16})
        self.assertTrue(a['upper_bound_missing'])
        self.assertEqual(a['per_tooth'][0]['force_interval_N'], [0.0, None])

    def test_total_force_rejects_uncalibrated_measurement(self):
        b = self.graph([100, 100])
        b['total_force_evidence']['kind'] = 'MEASURED'
        with self.assertRaises(ContractError):
            _force_bounds(b, {'tooth_fdi': 16})

    def test_graph_zero_and_duplicate_faults(self):
        for fault in ['zero', 'duplicate', 'orientation', 'reverse']:
            b = self.graph([100, 100])
            if fault == 'zero':
                b['total_force_interval_N'] = [0, 100]
            elif fault == 'duplicate':
                b['tooth_pairs'] *= 2
            elif fault == 'orientation':
                b['tooth_pairs'][0]['upper_fdi'] = 46
            else:
                b['tooth_pairs'][0]['gap_interval_mm'] = ['1', '0']
            with self.assertRaises(ContractError):
                _force_bounds(b, {'tooth_fdi': 16})

    def test_gap_tie_preserves_strict_positive(self):
        b = self.graph([100, 100])
        b['tooth_pairs'][1]['gap_interval_mm'] = ['0', '0']
        a = _force_bounds(b, {'tooth_fdi': 16})
        self.assertTrue(all((r['classification'] == 'MUST' for r in a['per_tooth'] if r['fdi'] != 17)))

    def test_contact_polygon_independent(self):
        xy = np.array([[0.0, 0], [1, 0], [0, 1]])
        f = np.array([[0, 1, 2]])
        pg = np.array([-0.1, 0.2, 0.2])
        rg = np.array([0.05] * 3)
        a = contact.compare(xy, f, pg, rg)
        b = quality.contact_map(xy, f, pg, rg)
        self.assertAlmostEqual(a['symdiff_mm2'], b['contact_symdiff_mm2'], places=12)
        shifted = pg + 0.2
        self.assertGreater(contact.compare(xy, f, shifted, rg)['symdiff_mm2'], a['symdiff_mm2'])

    def test_contact_summary_insufficient(self):
        w = reachability.sufficiency(contact)
        self.assertTrue(w['bit_identical'])
        self.assertEqual(w['summary_identity_error'], 0)
        self.assertGreater(w['downstream_difference_mm2'], 0)

    def test_geometry_identity_faults(self):
        for (k, v) in [('case_id', 'other'), ('frame_id', 'wrong'), ('pose_id', 'loaded_other'), ('unit', 'um')]:
            c = dict(self.c, **{k: v})
            with self.assertRaises(ContractError):
                analyze(self.b, c, {'propose_adjustment': False})

    def test_input_arrays_not_mutated(self):
        old = self.c['vertices_mm'].copy()
        a = analyze(self.b, self.c, {'propose_adjustment': False})
        self.assertTrue(np.array_equal(old, self.c['vertices_mm']))
        self.assertEqual(a['force_intervals']['status'], 'UNKNOWN')
        self.assertEqual(a['physical_validation'], 'NOT_RUN')

    def test_cap_wall_fault(self):
        c = dict(self.c, relief_cap_mm=1.0)
        with self.assertRaises(ContractError):
            analyze(self.b, c, {'propose_adjustment': False})

    def test_x82_interval_independent_qp_and_force_fault(self):
        s = json.loads((R / 'inputs/F5367_R4_state.json').read_text())
        i = s['predicted_fdi'].index(16)
        b = dict(case_id=s['case'], response_geometry_sha256=s['geometry_sha256'], pose_id='cal0', calibration_pose_id='cal0', matched_load_rate_history=True)
        a = query_height(b, s, {'operation': 'uniform_tooth_height', 'fdi': 16, 'height_change_mm': -0.02})
        (A, N, H, f0, w, hb, rho, eta) = height_force.validate(s)
        C = N @ np.linalg.inv(H) @ N.T + (np.eye(len(f0)) - N @ N.T) / 750
        h = np.zeros(len(f0))
        h[i] = -0.02
        g = -C @ f0
        o = minimize(lambda f: 0.5 * f @ C @ f + (g - h) @ f, f0, jac=lambda f: C @ f + g - h, bounds=[(0, None)] * len(f0), constraints={'type': 'eq', 'fun': lambda f: A.T @ f - w, 'jac': lambda f: A.T}, method='SLSQP', options={'ftol': 1e-12, 'maxiter': 1000})
        self.assertTrue(o.success)
        iv = np.array(a['force_interval_N'])
        self.assertTrue(np.all(o.x >= iv[:, 0]) and np.all(o.x <= iv[:, 1]))
        self.assertGreater(o.x[i] + 3, iv[i, 1])
        self.assertAlmostEqual(a['force_N'][i], 11.66883172, places=5)

    def test_x82_rejects_region_and_identity(self):
        s = json.loads((R / 'inputs/F5367_R4_state.json').read_text())
        b = dict(case_id=s['case'], response_geometry_sha256=s['geometry_sha256'], pose_id='cal0', calibration_pose_id='cal0', matched_load_rate_history=True)
        self.assertEqual(query_height(b, s, {'operation': 'regional_height'})['status'], 'UNKNOWN')
        for k in ['case_id', 'response_geometry_sha256', 'calibration_pose_id']:
            with self.assertRaises(ContractError):
                query_height(dict(b, **{k: 'other'}), s, {'operation': 'uniform_tooth_height'})

    def probes(self):
        B = np.eye(2)
        digest = hashlib.sha256(B.astype('<f8').tobytes()).hexdigest()
        J = np.array([[100.0, 0], [0, 100], [-100, 0], [0, -100]])
        f = np.array([10.0, 8, 10, 12])
        step = 0.01
        p = dict(schema='regional-actuation-probes-v1', units={'height': 'mm', 'force': 'N'}, observation_operator='calibrated_regional_axial_reactions', linear_support_closure=True, matched_load_rate_history=True, height_error_mm=0.0, baseline_force_N=f.tolist(), plus_force_N=(f + step * J.T).tolist(), minus_force_N=(f - step * J.T).tolist(), step_mm=step, raw_force_channel_bound_N=0.01, geometry_sha256='declared_fixture_geometry', basis_sha256=digest, case='simulation_fixture', frame='roof_frame', kind='SIMULATION')
        b = {k: p[k] for k in ['case', 'frame', 'geometry_sha256', 'basis_sha256']}
        return (b, p, B, J, f)

    def test_regional_target_and_fault(self):
        (b, p, B, J, f) = self.probes()
        a = regional_design(b, p, B, np.ones(2), np.full(2, 0.1))
        self.assertEqual(a['status'], 'CONDITIONAL_TARGET_DESIGN')
        truth = f + J @ np.asarray(a['coefficients_mm'])
        iv = np.array(a['at_selected_height']['force_interval_N'])
        self.assertTrue(np.all(truth >= iv[:, 0]) and np.all(truth <= iv[:, 1]))
        self.assertGreater(truth[0] + 3, iv[0, 1])

    def test_regional_missing_closure_and_basis_fault(self):
        (b, p, B, J, f) = self.probes()
        p['linear_support_closure'] = False
        self.assertEqual(regional_design(b, p, B, np.ones(2), np.full(2, 0.1))['status'], 'UNKNOWN')
        (b, p, B, J, f) = self.probes()
        with self.assertRaises(ContractError):
            regional_design(b, p, B[:, ::-1], np.ones(2), np.full(2, 0.1))

    def test_x54_deformable_scenario_and_ncp_fault(self):

        def patch(u, l, x):
            return dict(upper_fdi=u, lower_fdi=l, region=0, gap_mm=0.0, area_mm2=1.0, basis=[[0, 0, 1], [1, 0, 0], [0, 1, 0]], upper_point_mm=[x, 0, 0], lower_point_mm=[x, 0, 0])
        g = dict(case='simulation_fixture', patches=[patch(16, 46, -1.0), patch(26, 36, 1.0)], centers={'upper': {'16': [-1.0, 0, 0], '26': [1.0, 0, 0]}, 'lower': {'46': [-1.0, 0, 0], '36': [1.0, 0, 0]}})
        (a, _) = deformable.solve(g, mu=0.0, controls=True)
        self.assertTrue(a['numerics']['all_pass'])
        self.assertTrue(a['controls']['injected_force_error']['rejected'])
        self.assertTrue(a['controls']['frictionless_NNLS']['pass_gate'])

    def test_force_port_calibration_fault(self):
        o = dict(context=dict(time_relation='SIMULTANEOUS', unit='force_share_same_total', same_pose_id='x', common_denominator_id='T', incidence_sha256='x', evidence_kind='MEASURED'), upper_region_interval_share=['1/2', '1/2'], signed_net_crossing_interval_share=['1/10', '1/10'])
        self.assertEqual(force_port.project(o)['status'], 'UNKNOWN')
        o['context'].update(calibration_locator='simulation_test_only', gain_uncertainty_status='BOUNDED')
        self.assertEqual(force_port.project(o)['lower_region_interval_share'], ['3/5', '3/5'])
        o['signed_net_crossing_interval_share'] = ['1', '1']
        self.assertEqual(force_port.project(o)['status'], 'REJECTED')
if __name__ == '__main__':
    unittest.main()
