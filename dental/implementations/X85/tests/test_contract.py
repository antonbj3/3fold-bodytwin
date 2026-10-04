import sys, unittest, copy
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from common import *
from regional import fit_response, predict, inverse
from plugin import generate
from geometry import ball_access
from round2 import load_item, simulated_calibration
from round5 import full_crown_force_query
from height_uncertainty import fit_bounded_heights
from lab_compare import compare as lab_compare

class Contracts(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        (t, z, inn) = load_item(read(ROOT / 'PREREG_R2.json')['cohort'][0])
        B = np.array(np.load(read(ROOT / 'raw/R2_ROWS.json')[0]['artifact_path'])['basis'])
        (model, probes, *_) = simulated_calibration(t, z, inn, B)
        cls.task = dict(t, regional_calibration=model, actuation_basis=B, unadjusted_z=z, inner_z=inn, geometry_sha256=model['geometry_sha256'], basis_sha256=model['basis_sha256'])
        cls.probes = clean(probes)
        cls.model = model

    def test_exact_regional_summary_counterexample(self):
        r = read(ROOT / 'raw/R0_SUFFICIENCY.json')
        self.assertEqual(r['summary_identity_error'], 0)
        self.assertTrue(r['bitidentical_summary'])
        self.assertGreaterEqual(r['downstream_max_region_force_difference_N'], 0.1)

    def test_frame_fault_rejects(self):
        t = copy.deepcopy(self.task)
        t['frame'] = 'other'
        with self.assertRaises(ValueError):
            generate(t)

    def test_actual_geometry_fault_rejects(self):
        t = copy.deepcopy(self.task)
        t['unadjusted_z'] = np.array(t['unadjusted_z']) + 0.001
        with self.assertRaises(ValueError):
            generate(t)

    def test_basis_bytes_fault_rejects(self):
        t = copy.deepcopy(self.task)
        t['actuation_basis'] = np.array(t['actuation_basis'])
        t['actuation_basis'][0, 0] += 0.001
        with self.assertRaises(ValueError):
            generate(t)

    def test_unknown_observation_is_not_force(self):
        p = copy.deepcopy(self.probes)
        p['observation_operator'] = 'TScan_relative_percent'
        self.assertEqual(fit_response(p)['status'], 'UNKNOWN')

    def test_missing_linearity_and_history(self):
        for key in ['linear_support_closure', 'matched_load_rate_history']:
            p = copy.deepcopy(self.probes)
            p[key] = False
            self.assertEqual(fit_response(p)['status'], 'UNKNOWN')

    def test_unenclosed_height_error_abstains(self):
        p = copy.deepcopy(self.probes)
        p['height_error_mm'] = 0.001
        self.assertEqual(fit_response(p)['status'], 'UNKNOWN')

    def test_bounded_height_denominator_and_worlds(self):
        p = copy.deepcopy(self.probes)
        p['height_error_mm'] = p['step_mm']
        with self.assertRaises(ValueError):
            fit_bounded_heights(p, [0.0005, 0.0005])
        r = read(ROOT / 'raw/R6_SUMMARY.json')
        self.assertEqual(r['status'], 'FAIL')
        self.assertEqual(r['enclosed'], 0)
        self.assertTrue(r['zero_denominator_rejected'])
        r = read(ROOT / 'raw/R7_SUMMARY.json')
        self.assertEqual(r['enclosed'], r['worlds'])
        self.assertEqual(r['faults_rejected'], r['worlds'])

    def test_old_exact_height_bound_can_fail(self):
        p = copy.deepcopy(self.probes)
        s = p['step_mm']
        d = 0.001
        f0 = np.array(p['baseline_force_N'])
        J = np.array(self.model['J_N_per_mm'])
        p.update(raw_force_channel_bound_N=0.0, height_error_mm=d, plus_force_N=np.array([f0 + J[:, j] * (s - d) for j in range(2)]).tolist(), minus_force_N=np.array([f0 - J[:, j] * (s - d) for j in range(2)]).tolist())
        new = fit_bounded_heights(p)
        oldp = copy.deepcopy(p)
        oldp['height_error_mm'] = 0.0
        old = fit_response(oldp)
        a = np.array([-0.02, 0.0])
        truth = f0 + J @ a
        bad = np.array(predict(old, a)['force_interval_N'])
        good = np.array(predict(new, a)['force_interval_N'])
        self.assertTrue(np.any(truth < bad[:, 0] - 1e-07) | np.any(truth > bad[:, 1] + 1e-07))
        self.assertTrue(np.all(truth >= good[:, 0] - 1e-07) & np.all(truth <= good[:, 1] + 1e-07))

    def test_negative_force_error_rejects(self):
        p = copy.deepcopy(self.probes)
        p['raw_force_channel_bound_N'] = -0.05
        with self.assertRaises(ValueError):
            fit_response(p)

    def test_same_information_qp_and_faults(self):
        rows = read(ROOT / 'raw/R2_CONTROLS.json')
        for c in rows:
            for h in c.get('heldouts', []):
                self.assertTrue(h['enclosed'])
                self.assertLess(h['maximum_nominal_error_N'], 1e-05)
                self.assertTrue(h['injected_plus3_N_rejected'])
            if 'wall_fault_rejected' in c:
                self.assertTrue(c['wall_fault_rejected'])

    def test_tool_can_reject(self):
        c = read(ROOT / 'raw/R2_CONTROLS.json')[-1]
        self.assertTrue(c['false_milling_pass_rejected'])

    def test_original_plugin_pinned(self):
        for (p, d) in read(ROOT / 'INPUT_LOCK_R2.json')['files'].items():
            self.assertEqual(sha(p), d['sha256'])

    def test_whole_crown_missing_calibration(self):
        q = full_crown_force_query(dict(mesh_sha256='actual', frame='A'))
        self.assertEqual(q['status'], 'UNKNOWN')
        self.assertIsNone(q['force_interval_N'])
        with self.assertRaises(ValueError):
            full_crown_force_query(dict(mesh_sha256='actual', frame='A'), dict(geometry_sha256='other', frame='A'))

    def test_external_count_insufficiency_and_fault(self):
        e = read(ROOT / 'raw/EXTERNAL_MEASUREMENT.json')
        self.assertEqual(e['contact_count_identity_error'], 0)
        self.assertEqual(e['relative_force_difference_percentage_points'], 7.0)
        self.assertTrue(e['fault_plus1_signal_rejected'])

    def test_full_crown_geometric_controls(self):
        for c in read(ROOT / 'raw/R5_CONTROLS.json'):
            self.assertTrue(c['protected_vertices_identical'])
            self.assertTrue(c['force_binding_fault_rejected'])
            self.assertLess(c['contact_polygon_parity_error_mm2'], 1e-07)
            self.assertTrue(c['fault_plus3mm2_rejected'])

    def test_lab_consumer_rejects_fault_without_fit(self):
        frozen = dict(case='S', frame='R', geometry_sha256='G', force_interval_N=[[7.0, 9.0], [7.0, 9.0]])
        measured = dict(case='S', frame='R', geometry_sha256='G', unit='N', matched_load_rate_history=True, reaction_N=[8.0, 8.0], reaction_error_N=[0.1, 0.1])
        self.assertEqual(lab_compare(frozen, measured)['status'], 'PASS_FROZEN_OBSERVATION_GATE')
        measured['reaction_N'][0] += 3
        self.assertEqual(lab_compare(frozen, measured)['status'], 'FAIL_FROZEN_OBSERVATION_GATE')
        measured['geometry_sha256'] = 'other'
        with self.assertRaises(ValueError):
            lab_compare(frozen, measured)
if __name__ == '__main__':
    unittest.main()
