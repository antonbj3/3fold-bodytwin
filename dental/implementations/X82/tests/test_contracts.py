import copy, json, sys, unittest
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R / 'code'))
from height_force import predict, inverse, make_edit, project_reaction, at_load
from run_r0_r1 import conventional

class Contracts(unittest.TestCase):

    def setUp(self):
        self.s = json.loads((R / 'inputs/F5367_R4_state.json').read_text())
        self.i = self.s['predicted_fdi'].index(16)

    def test_invalid_uncertainties_before_unknown(self):
        for key in ['spectral_error_N_per_mm', 'baseline_error_N']:
            for value in [-1, np.nan, np.inf, True]:
                s = copy.deepcopy(self.s)
                s[key] = value
                s['complete_active_calibration'] = False
                with self.subTest(key=key, value=str(value)), self.assertRaises(ValueError):
                    predict(s, make_edit(s, self.i, -0.02))

    def test_reaction_rejects_reversed_boxes(self):
        for (fe, de) in [(-0.1, 0), (0, -0.1), (np.nan, 0), (0, np.inf)]:
            with self.assertRaises(ValueError):
                project_reaction(1, 2, 0.005, 0.01, de, fe, True, True, True)

    def test_identity_unit_domain_and_finite(self):
        e = make_edit(self.s, self.i, -0.02)
        for (key, value) in [('unit', 'um'), ('case', 'other'), ('geometry_sha256', 'bad')]:
            bad = dict(e)
            bad[key] = value
            with self.assertRaises(ValueError):
                predict(self.s, bad)
        for h in [0.051, np.nan]:
            with self.assertRaises(ValueError):
                predict(self.s, make_edit(self.s, self.i, h))

    def test_unknown_not_default_stiffness(self):
        for (field, value) in [('complete_active_calibration', False), ('radius_model', 'unsupported')]:
            s = dict(self.s)
            s[field] = value
            a = predict(s, make_edit(s, self.i, -0.02))
            self.assertEqual(a['status'], 'UNKNOWN')
            self.assertIsNone(a['force_interval_N'])
        a = predict(self.s, {**make_edit(self.s, self.i, -0.02), 'load_N': 315})
        self.assertEqual(a['status'], 'UNKNOWN')

    def test_negative_height_can_release_contact(self):
        a = predict(self.s, make_edit(self.s, self.i, -0.05))
        self.assertTrue(np.all(np.array(a['force_interval_N'])[:, 0] >= 0))
        self.assertGreaterEqual(min(a['force_N']), 0)

    def test_inverse_guarantee_and_rejection(self):
        a = inverse(self.s, 16, 11.252173809585972)
        self.assertEqual(a['status'], 'CONDITIONAL_TARGET_BAND')
        i = self.i
        (lo, hi) = a['at_selected_height']['force_interval_N'][i]
        self.assertGreaterEqual(lo, a['target_band_N'][0])
        self.assertLessEqual(hi, a['target_band_N'][1])
        self.assertEqual(inverse(self.s, 16, 1000)['status'], 'INFEASIBLE_IN_CALIBRATED_BOX')
        with self.assertRaises(ValueError):
            inverse(self.s, 16, -1)

    def test_analytic_radius_with_actual_response_perturbations(self):
        rng = np.random.default_rng(8201)
        N = np.array(self.s['basis_N'])
        H = np.array(self.s['H_N_per_mm'])
        f0 = np.array(self.s['baseline_force_N'])
        m = H.shape[0]
        for k in range(24):
            E = rng.normal(size=(m, m))
            E = (E + E.T) / 2
            E *= self.s['spectral_error_N_per_mm'] * 0.95 / np.linalg.norm(E, 2)
            e = rng.normal(size=m)
            e *= self.s['baseline_error_N'] * 0.95 / np.linalg.norm(e)
            truth = dict(self.s, H_N_per_mm=(H + E).tolist(), baseline_force_N=(f0 + N @ e).tolist())
            h = [-0.04, 0.04, -0.02, 0.02][k % 4]
            edit = make_edit(self.s, self.i, h)
            ans = predict(self.s, edit)
            fc = conventional(truth, np.array(edit['height_change_mm']))
            bounds = np.array(ans['force_interval_N'])
            with self.subTest(k=k):
                self.assertLessEqual(np.linalg.norm(fc - ans['force_N']), ans['joint_l2_error_radius_N'] + 1e-06)
                self.assertTrue(np.all(fc >= bounds[:, 0] - 1e-06) and np.all(fc <= bounds[:, 1] + 1e-06))
                inv = inverse(self.s, 16, 11.252173809585972)
                br = inv['possible_exact_target_height_bracket_mm']
                fl = conventional(truth, np.array(make_edit(truth, self.i, br[0])['height_change_mm']))[self.i]
                fr = conventional(truth, np.array(make_edit(truth, self.i, br[1])['height_change_mm']))[self.i]
                self.assertLessEqual(fl, inv['target_N'] + 0.0001)
                self.assertGreaterEqual(fr, inv['target_N'] - 0.0001)

    def test_load_column_is_new_information(self):
        s = json.loads((R / 'exports/F5367_load_state.json').read_text())
        a = predict(s, {**make_edit(s, self.i, -0.02), 'load_N': 315})
        self.assertEqual(a['status'], 'CONDITIONAL_PASS')
        self.assertAlmostEqual(sum(a['force_N']), 315, places=6)
        for (key, value) in [('load_response_error_l2_N_per_N', -0.01), ('load_change_limit_N', np.nan)]:
            bad = dict(s)
            bad[key] = value
            with self.assertRaises(ValueError):
                predict(bad, {**make_edit(s, self.i, -0.02), 'load_N': 315})
        bad = copy.deepcopy(s)
        bad['load_response_N_per_N'][0] *= 1.2
        with self.assertRaises(ValueError):
            at_load(bad, 315)
        self.assertEqual(predict(s, {**make_edit(s, self.i, -0.02), 'load_N': 320})['status'], 'UNKNOWN')

class CorridorContracts(unittest.TestCase):

    def setUp(self):
        self.s = json.loads((R / 'exports/F5367_corridor_state.json').read_text())

    def test_corridor_refuses_extrapolation_and_missing_closure(self):
        from crown_corridor import predict
        self.assertEqual(predict(self.s, 0.04)['status'], 'UNKNOWN')
        self.assertEqual(predict(dict(self.s, linear_full_contact_closure=False), 0)['status'], 'UNKNOWN')
        self.assertEqual(predict(self.s, 0, w_N=[315, 0, 0])['status'], 'UNKNOWN')

    def test_height_uncertainty_is_a_quotient_enclosure(self):
        from crown_corridor import predict, error_contract
        ss = dict(self.s, probe_step_error_mm=0.0001, query_height_error_mm=0.0001)
        a = predict(ss, -0.02)
        self.assertIsNotNone(a['force_interval_N'])
        b = np.array(ss['baseline_force_N'])
        j = np.array(ss['height_column_N_per_mm'])
        (eta, beta, _, radius) = error_contract(ss, b, j, 0.025, 0.01, -0.02)
        u = j / np.linalg.norm(j)
        jtrue = (0.025 * j + eta * u) / (0.025 - 0.0001)
        ftrue = b - eta * u + jtrue * (-0.02 - 0.0001)
        self.assertLessEqual(np.linalg.norm(ftrue - np.array(a['force_N'])), radius + 1e-12)
        with self.assertRaises(ValueError):
            predict(dict(ss, probe_step_error_mm=-0.001), 0)
        with self.assertRaises(ValueError):
            predict(dict(ss, probe_step_error_mm=0.025), 0)

    def test_inverse_interval_contains_reference_root(self):
        from crown_corridor import inverse
        a = inverse(self.s, 11.252173809585972)
        self.assertEqual(a['status'], 'CONDITIONAL_TARGET_BAND')
        old = json.loads((R / 'inputs/F5367_R4_state.json').read_text())
        N = np.array(old['basis_N'])
        H = np.array(old['H_N_per_mm'])
        i = old['predicted_fdi'].index(16)
        j = N @ H @ N[i]
        true_h = (a['target_N'] - old['baseline_force_N'][i]) / j[i]
        self.assertLessEqual(a['possible_exact_target_height_bracket_mm'][0], true_h)
        self.assertGreaterEqual(a['possible_exact_target_height_bracket_mm'][1], true_h)
if __name__ == '__main__':
    unittest.main()
