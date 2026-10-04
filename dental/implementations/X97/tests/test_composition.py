import json, copy, hashlib, unittest
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from occlusion_module import analyze, load_demo_input, geometry_sha256
from occlusion_module.api import ContractError, _prepare, bite_sha256
from occlusion_module._vendor import height_uncertainty, regional
R = Path(__file__).resolve().parents[1]

class CompositionTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.row = json.loads((R / 'inputs/cohort.json').read_text())[1]

    def setUp(self):
        (self.b, self.c) = load_demo_input(R, self.row)

    def packet(self):
        t = _prepare(self.b, self.c)
        x = t['vertices'][:, 0]
        u = (x - x.min()) / np.ptp(x)
        B = np.c_[t['taper'] * (1 - u), t['taper'] * u]
        J = np.array([[100.0, 0], [0, 100], [-100, 0], [0, -100]])
        f = np.array([10.0, 8, 10, 12])
        s = 0.01
        p = dict(schema='regional-actuation-probes-v1', units={'height': 'mm', 'force': 'N'}, observation_operator='calibrated_regional_axial_reactions', linear_support_closure=True, matched_load_rate_history=True, height_error_mm=0.0001, baseline_force_N=f.tolist(), plus_force_N=(f + s * J.T).tolist(), minus_force_N=(f - s * J.T).tolist(), step_mm=s, raw_force_channel_bound_N=0.01, geometry_sha256=geometry_sha256(self.c), basis_sha256=hashlib.sha256(B.astype('<f8').tobytes()).hexdigest(), case=self.b['case_id'], frame=self.b['frame_id'], pose_id=self.b['pose_id'], bite_sha256=bite_sha256(self.b), kind='SIMULATION')
        self.b['regional_actuation'] = dict(probes=p, basis=B, reaction_region_ids=['crownA', 'crownB', 'otherA', 'otherB'], reaction_tooth_fdi=[36, 36, 34, 31], final_coefficient_error_mm=[1e-05, 1e-05])
        return (p, J, f)

    def test_latest_local_constraint_parity(self):
        a = analyze(self.b, self.c, {'local_boundary': True})
        old = json.loads((R / 'inputs/X95_ROUND3.json').read_text())['rows'][1]['outputs']['regional']
        self.assertAlmostEqual(a['proposed_adjustment']['contact']['symdiff_mm2'], old['contact']['symdiff_mm2'], places=7)
        self.assertTrue(a['proposed_adjustment']['gates']['local_vertical_interval_order_pass'])

    def test_constraint_geometry_and_cap_faults(self):
        for kind in ['vertex', 'cap', 'margin', 'constraint']:
            c = copy.deepcopy(self.c)
            if kind == 'vertex':
                c['vertices_mm'][0, 0] += 0.001
            elif kind == 'cap':
                c['relief_cap_mm'] *= 0.9
            elif kind == 'margin':
                c['margin_z_mm'] += 0.01
            else:
                c['local_boundary_certificate']['b'][0] += 0.01
            with self.assertRaises(ContractError):
                analyze(self.b, c, {'local_boundary': True, 'propose_adjustment': False})

    def test_actual_local_constraint_injected_crossing(self):
        q = self.c['local_boundary_certificate']
        A = csr_matrix((q['data'], q['indices'], q['indptr']), shape=tuple(q['shape']))
        self.assertGreater(A.shape[0], 0)
        row = A.getrow(0)
        d = np.zeros(A.shape[1])
        j = int(row.indices[np.argmax(row.data)])
        d[j] = (q['b'][0] + 0.01) / row[0, j]
        self.assertGreater(float((A @ d - q['b']).max()), 1e-08)

    def test_main_regional_calibration_branch(self):
        (p, J, f) = self.packet()
        a = analyze(self.b, self.c, {'propose_adjustment': False, 'regional_force_design': True})['regional_force_design']
        self.assertEqual(a['status'], 'CONDITIONAL_TARGET_DESIGN')
        truth = f + J @ np.array(a['coefficients_mm'])
        iv = np.array(a['at_selected_height']['force_interval_N'])
        self.assertTrue(np.all(truth >= iv[:, 0]) and np.all(truth <= iv[:, 1]))
        self.assertGreater(truth[0] + 3, iv[0, 1])
        self.assertEqual(a['geometry_gates']['protected_identity_error_mm'], 0)

    def test_main_regional_identity_and_calibration_faults(self):
        for field in ['case', 'frame', 'pose_id', 'geometry_sha256', 'basis_sha256', 'bite_sha256']:
            (p, J, f) = self.packet()
            p[field] = 'other'
            with self.assertRaises(ContractError):
                analyze(self.b, self.c, {'propose_adjustment': False, 'regional_force_design': True})
        (p, J, f) = self.packet()
        p['kind'] = 'MEASURED'
        self.assertEqual(analyze(self.b, self.c, {'propose_adjustment': False, 'regional_force_design': True})['regional_force_design']['status'], 'UNKNOWN')

    def test_probe_denominator_fault(self):
        (p, J, f) = self.packet()
        p['height_error_mm'] = p['step_mm']
        with self.assertRaises(ValueError):
            analyze(self.b, self.c, {'propose_adjustment': False, 'regional_force_design': True})

    def test_regional_missing_observation_cannot_generate_force_design(self):
        a = analyze(self.b, self.c, {'propose_adjustment': False, 'regional_force_design': True})
        self.assertEqual(a['regional_force_design']['status'], 'UNKNOWN')

    def test_foreign_tooth_response_blocked_in_main(self):
        s = json.loads((R / 'inputs/F5367_R4_state.json').read_text())
        self.b['height_response'] = s
        self.b['force_context'] = dict(case_id=s['case'], response_geometry_sha256=s['geometry_sha256'], pose_id=self.b['pose_id'], calibration_pose_id=self.b['pose_id'], matched_load_rate_history=True, reference_crown_sha256=geometry_sha256(self.c))
        with self.assertRaises(ContractError):
            analyze(self.b, self.c, {'propose_adjustment': False, 'height_query': {'operation': 'uniform_tooth_height', 'fdi': 16, 'height_change_mm': -0.02}})

    def test_dependent_basis_rejected_even_with_matching_digest(self):
        (p, J, f) = self.packet()
        B = self.b['regional_actuation']['basis']
        B[:, 1] = B[:, 0]
        p['basis_sha256'] = hashlib.sha256(B.astype('<f8').tobytes()).hexdigest()
        with self.assertRaises(ContractError):
            analyze(self.b, self.c, {'propose_adjustment': False, 'regional_force_design': True})

    def test_deformable_scenario_must_be_explicit_and_bound(self):
        sc = dict(case_id=self.b['case_id'], frame_id=self.b['frame_id'], pose_id=self.b['pose_id'], geometry_sha256=geometry_sha256(self.c), evidence_kind='SIMULATION', geometry={'case': self.b['case_id']}, parameters={})
        self.b['deformable_scenario'] = sc
        with self.assertRaises(ContractError):
            analyze(self.b, self.c, {'propose_adjustment': False})
        sc['evidence_kind'] = 'MEASURED'
        with self.assertRaises(ContractError):
            analyze(self.b, self.c, {'propose_adjustment': False})
