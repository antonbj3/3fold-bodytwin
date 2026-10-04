"""Behavioral regression tests for the actual interval consumer boundary."""
import copy
import json
from pathlib import Path
import sys
import unittest
from fractions import Fraction as Q
from hum_grasp_certificate import certify, verify, source_request
from hum_robotport_adapter import validate_catalog, validate_interval, validate_base, validate_contact_port, compose_catalog

ROOT = Path(__file__).resolve().parents[1]
# Works both in this lane and in a locally applied BodyTwin patch.
CATALOG = ROOT/'ports/PORT_HUM_20261003.json'
if not CATALOG.is_file():
    CATALOG = Path(__file__).resolve().parent/'PORT_HUM_20261003.json'


class RobotPortTests(unittest.TestCase):
    def setUp(self):
        self.catalog = json.loads(CATALOG.read_text())

    def reject(self, mutate):
        mutate(self.catalog)
        with self.assertRaises(ValueError):
            validate_catalog(self.catalog)

    def test_native_field_layout_and_interval_payloads(self):
        self.assertTrue(validate_catalog(self.catalog)['ok'])

    def test_point_value_rejected_where_interval_required(self):
        self.reject(lambda c: c['ports'][0].update(computed_value=16))

    def test_point_value_in_auxiliary_quantity_rejected(self):
        self.reject(lambda c: c['ports'][0]['quantities'].update(vx=0))

    def test_nonfinite_float_and_boolean_endpoints_rejected(self):
        for x in [float('nan'), float('inf'), True, 12.0]:
            with self.subTest(x=x):
                p = copy.deepcopy(self.catalog['ports'][0]['computed_value'])
                p['lower'] = x
                with self.assertRaises(ValueError): validate_interval(p)

    def test_reversed_interval_rejected(self):
        self.reject(lambda c: c['ports'][0]['computed_value'].update(lower='21'))

    def test_force_unit_change_rejected(self):
        self.reject(lambda c: c['ports'][0]['computed_value'].update(unit='Pa'))

    def test_gap_unit_change_rejected(self):
        self.reject(lambda c: c['ports'][1]['inputs']['gap_finger0_m'].update(unit='mm'))

    def test_certificate_normal_force_corruption_rejected(self):
        self.reject(lambda c: c['ports'][0]['certificate']['payload']['bounds']['normal_force_per_finger'][0].update(upper='19'))

    def test_auxiliary_velocity_corruption_rejected(self):
        self.reject(lambda c: c['ports'][0]['quantities']['vy'].update(upper='1'))

    def test_input_mu_change_without_certificate_rejected(self):
        self.reject(lambda c: c['ports'][1]['inputs']['mu'].update(lower='3/5'))

    def test_missing_commit_or_source_rejected(self):
        self.reject(lambda c: c['ports'][0]['computed_origin'].update(commit=None))
        self.setUp()
        self.reject(lambda c: c['ports'][0]['computed_value'].update(source_refs=['unknown']))

    def test_physical_upgrade_rejected(self):
        self.reject(lambda c: c['ports'][1]['certificate'].update(physical_status='CERTIFIED'))

    def test_published_dispersion_upgrade_rejected(self):
        self.reject(lambda c: c['ports'][1]['inputs']['surface_Ra_m'].update(certificate_status='EXACT_CONDITIONAL_MODEL_BOUND'))

    def test_internal_result_cannot_be_labelled_independent_measurement(self):
        self.reject(lambda c: c['ports'][1].update(computed_provenance_class='INDEPENDENT_MEASUREMENT'))

    def test_unrelated_existing_port_cannot_receive_grasp_force(self):
        self.reject(lambda c: c['ports'][0].update(node_id='bodytwin:SURGX-P-TIP-FORCE-POINT-SENSING'))

    def test_worst_box_failure_not_hidden_by_nominal_hold(self):
        request = source_request()
        request['mu'] = [['3/10', '1/2']]*2
        cert = certify(request)
        self.assertEqual(cert['decision'], 'REGIME_DEPENDENT')
        self.assertTrue(cert['any_slip_witness'])

    def test_load_exceeding_entire_capacity_slips(self):
        request = source_request()
        request['free_velocity_m_s'] = ['0', '-1', '0']
        self.assertEqual(certify(request)['decision'], 'SLIPS_ENTIRE_BOX')

    def test_oblique_rotational_frame_has_no_license(self):
        request = source_request()
        request['J'][0] = [1, 1, 0]
        with self.assertRaises(ValueError): certify(request)

    def test_gap_box_refinement_preserves_nested_force_enclosures(self):
        outer = certify(source_request())
        request = source_request()
        request['gaps_m'] = [['-9/5000', '-7/5000']]*2
        inner = certify(request)
        a = outer['bounds']['tangential_capacity']; b = inner['bounds']['tangential_capacity']
        self.assertLessEqual(Q(a['lower']), Q(b['lower']))
        self.assertLessEqual(Q(b['upper']), Q(a['upper']))
        self.assertTrue(verify(inner))

    def test_native_bodytwin_validation_is_read_only(self):
        native = Path(self.catalog['base_catalog']['file']).parent
        self.assertEqual(validate_base(self.catalog, native)['original_records_unchanged'], 202)

    def test_composed_catalog_preserves_each_original_record(self):
        native = Path(self.catalog['base_catalog']['file']).parent
        original = json.loads((native/'PORT.json').read_text())
        view = compose_catalog(self.catalog, native)
        self.assertEqual(view['ports'][:202], original['ports'])
        self.assertEqual(view['summary'], original['summary'])
        self.assertEqual(len(view['ports']), 204)
        self.assertEqual(len(view['connected_ports']), 2)

    def test_actual_contact_port_v1_roundtrip_and_witness_replay(self):
        upstream = ROOT/'code/upstream'
        if upstream.is_dir(): sys.path.insert(0, str(upstream))
        try:
            import field_engine.contact_port_v1
        except ImportError:
            self.skipTest('Contact port v1 field_engine dependency required for this optional transport check')
        for p in self.catalog['ports']:
            text = (CATALOG.parent/p['contact_port_v1_file']).read_text()
            r = validate_contact_port(text, p['certificate']['payload'])
            self.assertEqual(r['branch_set_status'], 'UNKNOWN')
            self.assertEqual(r['physical_status'], 'UNKNOWN')


if __name__ == '__main__':
    unittest.main()
