"""One-time range registration; no arbitrary dental upper bounds."""
from common import *

def define():
    n = source('r3_net')
    contracts = {}
    structured = set('scan_surface region_id anatomy_surface pose jaw_pose sdf material_law muscle_force force_direction load_spectrum stress contact_traction tmj_reaction tool_path implant_geometry procedure_state material_history design_choice measured_response calibration_state artifact_hash replay_hash license_scope measurement_plan'.split())
    signed = set('signed_error contact_gap offset tooth_displacement apex_error nerve_distance mbl design_margin shape_gradient access_clearance asbuilt_error prediction_residual build_angle'.split())
    unknown = set('cbct_gray hu_ref isq cell_viability'.split())
    exact = {'bic': ('0', '1'), 'porosity': ('0', '1'), 'failure_risk': ('0', '1'), 'relative_tooth_force': ('0', '100')}
    for (k, v) in n['variables'].items():
        status = 'PARTIAL_DOMAIN'
        lo = '0'
        hi = None
        basis = 'Nonnegative magnitude by the variable definition; no supported finite upper bound.'
        if k in structured:
            status = 'NOT_SCALAR'
            lo = None
            basis = 'Categorical, field, tensor, vector, spectrum or mixed-unit geometry. Needs component and frame/observation operator, not one scalar interval.'
        elif k in signed:
            status = 'UNKNOWN'
            lo = None
            basis = 'Signed/context-dependent quantity; no universal finite bounds available in admitted references.'
        elif k in unknown:
            status = 'UNKNOWN'
            lo = None
            basis = 'Scanner/assay/instrument convention and regime required; no finite sourced interval. Viability assays may be normalized ratios exceeding 100, not count fractions.'
        elif k in exact:
            status = 'MATHEMATICAL_DOMAIN'
            (lo, hi) = exact[k]
            basis = 'Subset measure / total measure is in [0,1]; percentage share is 100 times that. Does not establish empirical plausibility.'
        elif k in ['pulp_temperature', 'drill_temperature']:
            status = 'PARTIAL_DOMAIN'
            lo = '-273.15'
            basis = 'Absolute temperature cannot be negative (Celsius to kelvin definition); no dental thermal envelope supplied.'
        contracts[k] = {'status': status, 'unit': v['unit'], 'lower': lo, 'upper': hi, 'inclusive': [True, True], 'source': {'kind': 'closed_form' if lo is not None else 'external_review', 'locator': entry('r3_net')['original_path'] + '#/variables/' + k, 'snapshot_sha256': entry('r3_net')['sha256'], 'compared_quantity': v['desc'], 'refutes_us': False}, 'basis': basis, 'resolution_level': v['resolution_level'], 'scope': 'Variable definition only; numeric ingestion abstains without finite scope unless mathematical domain is the claimed output.', 'finite_empirical_range_known': False, 'replacement_measurement': 'Matched range/reference for ' + v['desc'] + ' with population, protocol, time, region, units and uncertainty.', 'cases': []}
    obs = read('SOURCE_OBSERVATIONS.json')['records']
    (s, c, ch) = obs
    contracts['precementation_internal_gap'] = {'status': 'SCOPED_SOURCE_ENVELOPE', 'unit': 'um', 'lower': None, 'upper': None, 'inclusive': [True, True], 'source': c['reference'], 'basis': 'No universal film bound; cases hold source group means at the stated pre-cementation region.', 'resolution_level': 'PER_SURFACE_REGION', 'scope': 'Pre-cementation internal gap; not assembled cement film or marginal opening.', 'finite_empirical_range_known': False, 'replacement_measurement': 'Matched cemented film and dry/replica repeated measurements on same copings.', 'cases': []}

    def case(var, id, lo, hi, quantity, protocol, level, role, ref, basis, context):
        r = contracts[var]
        r['status'] = 'SCOPED_SOURCE_ENVELOPE'
        r['source'] = ref
        r['basis'] = 'Finite only inside a listed source context; target-system plausibility UNKNOWN.'
        r['cases'].append({'id': id, 'lower': str(lo), 'upper': str(hi), 'inclusive': [True, True], 'unit': r['unit'], 'quantity': quantity, 'protocol_id': protocol, 'resolution_level': level, 'value_role': role, 'source': ref, 'basis': basis, 'context': context, 'statistical_coverage': 'Observed support / prescribed source range, not confidence or prediction interval', 'out_of_range_meaning': 'Reject or quarantine for this source scope; does not prove physical impossibility elsewhere.'})
    sp = s['protocol']['id']
    cp = c['protocol']['id']
    hp = ch['protocol']['id']
    for sequence in ['first', 'tenth']:
        (lo, hi) = s['reported'][sequence + '_range_N']
        case('preload', 'SAGHEB_' + sequence.upper(), lo, hi, 'axial_preload', sp, 'POPULATION', 'reported_group_mean', s['reference'], 'A source group mean lies in its measured specimen range.', {'sequence': sequence, 'fixture': 'separated_0.10_mm', 'surface': 'carbon_coated_dry'})
        case('preload', 'SAGHEB_' + sequence.upper() + '_SPECIMEN', lo, hi, 'axial_preload', sp, 'PER_TOOTH', 'reported_specimen_value', s['reference'], 'Observed endpoint specimen support only; no future specimen guarantee.', {'sequence': sequence, 'fixture': 'separated_0.10_mm', 'surface': 'carbon_coated_dry'})
    case('tightening_torque', 'SAGHEB_TORQUE', 25, 25, 'applied_screw_torque', sp, 'PER_TOOTH', 'prescribed_setpoint', s['reference'], 'Exact nominal protocol setpoint, not measurement uncertainty or general torque range.', {'fixture': 'separated_0.10_mm'})
    case('precementation_internal_gap', 'CUNALI_AMANN_MO', c['reported']['ct_mean_um'], c['reported']['replica_mean_um'], 'mid_occlusal_internal_gap', cp, 'POPULATION', 'reported_group_mean', c['reference'], 'Envelope of the two audited group means only. No raw individual min/max known.', {'material': 'Amann_Ceramill_Zi', 'region': 'mid_occlusal', 'state': 'precementation', 'methods': ['dry_micro_CT', 'PVS_replica']})
    for group in ['TX', 'EV']:
        rows = [x for x in ch['reported']['fatigue_rows'] if x['group'] == group]
        case('fatigue_load', 'CHOI_' + group, min((x['max_N'] for x in rows)), max((x['max_N'] for x in rows)), 'fatigue_peak_load', hp, 'PER_TOOTH', 'prescribed_peak', ch['reference'], 'Range of explicitly tested peak levels only; a continuous range does not license interpolated fatigue outcomes.', {'system': group, 'R': '0.1', 'frequency_Hz': 15, 'runout_cycles': 5000000})
        raw = [x[group] for x in ch['reported']['static_individual_N']]
        case('fracture_mean', 'CHOI_STATIC_' + group, min(raw), max(raw), 'static_maximum_breaking_load', hp, 'POPULATION', 'reported_group_mean', ch['reference'], 'Reported mean must lie within its three-specimen measured range.', {'system': group, 'rate_mm_min': 1, 'failure_definition': 'fracture_or_deformation'})
    write('RANGE_CONTRACTS.json', {'schema': 'scoped-plausible-ranges-r4', 'frozen_utc': utc(), 'claim_type': 'capability', 'variables': contracts, 'policy': 'Range completeness and scientific validation are separate. UNKNOWN/PARTIAL_DOMAIN/NOT_SCALAR reject scalar admission; source envelopes apply only to exactly matched source contexts. No empirical universal bounds invented.'})
    (H / 'RANGE_CONTRACTS.json.sha256').write_text(sha(H / 'RANGE_CONTRACTS.json') + '\n')
    print({'variables': len(contracts), 'finite_source_variables': sum((bool(x['cases']) for x in contracts.values())), 'finite_source_cases': sum((len(x['cases']) for x in contracts.values()))})
if __name__ == '__main__':
    define()
