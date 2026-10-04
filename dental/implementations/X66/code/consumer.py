"""Source-observation ports. They do not infer new-device clinical performance."""
from common import *

def query(lid, q):
    data = read('raw/observations.json')
    r1 = read('R1_RESULTS.json')

    def unknown(reason):
        return {'status': 'UNKNOWN', 'reason': reason, 'prediction': None}
    if lid == 'L01':
        return unknown('Primary ramp-time table absent; static0.5s digitization cannot identify10s response')
    if lid == 'L02':
        if q.get('quantity') != 'aged_cFDP_benchmark' or q.get('geometry_class') != 'cantilever_FDP' or q.get('aging_cycles') != 1200000:
            return unknown('Quantity/geometry/history mismatch')
        r = next((r for r in data[lid] if r['group'] == q.get('group')), None)
        if r is None:
            return unknown('No source group')
        if q.get('asserted_early_failures', r['early_failures']) != r['early_failures']:
            return unknown('Aging failures erased or corrupted')
        return {'status': 'SOURCE_OBSERVATION', 'group': r['group'], 'surviving_fracture_mean_N': r['mean_N'], 'SD_N': r['SD_N'], 'n_initial': r['initial_n'], 'n_surviving': r['survivor_n'], 'observed_aging_failures': r['early_failures'], 'mean_only_screen_ge350': r['mean_N'] >= 350, 'history_and_mean_screen': r['early_failures'] == 0 and r['mean_N'] >= 350, 'new_specimen_risk': None, 'clinical_prediction': False, 'source_locator': r['source_locator'], 'resolution_level': 'POPULATION'}
    if lid == 'L03':
        if q.get('quantity') != 'post_cementation_vertical_marginal_gap' or q.get('cementation_force_N') != 20 or q.get('cementation_duration_s') != 60 or (q.get('state') != 'cemented'):
            return unknown('Measurement state/force/duration mismatch')
        r = next((r for r in data[lid] if r['group'] == q.get('group') and r['site'] == q.get('site')), None)
        if r is None:
            return unknown('No exact source site')
        if q.get('n_independent_copings', 8) != 8:
            return unknown('Four repeated sites cannot become32 independent copings')
        return {'status': 'SOURCE_OBSERVATION', 'mean_gap_um': r['mean_um'], 'SD_um': r['SD_um'], 'n_copings': 8, 'site': r['site'], 'new_coping_maximum_gap': None, 'per_coping_site_covariance': 'UNKNOWN', 'source_locator': r['source_locator'], 'resolution_level': 'PER_SURFACE_REGION'}
    if lid == 'L04':
        if q.get('quantity') != 'nominal_normalized_ISO_survival_observation' or q.get('angle_deg') != 30 or q.get('R') != 0.1 or (q.get('horizon_cycles') != 5000000) or (q.get('frequency_Hz') != 15):
            return unknown('Exact normalized benchmark protocol required')
        if q.get('design') not in data[lid]['designs']:
            return unknown('Design outside source cohort')
        if q.get('absolute_force_N') is not None:
            return unknown('Per-design absolute monotonic Fult missing; nominal levels were rounded10N')
        if q.get('use_printed_fit', False):
            return unknown('Printed coefficients and endpoint mutually inconsistent; no silent refit')
        frac = q.get('nominal_load_fraction')
        if frac == 0.1:
            return {'status': 'SOURCE_OBSERVATION', 'outcome': 'Four_of_four_runout_at_nominal_10pct_Fult', 'failure_count': 0, 'n': 4, 'group': q['design'], 'load_rounding_N': 10, 'population_failure_probability_upper95pct': r1[lid]['failure_probability_upper95pct_at_four_of_four_runout'], 'binomial_iid_closure': True, 'general_strength_certificate': False, 'source_locator': read('SOURCE_CONTRACTS.json')[lid]['source_path'] + '#/article/body/sec[3]/sec[2]/sec[2]/p', 'resolution_level': 'POPULATION'}
        if frac == 0.15:
            return {'status': 'SOURCE_OBSERVATION', 'outcome': 'At_least_one_failure_at_nominal_15pct_Fult', 'failure_count_exact': None, 'n': 4, 'general_strength_certificate': False, 'resolution_level': 'POPULATION'}
        return unknown('Unmeasured nominal load class; bracket is not an interpolation law')
    if lid == 'L05':
        if q.get('quantity') == 'full_luting_temperature':
            return unknown('Exothermic heat/contact/whole-tooth geometry unknown; radiant transmission alone cannot determine total temperature')
        if q.get('quantity') != 'fixed_closure_slab_radiation_only' or q.get('dentin_mm') != 1 or q.get('ceramic_mm') != 3.5 or (q.get('duration_s') != 40):
            return unknown('Only stated fixed-property1mm slab40s envelope implemented')
        return {'status': 'MODEL_CONDITIONAL_RANGE', 'deltaT_C_range': r1[lid]['native_slab_eta0to1_deltaT_C_range'], 'measured_no_luting_deltaT_C': 4.1, 'closure': 'fixed-property,adiabatic slab;absorption in[0,1]', 'floating_and_continuum_rigorous_enclosure': 'MISSING inR1', 'clinical_prediction': False, 'resolution_level': 'PER_TOOTH'}
    return unknown('Unknown link')
