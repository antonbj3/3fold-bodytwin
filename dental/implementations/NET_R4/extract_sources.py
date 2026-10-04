"""Replay primary rows from audited local snapshots; no network or fitting."""
import re, xml.etree.ElementTree as ET
from decimal import Decimal as D
from common import *

def txt(el):
    return ' '.join(''.join(el.itertext()).split())

def source_ref(k, locator, doi):
    r = entry(k)
    return {'source_id': k, 'source_file': r['original_path'], 'source_sha256': r['sha256'], 'snapshot': r['snapshot'], 'locator': locator, 'url': doi}

def extract(write_output=True):
    verify()
    s = source_path('sagheb').read_text()
    m = re.search('initially\\s+([\\d.]+)\\s*N\\s*±\\s*([\\d.]+)\\s*\\(range ([\\d.]+) to ([\\d.]+)\\) to ([\\d.]+)\\s*N\\s*±\\s*([\\d.]+)\\s*\\(range ([\\d.]+) to ([\\d.]+)\\)', s)
    assert m, 'SAGHEB_ABSTRACT_LAYOUT'
    x = list(m.groups())
    delta = D(x[0]) - D(x[4])
    sag = {'id': 'SAGHEB_PRELOAD', 'edge_id': 'D-E-PRELOAD-HISTORY', 'variable': 'preload', 'claim_type': 'information_link', 'quantity': 'axial_preload_first_minus_tenth_group_mean', 'unit': 'N', 'resolution_level': 'POPULATION', 'anatomical_support_level': 'PER_TOOTH', 'timescale': 'HANDOVER', 'reference': source_ref('sagheb', 'Abstract Results; Methods paragraph: 25 unused units, 25 Ncm, ten tightenings, hold 60 s; Figure 5 (Table 1 does NOT contain the endpoint means)', 'https://doi.org/10.1186/s40729-023-00473-3'), 'audit_corrections': ['C011', 'C013', 'C016'], 'population': {'n': 25, 'unit': 'implant-abutment-screw complexes', 'system': 'Nobel Replace Select Tapered / Temporary Abutment Non-engaging / carbon-coated screw'}, 'protocol': {'id': 'SAGHEB_2023_DRY_25NCM', 'torque_Ncm': 25, 'repetitions': 10, 'hold_seconds': 60, 'lubrication': 'none; dry', 'sampling_interval_s': '0.022', 'comparison_states': [1, 10]}, 'reported': {'first_mean_N': x[0], 'first_SD_N': x[1], 'first_range_N': x[2:4], 'tenth_mean_N': x[4], 'tenth_SD_N': x[5], 'tenth_range_N': x[6:8]}, 'derived': {'difference_N': str(delta), 'relative_loss_percent': str(delta / D(x[0]) * 100)}, 'uncertainty': 'Group SD and endpoint observed ranges; paired difference SD/covariance absent. No individual loss bound or target-system prediction.', 'supports': ['published_group_preload_change'], 'does_not_support': ['assembly_gap_prediction', 'friction_cause_identification', 'relaxation_after_tool_release', 'new_system_transfer'], 'full_edge_validated': False, 'refutes_us': True}
    c = source_path('cunali').read_text()
    table = c[c.index('Table 1. Means'):]
    replica = re.search('Replica\\s+77\\.90±12\\.21 a\\s+83\\.63±14\\.43 a\\s+162\\.65±39\\.88 b\\s+([\\d.]+)±([\\d.]+)', table)
    ct = re.search('Micro-CT\\s+68\\.73±8\\.86 a\\s+81\\.00±15\\.84 a\\s+110\\.83±14\\.08 a\\s+([\\d.]+)±([\\d.]+)', table)
    assert replica and ct, 'CUNALI_TABLE_LAYOUT'
    cun = {'id': 'CUNALI_INTERNAL_GAP', 'edge_id': 'D-E-REPLICA-CT', 'variable': 'precementation_internal_gap', 'claim_type': 'information_link', 'quantity': 'mid_occlusal_internal_gap_replica_minus_dry_CT_group_mean', 'unit': 'um', 'resolution_level': 'POPULATION', 'anatomical_support_level': 'PER_SURFACE_REGION', 'region': 'mid_occlusal', 'timescale': 'SIMULTANEOUS', 'observation_schedule': 'Sequential measurement sessions, nominal same pre-cementation restoration; not synchronized acquisition', 'reference': source_ref('cunali', 'Table 1, p470, Amann row, Mid-occlusal wall column; Methods pp468-469', 'https://doi.org/10.1590/0103-6440201601531'), 'audit_corrections': ['C017', 'C018', 'C019', 'C020'], 'population': {'n': 10, 'unit': 'Amann Ceramill Zi zirconia copings', 'design': 'lower first molar metal master; same copings examined with both methods'}, 'protocol': {'id': 'CUNALI_2017_AMANN_MO', 'replica': 'Light-body PVS; firm finger seating; 5 minutes before removal; four section readings averaged per region', 'micro_CT': 'Copings cleaned with alcohol and reseated without cement; SkyScan 1172; 100 kV, 100 uA, pixel 13 um', 'cementation': 'not performed'}, 'reported': {'replica_mean_um': replica[1], 'replica_SD_um': replica[2], 'ct_mean_um': ct[1], 'ct_SD_um': ct[2]}, 'derived': {'difference_um': str(D(replica[1]) - D(ct[1]))}, 'uncertainty': 'Marginal group SDs, not instrument noise; paired specimen differences and reseating repeatability unavailable.', 'supports': ['published_regional_method_state_contrast'], 'does_not_support': ['marginal_gap_prediction', 'cemented_film_prediction', 'instrument_noise_SD', 'individual_method_agreement'], 'full_edge_validated': False, 'refutes_us': True}
    root = ET.parse(source_path('choi'))
    t2 = root.find('.//table-wrap[@id="materials-12-03264-t002"]')
    static = []
    for row in t2.findall('.//tbody/tr')[:3]:
        cells = [txt(el) for el in row.findall('td')]
        static.append({'TX': int(cells[1].split()[0]), 'EV': int(cells[3].split()[0])})
    t3 = root.find('.//table-wrap[@id="materials-12-03264-t003"]')
    fatigue = []
    group = None
    for row in t3.findall('.//tr'):
        cells = [txt(el) for el in row.findall('td')]
        if len(cells) == 1 and cells[0].startswith(('TX', 'EV')):
            group = cells[0][:2]
        if len(cells) == 4 and cells[0].isdigit():
            (low, high) = map(int, cells[1].split('–'))
            cycles = [int(v.strip().replace(',', '')) for v in cells[2].split(';')]
            fatigue.append({'group': group, 'loading_level_percent': int(cells[0]), 'min_N': low, 'max_N': high, 'cycles': cycles, 'events': [v < 5000000 for v in cycles], 'censoring': '5,000,000 is right-censored runout, not failure', 'resolution_level': 'PER_TOOTH', 'unit_of_analysis': 'implant specimen (not a patient tooth)', 'table_locator': f'Table 3/{group}/{cells[0]}%/three reported specimens'})
    cho = {'id': 'CHOI_FATIGUE', 'edge_id': 'D-E-K33', 'variable': 'fatigue_load', 'claim_type': 'information_link', 'quantity': 'protocol_specific_load_cycles_failure_or_runout', 'unit': 'N', 'resolution_level': 'PER_TOOTH', 'timescale': 'HANDOVER', 'reference': source_ref('choi', 'Tables 2 and 3; Methods: 30 degrees, 11-mm center height, R=0.1, 15 Hz, 5 million cycles', 'https://doi.org/10.3390/ma12193264'), 'audit_corrections': ['C027', 'C028', 'C029', 'C030'], 'population': {'static_n_per_system': 3, 'fatigue_n_per_system': 12, 'systems': {'TX': 'OsseoSpeed TX diameter4.0 mm length11 mm', 'EV': 'OsseoSpeed EV diameter4.2 mm length11 mm'}, 'material': 'commercially pure grade4 Ti fixtures and abutments'}, 'protocol': {'id': 'CHOI_2019_TX_EV', 'load_angle_deg': 30, 'loading_center_height_mm': 11, 'R': 0.1, 'frequency_Hz': 15, 'runout_cycles': 5000000, 'static_loading_speed_mm_min': 1, 'standard_reported': 'ISO 14801:2013; historical source protocol, not a current standards recommendation'}, 'reported': {'static_individual_N': static, 'static_mean_N': {'TX': 711, 'EV': 791}, 'static_SD_N': {'TX': 36, 'EV': 58}, 'fatigue_rows': fatigue}, 'derived': {'nominal_40percent_static_N': {'TX': '284.4', 'EV': '316.4'}}, 'uncertainty': '3 specimens per tested load/system; runout is censoring; no population fatigue limit or interpolation between tested loads identified. Static and fatigue specimens are different.', 'supports': ['published_specimen_load_cycle_outcomes'], 'does_not_support': ['universal_40percent_rule', 'measured_root_defect_to_fatigue_law', 'same_specimen_static_fatigue_pairs', 'new_system_transfer'], 'full_edge_validated': False, 'refutes_us': True}
    assert len(fatigue) == 8 and len(static) == 3
    out = {'schema': 'audited-primary-observations-r4', 'claim_type': 'information_link', 'records': [sag, cun, cho], 'source_status': 'Published observations already existed; R4 replays and connects them. No newly measured physical data.'}
    sag['protocol']['fixture_separation_mm'] = '0.10'
    sag['protocol']['fixture_condition'] = 'Implant and abutment held out of contact for load-cell measurement.'
    cho['protocol'].pop('loading_center_height_mm', None)
    cho['protocol']['loading_center_height_mm'] = 'UNKNOWN: Figure 2 dimensions not extracted from local XML'
    cho['reference']['locator'] = 'Tables 2 and 3; Methods: 30 degrees, R=0.1, 15 Hz, 5 million cycles. Fixture lever dimension not verified.'
    cho['protocol']['environment'] = 'atmospheric; 20 degC +/-5 per Methods'
    if write_output:
        write('SOURCE_OBSERVATIONS.json', out)
    return out
if __name__ == '__main__':
    r = extract()
    print([(x['id'], x.get('derived')) for x in r['records']])
