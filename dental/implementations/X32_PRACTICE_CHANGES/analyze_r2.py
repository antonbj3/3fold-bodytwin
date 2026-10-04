import json, copy, hashlib, xml.etree.ElementTree as E, csv
from pathlib import Path
from extract_reviews import extract
R = Path(__file__).resolve().parent
EXPECTED_CAUSES = ['Loss of retention', 'Tooth loss', 'Tooth fracture', 'Replace old crown by new one due to change situation', 'Recurrent caries', 'Periapical destruction', 'Excessive bone loss', 'Major chipping', 'Mobility', 'Loss of tooth vitality', 'Esthetic dysfunction', 'Framework fracture', 'Tooth loss unknown reason', 'Poor fitting']

def gate(rows, total, source_quantity):
    assert source_quantity == 'Reasons for the single crown failures', 'wrong endpoint'
    assert [x['cause'] for x in rows] == EXPECTED_CAUSES, 'cause identity mismatch'
    assert sum((x['n'] for x in rows)) == total == 230, 'denominator mismatch'
    assert all((x['n'] >= 0 for x in rows)), 'negative count'
    return True

def main():
    tree = E.parse(R / 'sources/PMC9546353.xml').getroot()
    t = tree.find(".//table-wrap[@id='eos12871-tbl-0005']")
    rows = []
    total = None
    caption = ''.join(t.find('caption').itertext()).strip()
    for tr in t.findall('.//tbody/tr'):
        a = [''.join(x.itertext()).strip() for x in tr]
        if a[0] == 'Total':
            total = int(a[1])
            continue
        (cause, n, percentage) = a
        family = 'X1b_framework_precursor' if cause == 'Framework fracture' else 'X13_X14_fit_precursor' if cause == 'Poor fitting' else 'UNMODELLED_CLINICAL_FAILURE'
        rows.append(dict(cause=cause, n=int(n), reported_pct=float(percentage), mapped_family=family, resolution_level='POPULATION', exact_mechanism_validation='UNKNOWN', source='https://pmc.ncbi.nlm.nih.gov/articles/PMC9546353/#eos12871-tbl-0005'))
    assert gate(rows, total, caption)
    faults = {}
    for (name, mutator) in [('wrong_count', lambda a: a[0].update(n=a[0]['n'] + 1)), ('tooth_as_crown_fracture', lambda a: a[2].update(cause='Framework fracture'))]:
        bad = copy.deepcopy(rows)
        mutator(bad)
        try:
            gate(bad, total, caption)
        except AssertionError:
            faults[name + '_rejected'] = True
    try:
        gate(rows, total, 'Biological and technical complications for single crowns')
    except AssertionError:
        faults['overlapping_complications_as_disjoint_failures_rejected'] = True

    def endpoint_gate(fraction, validated):
        assert fraction is None or validated, 'No matched clinical calibration'
    try:
        endpoint_gate(4 / 230, False)
    except AssertionError:
        faults['family_share_as_clinical_coverage_rejected'] = True
    assert all(faults.values()) and len(faults) == 4
    with (R / 'raw/failure_causes.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    freq = extract()
    framework = sum((x['n'] for x in rows if x['cause'] == 'Framework fracture'))
    fit = sum((x['n'] for x in rows if x['cause'] == 'Poor fitting'))
    out = {'round': 'R2', 'claim_type': 'information_link', 'external_referent': {'kind': 'independent_measurement', 'locator': rows[0]['source'], 'compared_quantity': 'Disjoint attributed causes of 230 failed crowns; 1037 crowns/401 patients, failures in149 patients; mean followup134.8 months', 'refutes_us': True}, 'cohort': {'all_crowns': 1037, 'all_patients': 401, 'failed_crowns': total, 'patients_with_failed_crowns': 149, 'followup_mean_months': 134.8, 'material_scope': 'historical mixed materials; not modern monolithic zirconia or implants', 'resolution_level': 'POPULATION'}, 'X1b_family_matched_failures': framework, 'X1b_family_share_percent': 100 * framework / total, 'X1b_exact_mechanism_mapping_bound_percent': [0, 100 * framework / total], 'X1b_plus_fit_precursor_family_upper_percent': 100 * (framework + fit) / total, 'validated_clinical_failure_fraction': None, 'validated_clinical_endpoint_models': 0, 'sampling_interval': 'UNKNOWN: cause-to-patient assignments unavailable, so no patient bootstrap or iid interval', 'mapping_bound_type': 'finite-cohort epistemic mapping bound, not confidence interval', 'review_frequency_cells': len(freq), 'review_exclusions': {'candidate_cells': 80, 'retained': 40, 'excluded': 40, 'excluded_fraction': 0.5, 'reasons': {'other_materials_outside_declared_table_subset': 30, 'Sailer_zirconia_material_reclassification_correction_fulltext_unverified': 10}}, 'fault_injections': faults, 'attribution_control': 'Independent table-row sum =230 and exact framework count =4; no method contest', 'limits': ['Crown frame fracture label lacks intaglio tensile initiation and material match', 'Failure reason attribution retrospective; not independently examined at follow-up', 'Different review denominators prevent inferred unique failure fraction', 'Mechanical group mean/fit precursor does not predict retention, caries, chipping, tooth loss or biological events']}
    (R / 'rounds/RESULTS_R2.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: out[k] for k in ['X1b_family_matched_failures', 'X1b_family_share_percent', 'X1b_exact_mechanism_mapping_bound_percent', 'X1b_plus_fit_precursor_family_upper_percent', 'review_frequency_cells', 'fault_injections']}, indent=2))
if __name__ == '__main__':
    main()
