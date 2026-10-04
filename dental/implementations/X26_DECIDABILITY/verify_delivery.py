import csv
import hashlib
import json
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from query_force import query
ROOT = Path(__file__).resolve().parent

def facit_ok(rows, tol):
    return len(rows) == 4 and all((abs(x['source_mean_um'] - x['export_mean_um']) <= tol and abs(x['source_sd_um'] - x['export_sd_um']) <= tol for x in rows))

def main():
    for sidecar in ROOT.glob('PREREG*.sha256'):
        data = sidecar.with_suffix('.json')
        assert hashlib.sha256(data.read_bytes()).hexdigest() == sidecar.read_text().strip(), str(data) + ' frozen hash drift'
    for sidecar in ROOT.glob('FROZEN_PREDICTIONS*.sha256'):
        data = sidecar.with_suffix('.json')
        assert hashlib.sha256(data.read_bytes()).hexdigest() == sidecar.read_text().strip(), str(data) + ' frozen prediction hash drift'
    pr = json.loads((ROOT / 'PREREG_DELIVERY.json').read_text())
    assert hashlib.sha256((ROOT / 'PREREG_DELIVERY.json').read_bytes()).hexdigest() == (ROOT / 'PREREG_DELIVERY.sha256').read_text().strip()
    m = json.loads((ROOT / 'EXTERNAL_MEASUREMENT_MANIFEST.json').read_text())
    assert hashlib.sha256((ROOT / m['local']).read_bytes()).hexdigest() == m['sha256']
    doc = ET.parse(ROOT / m['local']).getroot()
    table = next((t for t in doc.iter('table-wrap') if t.attrib.get('id') == 'TAB2'))
    with (ROOT / 'inputs/X13_facit.csv').open() as f:
        export = [x for x in csv.DictReader(f) if x['study'] == 'PMC10246932' and x['source_table'] == 'TAB2']
    checks = []
    for r in list(table.iter('tr'))[1:]:
        cells = [''.join(x.itertext()).strip() for x in r.findall('td')]
        a = next((x for x in export if x['arm'] == cells[0]))
        checks.append({'arm': cells[0], 'source_mean_um': float(cells[1]), 'source_sd_um': float(cells[2]), 'export_mean_um': float(a['measured_mean_um']), 'export_sd_um': float(a['reported_sd_um']), 'resolution_level': 'POPULATION', 'physical_point_bound': 'UNKNOWN', 'doi': m['doi'], 'table': 'TAB2'})
    tol = pr['metrics']['mean_and_sd_absolute_tolerance_um']
    bad = [dict(checks[0], export_mean_um=checks[0]['export_mean_um'] + 1)] + checks[1:]
    facit = {'external_referent': pr['external_referent'], 'rows': checks, 'export_matches_primary': facit_ok(checks, tol), 'injected_1um_rejected': not facit_ok(bad, tol), 'dropout': {'source_groups': 4, 'kept': 4, 'rejected': 0}, 'reported_means_equal_CAD_settings': all((abs(x['source_mean_um'] - float(x['arm'].split()[0][1:])) <= tol for x in checks)), 'restriction': 'These dry marginal group means are not local field bounds or a physical error floor for every crown.'}
    (ROOT / 'EXTERNAL_FACIT_CHECK.json').write_text(json.dumps(facit, indent=2) + '\n')
    assert facit['export_matches_primary'] and facit['injected_1um_rejected']
    empty_rejected = False
    try:
        query(ROOT / 'FORCE_MEASUREMENT_TEMPLATE.json')
    except ValueError as e:
        empty_rejected = 'NOT_RUN_MISSING_MEASUREMENT' in str(e)
    r3 = json.loads((ROOT / 'RESULTS_R3.json').read_text())
    r = r3['rows'][1]
    sample = {'case': r['case'], 'fdi': r['fdi'], 'basis_sha256': r['basis_sha256'], 'force_fractions': r['center'], 'epsilon_fraction': 0.2, 'total_force_N': 100, 'provenance_locator': 'our_own_fixture: X18 independent mixed-load test vector; not physical pressure'}
    file = ROOT / 'raw/pressure_port_own_fixture.json'
    file.write_text(json.dumps(sample, indent=2) + '\n')
    result = query(file)
    wrong = dict(sample, basis_sha256='0' * 64)
    file_bad = ROOT / 'raw/pressure_port_wrong_hash.json'
    file_bad.write_text(json.dumps(wrong, indent=2) + '\n')
    wrong_rejected = False
    try:
        query(file_bad)
    except ValueError:
        wrong_rejected = True
    report = {'external_facit': facit['export_matches_primary'], 'external_injection_rejected': facit['injected_1um_rejected'], 'empty_measurement_rejected': empty_rejected, 'wrong_basis_hash_rejected': wrong_rejected, 'pressure_fixture_conditional_upper_MPa': result['upper_MPa'], 'pressure_fixture_physical_status': result['physical_status'], 'fixture_status': 'our_own_fixture, not measured pressure', 'all_round_gates': True}
    for name in ['RESULTS_R2.json', 'RESULTS_R3.json', 'RESULTS_R4.json']:
        d = json.loads((ROOT / name).read_text())
        report['all_round_gates'] &= all(d['gates'].values()) and all(d['mutation_rejections'].values())
    (ROOT / 'VERIFICATION.json').write_text(json.dumps(report, indent=2) + '\n')
    assert empty_rejected and wrong_rejected and (result['physical_status'] == 'UNKNOWN') and report['all_round_gates']
    print('Verified primary measurement cells and rejecting force-measurement port')
if __name__ == '__main__':
    main()
