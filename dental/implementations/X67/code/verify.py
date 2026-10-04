import numpy as np, re, copy, xml.etree.ElementTree as ET
from common import *
from r4 import check
from bench_calibrate import check_metadata

def main():
    for k in ['R1', 'R2', 'R3', 'R4']:
        verify('PREREG_' + k + '.json')
    for n in ['FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R2.json', 'FROZEN_PREDICTIONS_R3.json']:
        verify(n)
    for d in load('FROZEN_PREDICTIONS.json')['predictions']:
        assert sha(ROOT / d['path']) == d['sha256']
    sources = load('raw/source_manifest.json')
    for d in sources:
        assert sha(d['path']) == d['sha256']
    rec = load('raw/measurements.json')
    roots = {d['study']: ET.parse(d['path']).getroot() for d in sources}
    for d in rec:
        (row, col) = [int(v) for v in re.findall('\\[(\\d+)\\]', d['table_locator'])]
        tab = next((e for e in roots[d['study']].findall('.//table-wrap') if e.get('id') == d['table_id']))
        cell = list(tab.findall('.//tr')[row - 1])[col - 1]
        assert d['raw_cell'] == ' '.join(''.join(cell.itertext()).split())
    p = load('FROZEN_PREDICTIONS_R2.json')
    assert p['basis_sha256'] == sha(ROOT / 'raw/source_basis.npz')
    assert p['dataset_sha256'] == sha(ROOT / 'raw/R2_dataset.json')
    r1 = load('attempts/R1_results.json')
    r2 = load('attempts/R2_results.json')
    r3 = load('attempts/R3_results.json')
    r4 = load('attempts/R4_results.json')
    assert r1['sufficiency']['identity_error_J'] == 0
    assert all((v['identity_error_C'] == 0 for v in r3['sufficiency'].values() if isinstance(v, dict) and 'identity_error_C' in v))
    assert r1['implied_eta_for_scalar_fit'] > 1 and r1['faults']['primary_plus10_C_rejected']
    assert not r2['transfer_gate_pass'] and (not r2['physical_primary_evidence_gate_pass'])
    assert r2['faults']['no_heldout_leakage'] and r2['faults']['all_injected_plus10_rejected']
    assert all((d['unbounded_control_fault_rejected'] for d in r2['folds']))
    assert r3['technical_gate_pass'] and r3['independent_physical_measurements'] == 0
    assert all((v['rejected'] for v in r3['faults'].values()))
    certificates = load('raw/R4_certificates.json')
    for d in certificates:
        A = np.array(d['A_ub'])
        b = np.array(d['b_ub'])
        for n in ['lower_certificate', 'negative_upper_certificate']:
            assert check(A, b, d[n])
        wrong = copy.deepcopy(d['lower_certificate'])
        from fractions import Fraction
        wrong['lower_bound_rational'] = str(Fraction(wrong['lower_bound_rational']) + Fraction(1, 10))
        assert not check(A, b, wrong)
    from bench_calibrate import read_csv
    src = (ROOT / 'raw/heater_fixture.csv').read_text()
    fault = ROOT / 'raw/units_fault.csv'
    fault.write_text(src.replace('time_s', 'time_ms', 1))
    try:
        read_csv(fault)
        raise AssertionError('unit fault accepted')
    except ValueError:
        pass
    meta = load('raw/heater_meta.json')
    bad = load('raw/drill_meta.json')
    bad['mechanical_work_provenance'] = 'motor_limit'
    try:
        check_metadata(meta, bad)
        raise AssertionError('motor limit used as measured work')
    except ValueError:
        pass
    result = load('results.json')
    assert result['review_state'] == 'PENDING_INDEPENDENT_REVIEW' and isinstance(result['external_referent']['refutes_us'], bool)
    for d in result['artifact_manifest']:
        assert sha(ROOT / d['path']) == d['sha256']
    assert sha(ROOT / 'results.json') == load('GRAPH_FEEDBACK.json')['sha256']
    code_manifest = [dict(path=str(p.relative_to(ROOT)), sha256=sha(p), bytes=p.stat().st_size) for p in sorted((ROOT / 'code').glob('*.py'))]
    code_manifest.append(dict(path='run_all.sh', sha256=sha(ROOT / 'run_all.sh'), bytes=(ROOT / 'run_all.sh').stat().st_size))
    dump('CODE_MANIFEST.json', code_manifest)
    dump('VERIFY.json', dict(status='PASS', timestamp_utc=now(), checks=['four immutable PREREGs', 'frozen prediction bytes and retrospective disclosure', 'all primary table cell locators and source hashes', 'basis/dataset freeze', 'exact summary identities', 'actual fault injections and measured-work provenance', 'all exact rational LP dual certificates', 'result/feedback binding', 'artifact and code manifest'], scientific_validation='UNKNOWN; no physical measurements', certificates_verified=2 * len(certificates), known_procedural_limit='Only shared R1 decomposition predates all numerics; dedicated mathematical leaf expansion saved after R2-R4, with their operation/assumption contracts frozen beforehand.'))
    state('DEMO_VERIFIED', 'One-command verification PASS; temperature/dose/benklass physical claims UNKNOWN', 'Register pending graph feedback and leave next actual measurement construction')
    print('VERIFY PASS: sources, freezes, scientific negative gates, contracts and 12 exact LP certificates')
if __name__ == '__main__':
    main()
