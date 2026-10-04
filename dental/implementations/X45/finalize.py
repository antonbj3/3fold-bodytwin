"""Package exact scientific outputs; timing can vary while prediction hashes stay frozen."""
import datetime as dt, hashlib, json, platform, re, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(p, o):
    p.write_text(json.dumps(o, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def strip_timing(o):
    if isinstance(o, dict):
        return {k: strip_timing(v) for (k, v) in o.items() if k not in ['kernel_wall_seconds', 'cost']}
    if isinstance(o, list):
        return [strip_timing(v) for v in o]
    return o

def main():
    r = json.loads((ROOT / 'results.json').read_text())
    coupon = json.loads((ROOT / 'raw/COUPON_CONSUMER.json').read_text())
    coupon['external_referent']['refutes_us'] = False
    r['coupon_consumer'] = coupon
    r['control_repair'] = json.loads((ROOT / 'raw/CONTROL_REPAIR.json').read_text())
    r['external_referent'].append(coupon['external_referent'])
    r['artifact_resolution_contract'] = 'Every empirical edge is POPULATION; field/export geometry is PER_POINT/PER_SURFACE_REGION. PHENOMENOLOGICAL debt is explicit in CHAIN_PORTS.json.'
    r['current_capability'] = 'Two bounded history chains to explicit research consumers, with no full Northstar EXIST promotion'
    r['retained_control_weakness'] = {'construction': 'K34_WEAR_R1', 'undetected_source_mutations': 1, 'tested': 15, 'reason': 'corrupted measurement can move toward an already wrong model', 'replacement_controls': '52/52 specified returned-port faults detected,30 exact-rational enclosure checks'}
    r['attrition'].update({'ISQ_R1_rejected': 6, 'ISQ_R1_tested': 8, 'wear_R1_rejected': 15, 'wear_R1_tested': 15, 'wear_R2_rejected': 7, 'wear_R2_tested': 15, 'source_plus10mg_mutation_undetected': 1, 'source_plus10mg_mutations': 15})
    r['cost']['wall_seconds_full_shell_replay'] = 'Individual subprocess timings in ' + sys.argv[1] + '; startup, sources, controls and figure included there. Physical fallback NOT_RUN.'
    code = [p for p in ROOT.glob('*.py')] + [ROOT / 'run_all.sh']
    r['code_sha256'] = {p.name: sha(p) for p in sorted(code)}
    r['frozen_prediction_files'] = {p.name: sha(p) for p in sorted(ROOT.glob('FROZEN_PREDICTIONS_*.json'))}
    r['preregistration_files'] = {p.name: sha(p) for p in sorted(ROOT.glob('PREREG_*.json'))}
    r['artifact_bytes'] = {str(p.relative_to(ROOT)): p.stat().st_size for p in ROOT.rglob('*') if p.is_file() and 'replays' not in p.parts}
    r['large_arrays'] = []
    r['source_screening_file'] = 'SOURCE_SCREENING.json'
    r['plot'] = {'path': 'figures/chains.png', 'sha256': sha(ROOT / 'figures/chains.png')}
    r['runtime_versions'] = {'python': platform.python_version(), 'numpy': __import__('numpy').__version__, 'scipy': __import__('scipy').__version__, 'mpmath': __import__('mpmath').__version__, 'figure_runtime': '/usr/bin/python3 -s; distro NumPy/matplotlib avoids user NumPy2 ABI mismatch'}
    write(ROOT / 'results.json', r)
    payload = strip_timing({k: r[k] for k in ['schema', 'review_state', 'claim_type', 'external_referent', 'constructed_chains', 'full_Northstar_chain_promotions', 'attempts', 'sufficiency_tests', 'source_manifest', 'parser_controls', 'coupon_consumer', 'control_repair', 'retained_control_weakness', 'code_sha256', 'frozen_prediction_files', 'preregistration_files', 'attrition']})
    snapshot = ROOT / 'RESULTS_FOR_REVIEW.json'
    if snapshot.exists():
        old = json.loads(snapshot.read_text())
        if old != payload:
            raise ValueError('scientific review snapshot drift; preserve this failure and version the changed construction')
    else:
        write(snapshot, payload)
    manifest = {'primary_inputs': r['source_manifest'], 'own_code': r['code_sha256'], 'scientific_snapshot': {'path': str(snapshot), 'sha256': sha(snapshot)}, 'frozen_predictions': r['frozen_prediction_files'], 'local_artifacts': {str(p.relative_to(ROOT)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(ROOT.glob('raw/*.json'))}, 'no_large_arrays': True, 'intermediate_limit_bytes': 3000000000}
    write(ROOT / 'SOURCE_MANIFEST.json', manifest)
    if r['attempts']['K31_R1']['passes'] != 2 or r['attempts']['K34_WEAR_R1']['passes'] != 0 or r['attempts']['K34_WEAR_R2']['passes'] != 8:
        raise ValueError('negative predecessor drift')
    if r['attempts']['K31_R2']['passes'] != 6 or r['attempts']['K34_WEAR_R3']['passes'] != 15:
        raise ValueError('successor gate failure')
    if not all((c['pass'] for c in r['attempts']['K34_LTD_R1']['correctness_control'])):
        raise ValueError('LTD arithmetic failed')
    if not r['control_repair']['all_detected'] or not r['control_repair']['all_enclosed']:
        raise ValueError('repaired control failed')
    current = json.loads((ROOT / 'CURRENT_WORK_STATE.json').read_text())
    current.update(status='ROUND_COMPLETE_PENDING_INDEPENDENT_REVIEW', updated_utc=dt.datetime.now(dt.timezone.utc).isoformat(), latest_gate={'K31': '6/6 held observation means; BIC/stiffness REFUSED', 'K34': '8 exposure states;15/15 conditional wear brackets; rate closures REFUTED', 'controls': '52/52 specified new faults caught; original1/15 weak mutation retained'}, next_operation='Independent review, then same-specimen temporal RFA/stiffness/BIC and matched wear-depth/volume acquisition. Full K16/K21/K48 remain data-blocked.')
    write(ROOT / 'CURRENT_WORK_STATE.json', current)
    print('Scientific snapshot', sha(snapshot), 'frozen predictions', len(r['frozen_prediction_files']))
if __name__ == '__main__':
    main()
