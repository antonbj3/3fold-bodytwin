"""One-command cache-aware reconstruction, preserving every frozen prediction/gate."""
import sys, os, time, json, subprocess, datetime
from pathlib import Path
P = Path(__file__).resolve().parents[1]
os.chdir(P)
os.environ['PYTHONPATH'] = str(P / 'code')
ledger = P / 'raw/RUN_LEDGER.jsonl'
ledger.parent.mkdir(exist_ok=True)

def run(args, log, system=False):
    env = os.environ.copy()
    if system:
        env['PYTHONNOUSERSITE'] = '1'
    st = time.perf_counter()
    with open(P / 'raw' / log, 'a') as fh:
        r = subprocess.run(args, stdout=fh, stderr=subprocess.STDOUT, env=env)
    with ledger.open('a') as fh:
        fh.write(json.dumps(dict(command=args, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), wall_s=time.perf_counter() - st, exit_code=r.returncode, log='raw/' + log)) + '\n')
    print(' '.join(args), 'exit', r.returncode, 's', round(time.perf_counter() - st, 2), flush=True)
    if r.returncode:
        raise SystemExit(r.returncode)
py = sys.executable
if not (P / 'PREREG_R1.json').exists():
    run([py, 'code/preregister.py'], 'preregister.log')
if len(list((P / 'raw/cases').glob('*.json'))) < 993:
    run([py, 'code/contact.py'], 'batch_R1.log')
if not (P / 'FROZEN_PREDICTIONS.json').exists():
    run([py, 'code/report_capability.py', 'fit', '--round', 'R1'], 'fit_R1.log')
if not (P / 'raw/RESULTS_R1.json').exists():
    run([py, 'raw/failed_checks/report_capability_R1_frozen.py', 'evaluate', '--round', 'R1'], 'evaluate_R1.log')
if not (P / 'PREREG_R2.json').exists():
    run([py, 'code/preregister_r2.py'], 'preregister.log')
if not (P / 'raw/OPPOSITION_FEATURES_R2.json').exists():
    run([py, 'code/opposition.py', 'extract'], 'opposition_R2.log')
if not (P / 'FROZEN_PREDICTIONS_R2.json').exists():
    run([py, 'code/report_capability.py', 'fit', '--round', 'R2'], 'fit_R2.log')
if not (P / 'raw/RESULTS_R2.json').exists():
    run([py, 'code/report_capability.py', 'evaluate', '--round', 'R2'], 'evaluate_R2.log')
if not (P / 'raw/RESULTS_TOOTHWISE_R2.json').exists():
    run([py, 'code/opposition.py', 'evaluate'], 'toothwise_R2.log')
if not (P / 'PREREG_R3.json').exists():
    run([py, 'code/preregister_r3.py'], 'preregister.log')
if not (P / 'raw/RESULTS_R3.json').exists():
    run([py, 'code/certificates.py'], 'certificates_R3.log')
if not (P / 'PREREG_R4.json').exists():
    run([py, 'code/preregister_r4.py'], 'preregister.log')
if not (P / 'raw/RESULTS_R4.json').exists():
    run([py, 'code/refine_certificates.py'], 'refine_R4.log')
if not (P / 'PREREG_R5.json').exists():
    run([py, 'code/preregister_r5.py'], 'preregister.log')
if not (P / 'raw/CORRECTED_REPORTS_R5.json').exists():
    run([py, 'code/repair_reports.py', 'extract'], 'repair_R5.log')
if not (P / 'FROZEN_PREDICTIONS_R5.json').exists():
    run([py, 'code/repair_reports.py', 'fit_freeze'], 'repair_R5.log')
if not (P / 'raw/RESULTS_R5.json').exists():
    run([py, 'code/repair_reports.py', 'evaluate'], 'evaluate_R5.log')
if not (P / 'PREREG_R6.json').exists():
    run([py, 'code/scoped_reports.py', 'prereg'], 'scopes_R6.log')
if not (P / 'raw/RESULTS_R6.json').exists():
    run([py, 'code/scoped_reports.py', 'run'], 'scopes_R6.log')
if not (P / 'FROZEN_MEASUREMENT_PANEL.json').exists():
    run([py, 'code/compare_measurement.py', '--prepare'], 'measurement_prepare.log')
run([py, 'code/verify.py'], 'verification.log')
run([py, 'code/package_results.py'], 'package.log')
run([py, 'code/extend_package.py'], 'package_extended.log')
run(['/usr/bin/python3', 'code/figure.py'], 'figure.log', True)
run(['/usr/bin/python3', 'code/certificate_figure.py'], 'figure.log', True)
c = json.load(open(P / 'raw/INPUT_MANIFEST.json'))['numerical_panel'][0]
run([py, 'code/demo_case.py', '--case', c], 'demo.log')
run([py, 'code/final_docs.py'], 'final_docs.log')
print('Complete: README_DEMO.md, results.json, figures/, raw/TOOTH_PAIR_MAP.csv. PENDING_INDEPENDENT_REVIEW.', flush=True)
