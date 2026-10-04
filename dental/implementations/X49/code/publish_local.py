"""Generate local reader artifacts from saved numerical output; no publishing."""
import html
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parents[1]

def main():
    active = json.loads((HERE / 'ACTIVE_ROUND.json').read_text())
    r = json.loads((HERE / ('rounds/' + active['round'] + '.json')).read_text())
    (HERE / 'results.json').write_text(json.dumps(r, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    run = Path(r['run_directory'])
    rows = []
    for case in r['cases']:
        report = json.loads((run / (case['case'] + '.json')).read_text())
        (HERE / 'reports' / (case['case'] + '.json')).write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
        for (name, rule) in report['rules'].items():
            witness = rule.get('witness', rule.get('witnesses', []))
            rows.append('<tr><td>' + html.escape(case['case']) + '</td><td>' + html.escape(name) + '</td><td>' + rule['status'] + '</td><td>' + html.escape(rule.get('reason', '')) + '</td><td><pre>' + html.escape(json.dumps(witness, ensure_ascii=False)) + '</pre></td></tr>')
    doc = '<!doctype html><html lang="en"><meta charset="utf-8"><title>X49 design gate</title><style>body{font:15px system-ui;margin:24px;max-width:1400px}table{border-collapse:collapse;width:100%}td,th{padding:8px;border:1px solid #ccc;vertical-align:top}pre{white-space:pre-wrap;max-width:450px}</style><h1>Independent geometry gate</h1><p>Digital surface decisions with explicit scope. Physical fit, cement, machine motion and clinical eligibility remain UNKNOWN.</p><p>' + str(r['passed_tests']) + '/' + str(r['total_tests']) + ' frozen controls. UNKNOWN slots: ' + str(r['rule_counts']['UNKNOWN']) + '/' + str(sum(r['rule_counts'].values())) + '.</p><img width="100%" src="../figures/design_gate.png"><table><tr><th>Case</th><th>Rule</th><th>Status</th><th>Reason</th><th>Surface witness</th></tr>' + ''.join(rows) + '</table></html>'
    (HERE / 'reports/index.html').write_text(doc)
    print('Saved results.json, reports/index.html and per-case JSON')
if __name__ == '__main__':
    main()
