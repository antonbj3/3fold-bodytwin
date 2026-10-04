import html
import json

def html_report(out):
    esc = lambda x: html.escape(str(x))
    rows = []
    for (name, r) in out.get('rules', {}).items():
        measurement = {k: v for (k, v) in r.items() if k not in ('status', 'reason', 'resolution', 'manual_adjustment_minutes', 'manual_adjustment_source')}
        rows.append('<tr><td>' + esc(name) + '</td><td class="' + esc(r['status']) + '">' + esc(r['status']) + '</td><td>' + esc(r['reason']) + '<pre>' + esc(json.dumps(measurement, indent=2)) + '</pre></td><td>UNKNOWN</td></tr>')
    return '<!doctype html><html lang="en"><meta charset="utf-8"><title>STL design gate</title>\n<style>body{font:16px system-ui;max-width:1200px;margin:36px auto;padding:0 18px;color:#182e39}table{border-collapse:collapse;width:100%}td,th{padding:12px;border:1px solid #cbd5dc;text-align:left;vertical-align:top}pre{white-space:pre-wrap;max-width:680px;font-size:12px}.PASS{color:#18724b}.FAIL{color:#a42231}.UNKNOWN{color:#936215}</style>\n<h1>STL design gate — ' + esc(out['verdict']) + '</h1><p>Geometric decisions at the declared pose and scope. Complete restoration eligibility: UNKNOWN. Thresholds describe the named product IFU or the supplied lab protocol.</p>\n<table><tr><th>Rule</th><th>Decision</th><th>Measurement, point witnesses and scope</th><th>Estimated manual adjustment minutes per FAIL</th></tr>' + ''.join(rows) + '</table><details><summary>Input provenance and full JSON</summary><pre>' + esc(json.dumps(out, indent=2)) + '</pre></details></html>'
