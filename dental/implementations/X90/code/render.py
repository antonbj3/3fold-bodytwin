"""Self-contained HTML: no scripts, remote fonts, assets or network requests."""
import base64
import html
import json
from pathlib import Path
from common import read, dump, sha, validate_decisions

def esc(value):
    return html.escape(str(value), quote=True)

def short(value):
    if value is None:
        return 'UNKNOWN'
    if isinstance(value, float):
        return f'{value:.6g}'
    if isinstance(value, dict):
        if 'FDI' in value:
            return 'FDI ' + str(value['FDI']) + ' · ' + value['design_verdict']
        if 'patient_margin_mm' in value:
            return 'Patient margin UNKNOWN · conditional population budgets below'
        return value.get('status', value.get('decision', 'See numerical evidence'))
    return str(value)

def render(out):
    out = Path(out)
    case = read(out / 'case.json')
    validate_decisions(case['decisions'])
    fig = base64.b64encode((out / 'case_figure.png').read_bytes()).decode()
    text = ['<!doctype html><html lang="en"><head><meta charset="utf-8">', '<meta name="viewport" content="width=device-width,initial-scale=1">', '<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; img-src data:; style-src \'unsafe-inline\'">', '<title>Dental forskningsfall ' + esc(case['patient']) + '</title>', '<style>body{font:16px/1.55 system-ui,sans-serif;color:#172a39;background:#f4f6f8;margin:0}main{max-width:1080px;margin:auto;padding:32px}h1{font-size:32px}h2{font-size:22px}.card{background:white;border:1px solid #d4dfe5;border-radius:10px;padding:20px;margin:18px 0}.status{font-weight:700;color:#9c4521}.value{font-size:21px}table{border-collapse:collapse;width:100%;font-size:14px}td,th{text-align:left;padding:10px;border-bottom:1px solid #dce3e7}img{width:100%;height:auto}pre{white-space:pre-wrap;overflow-wrap:anywhere;font:12px/1.4 monospace}code{overflow-wrap:anywhere}.muted{color:#4d6270}.tag{background:#e4edf1;padding:3px 7px;border-radius:4px;font-size:13px}@media print{body{background:white}main{padding:0}.card{break-inside:avoid}}</style></head><body><main>', '<h1>Fallrapport · ' + esc(case['patient']) + '</h1>', '<p> The delivered scan can now be followed to a research chandelier, design control and a test contract in a command. The original surface is geometric reference. Force , pulp and nerve require the measurements specified per decision. </p>', '<p class="status">Lab release: ABSTAIN · physical measurement: not performed · independent adapter review pending</p>', "<p> For a prosthetic researcher's attempt. Virtual preparation and predicted tooth owners. No clinical recommendations. </p>", '<img alt=" Measured tooth and generated research crown in the same coordinate frame" src="data:image/png;base 64,' + fig + '">']
    for row in case['decisions']:
        text += ['<section class="card" id="' + esc(row['id']) + '"><h2>' + esc(row['title']) + '</h2>', '<div class="value">' + esc(short(row['value'])) + ' <span class="muted">' + esc(row['unit']) + '</span></div>', '<p><span class="tag">' + esc(row['evidence']) + '</span> <span class="tag">' + esc(row['resolution']) + '</span> <span class="tag">' + esc(row['timescale']) + '</span> <strong>' + esc(row['status']) + '</strong></p>', '<p> <strong> Uncertainty : </strong> ' + esc(row['uncertainty']) + '</p>', '<p> <strong> What would change the decision: </strong> </p> <ul>']
        for c in row['would_change']:
            text.append('<li>' + esc(c['measurement']) + ' <span class="muted">(' + esc(c['edge']) + ', ' + esc(c['resolution']) + ', ' + esc(c['timescale']) + ')</span></li>')
        text.append('</ul>')
        if row['id'] == 'force':
            text.append('<table> <tr> <th> Predicted FDI </th> <th> Force interval [ N ] </th> <th> Status </th> </tr>')
            for t in row['per_tooth']:
                text.append('<tr><td>' + esc(t['fdi']) + '</td><td>' + esc(short(t['force_interval_N'])) + '</td><td>' + esc(t['evidence']) + '</td></tr>')
            text.append('</table>')
        if row['id'] == 'nerve':
            text.append('<table> <thead> <tr> <th> Guide </th> <th> Model Case </th> <th> Guide Error Budget [ mm ] </th> <th> Patient margin </th> </tr> </thead> <tbody>')
            for g in row['value']['population_guide_error_budgets']:
                text.append('<tr><td>' + esc(g['guide']) + '</td><td>' + esc(g['target']) + '</td><td>' + esc(f"{g['guide_budget_mm']:.4f}") + '</td><td>UNKNOWN</td></tr>')
            text.append("</tbody> </table> <p> Model target is conditional on the source's sample elements as population disclosure. Empirical physical coverage is unknown . </p>")
        text += ['<details> <summary> Full value , sources and measurement contract </summary> <pre>' + esc(json.dumps(row, ensure_ascii=False, sort_keys=True, indent=2)) + '</pre></details></section>']
    facit = case['external_referent']
    text += ['<section class="card"> <h2> Compared to the original surface </h2>', '<p> 8192 area distributed probes per surface. PER_TOOTH - summary of point-to-triangle distance. Scanner error and rigorous sampling/floating point enclosure missing. </p>', '<table> <tr> <th> Quantity [ mm ] </th> <th> Generated crown </th> <th> Direct rigid control </th> </tr>']
    for key in ['p95_mm', 'rms_mm', 'sampled_max_mm']:
        text.append('<tr><td>' + esc(key) + '</td><td data-metric="' + key + '">' + repr(facit['measured_original'][key]) + '</td><td>' + repr(facit['equally_informed_rigid_control'][key]) + '</td></tr>')
    text += ["</table> <p> The source's FDI labels here are X11 predictions, not independent annotation. The result of the check does not mean an algorithm advantage. </p> </section>", '<section class="card"> <h2> Frozen laboratory predictions </h2> <p> Frozen prior to the numeric original surface test. This is a retrospective reconstruction; the original archive can be read by the calculation process and no blind prediction test is claimed. Physical measurement has not been performed. Intended CAD - geometry can be measured against a manufactured surface; physical load, film and lifetime lack yet matched model. </p> <p> SHA256 : <code>' + esc(case['frozen_predictions_sha256']) + '</code></p><pre>' + esc((out / 'FROZEN_PREDICTIONS.json').read_text()) + '</pre></section>', '<section class="card"> <h2> Provenance and limitations </h2> <p> Bits2Bites : CC BY - NC - SA according to the dataset briefing. Local private research. Inherited X11 training basis has its own conditions. No new right to publishing is assumed. </p> <p> Digital geometry identity: <code>' + esc(case['source_geometry_sha256']) + '</code></p><pre>' + esc(json.dumps(case['source_members'], ensure_ascii=False, sort_keys=True, indent=2)) + '</pre>', '<details> <summary> Full machine-readable casework </summary> <pre id="case-json">' + esc(json.dumps(case, ensure_ascii=False, sort_keys=True, indent=2)) + '</pre></details></section>', '</main></body></html>']
    (out / 'report.html').write_text('\n'.join(text) + '\n')

def verify_report(out):
    import re
    out = Path(out)
    case = read(out / 'case.json')
    text = (out / 'report.html').read_text()
    for key in ['p95_mm', 'rms_mm', 'sampled_max_mm']:
        match = re.search('data-metric="' + key + '">([^<]+)', text)
        if not match or float(match.group(1)) != case['external_referent']['measured_original'][key]:
            raise ValueError('HTML_VALUE_NOT_BOUND: ' + key)
    if re.search('(?:src|href)=["\\\'](?:https?:|file:|/)', text):
        raise ValueError('HTML_EXTERNAL_ASSET')
    if '<script' in text.lower():
        raise ValueError('HTML_ACTIVE_SCRIPT')
    expected = (out / 'FROZEN_PREDICTIONS.json.sha256').read_text().strip()
    if sha(out / 'FROZEN_PREDICTIONS.json') != expected or expected != case['frozen_predictions_sha256']:
        raise ValueError('PREDICTION_FREEZE_HASH')
    return True
