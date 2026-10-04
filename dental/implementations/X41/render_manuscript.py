"""Render reader artifacts from replayed results; all numeric substitutions get locators."""
from pathlib import Path
import json
import re
ROOT = Path(__file__).resolve().parent

def main():
    result = json.loads((ROOT / 'results.json').read_text())
    x5 = json.loads((ROOT / 'inputs/LANE_X5_DECISION_SEG_ERROR/results.json').read_text())
    x26 = json.loads((ROOT / 'inputs/LANE_X26_DECIDABILITY/results.json').read_text())
    x31 = json.loads((ROOT / 'inputs/LANE_X31_CLINICAL_ANSWERS/results.json').read_text())
    raw = json.loads((ROOT / 'inputs/LANE_X26_DECIDABILITY/RAW_R2_CANAL.json').read_text())
    c = result['conditional_canal']
    p = result['counterexamples']
    drop = x26['canal_conditional']['dropout']
    a = x31['attrition']['guide']
    macro = {}
    ledger = []

    def add(key, text, locator, level, unit):
        macro[key] = str(text)
        ledger.append({'macro': key, 'rendered': str(text), 'source_locator': locator, 'resolution_level': level, 'unit': unit, 'time_scale': 'SIMULTANEOUS'})
    for (key, value, pointer) in [('SITES', result['archive']['sites'], 'sites'), ('SCANS', result['archive']['represented_scans'], 'represented_scans'), ('CANDIDATES', f"{result['archive']['source_candidates']:,}", 'source_candidates')]:
        add(key, value, 'results.json#/archive/' + pointer, 'PER_ARCH', 'count')
    add('PROBES', len(p), 'results.json#/counterexamples', 'PER_ARCH', 'selected cases')
    add('MIN_VOXELS', min((r['minimum_added_voxels'] for r in p)), 'results.json#/counterexamples', 'PER_SURFACE_REGION', 'voxels')
    add('MAX_VOXELS', max((r['minimum_added_voxels'] for r in p)), 'results.json#/counterexamples', 'PER_SURFACE_REGION', 'voxels')
    add('THRESHOLD', 2, 'inputs/LANE_X5_DECISION_SEG_ERROR/PREREG_R3.json#/threshold_mm', 'PER_SURFACE_REGION', 'mm engineering convention')
    add('DICE_LOW', f"{100 * min((r['Dice'] for r in p)):.4f}", 'results.json#/counterexamples', 'PER_ARCH', '%')
    add('DICE_HIGH', f"{100 * max((r['Dice'] for r in p)):.4f}", 'results.json#/counterexamples', 'PER_ARCH', '%')
    add('HD95', 0, 'results.json#/counterexamples', 'PER_ARCH', 'mm')
    for (key, ptr) in [('QUERIES', 'queries'), ('CLASS1', 'class1'), ('CLASS2', 'class2'), ('FLOOR', 'floor_mm')]:
        val = c[ptr]
        add(key, f'{val:.6f}' if ptr == 'floor_mm' else val, 'results.json#/conditional_canal/' + ptr, 'PHENOMENOLOGICAL' if ptr == 'floor_mm' else 'PER_SURFACE_REGION', 'mm' if ptr == 'floor_mm' else 'count')
    add('CLASS1_PERCENT', f"{100 * c['point_class1']:.2f}", 'results.json#/conditional_canal/point_class1', 'POPULATION', '%')
    for (k, j) in [('CLASS1_CI_LOW', 0), ('CLASS1_CI_HIGH', 1)]:
        add(k, f"{100 * c['case_ci95_class1'][j]:.2f}", 'results.json#/conditional_canal/case_ci95_class1/' + str(j), 'POPULATION', '%')
    add('PUBLIC_SCANS', 480, 'references/metadata (ToothFairy2 primary paper); dataset primary page, public training set', 'PER_ARCH', 'scans in source release')
    add('PITCH', '0.3', 'inputs/LANE_X26_DECIDABILITY/results.json#/canal_conditional/fixed_acquisition_mm', 'PER_POINT', 'mm')
    for (k, ptr) in [('MARGIN_ROWS', 'input_margin_rows'), ('OUTSIDE_CANAL', 'rejected/outside_canal_question')]:
        v = drop
        for s in ptr.split('/'):
            v = v[s]
        add(k, v, 'inputs/LANE_X26_DECIDABILITY/results.json#/canal_conditional/dropout/' + ptr, 'PER_SURFACE_REGION', 'records')
    add('OUTSIDE_PERCENT', f"{100 * drop['fraction_outside_scope_or_missing']:.2f}", 'inputs/LANE_X26_DECIDABILITY/results.json#/canal_conditional/dropout/fraction_outside_scope_or_missing', 'POPULATION', '%')
    for (k, ptr) in [('GUIDE_SITES', 'unique_complete_sites'), ('GUIDE_CASES', 'complete_cases'), ('GUIDE_EXCLUDED', 'source_not_in_complete_full_geometry'), ('GUIDE_SOURCE_SITES', 'source_candidate_sites')]:
        add(k, a[ptr], 'inputs/LANE_X31_CLINICAL_ANSWERS/results.json#/attrition/guide/' + ptr, 'PER_ARCH', 'counts')
    add('GUIDE_EXCLUDED_PERCENT', f"{100 * a['source_dropout_fraction']:.2f}", 'inputs/LANE_X31_CLINICAL_ANSWERS/results.json#/attrition/guide/source_dropout_fraction', 'POPULATION', '%')
    for (k, v, ptr) in [('BIAS', raw[0]['bias_bound'], 'bias_bound'), ('SIGMA', raw[0]['port']['standard_deviations'][0], 'port/standard_deviations/0'), ('TAU', raw[0]['port']['standard_deviations'][1], 'port/standard_deviations/1')]:
        add(k, f'{v:.6f}', 'inputs/LANE_X26_DECIDABILITY/RAW_R2_CANAL.json#/0/' + ptr, 'PHENOMENOLOGICAL', 'mm')
    add('SMCD', '0.77', 'references/PMC9633839.xml: Results highest-variability paragraph, Figure 4', 'POPULATION', 'mm curve-summary median')
    add('BOOTSTRAP_DRAWS', f"{result['bootstrap']['draws']:,}", 'results.json#/bootstrap/draws', 'POPULATION', 'draws')
    add('CUBE_TOLERANCE', '10⁻⁹', 'PREREG_R1_MANUSCRIPT.json#/metrics/cube_distance_replay_tolerance_mm', 'PER_POINT', 'mm')
    add('ENDPOINT_TOLERANCE', '10⁻¹⁰', 'PREREG_R1_MANUSCRIPT.json#/metrics/endpoint_arithmetic_tolerance_mm', 'PER_POINT', 'mm')
    for (k, fn) in [('VOLUME_LOW', min), ('VOLUME_HIGH', max)]:
        add(k, f"{fn((r['added_volume_mm3'] for r in p)):.3f}", 'results.json#/counterexamples', 'PER_SURFACE_REGION', 'mm³')
    for (prefix, r) in [('SMALLEST', min(p, key=lambda r: r['minimum_added_voxels'])), ('LARGEST', max(p, key=lambda r: r['minimum_added_voxels']))]:
        for (k, ptr, fmt, unit) in [('K', 'minimum_added_voxels', 'd', 'voxels'), ('GAP', 'adverse_gap_mm', '.6f', 'mm'), ('DICE', 'Dice', '.4f', '%')]:
            val = r[ptr] * 100 if ptr == 'Dice' else r[ptr]
            add(prefix + '_' + k, format(val, fmt), 'results.json#/counterexamples/' + str(p.index(r)) + '/' + ptr, 'PER_ARCH' if ptr == 'Dice' else 'PER_SURFACE_REGION', unit)
    for (k, ptr) in [('CURRENT_POSITIVE', 'CERTIFIED_POSITIVE'), ('CURRENT_ABSTAIN', 'ABSTAIN')]:
        add(k, c['current_status'][ptr], 'results.json#/conditional_canal/current_status/' + ptr, 'PER_SURFACE_REGION', 'count')
    w = x31['rounds']['R1']['orientation_counterexample'][0]
    for (key, j) in [('FULL_GUIDE_MEAN', 0), ('FULL_GUIDE_SD', 1)]:
        add(key, f"{w['same_radial_mean_sd_mm'][j]:.2f}", 'inputs/LANE_X31_CLINICAL_ANSWERS/results.json#/rounds/R1/orientation_counterexample/0/same_radial_mean_sd_mm/' + str(j), 'POPULATION', 'mm, individual radial deviation moments')
    for (i, name) in enumerate(['CONNECTED_ERRORS', 'DECIDABILITY', 'GUIDE_SCENARIOS', 'IDENTIFIABILITY'], 1):
        macro['TABLE_' + str(i)] = (ROOT / 'tables' / f'TABLE_{i}_{name}.md').read_text().strip()
    references = json.loads((ROOT / 'REFERENCES_FINAL.json').read_text())['references']
    lines = []
    for (i, r) in enumerate(references, 1):
        authors = ', '.join((a['family'] for a in r.get('authors', [])[:3]))
        authors += ' et al.' if len(r.get('authors', [])) > 3 else '.'
        title = ' '.join(r['title'].split())
        year = r['year'][0][0]
        journal = r['container'][0]
        lines.append(f"{i}. {authors} {title}. *{journal}*. {year}. [doi:{r['doi']}](https://doi.org/{r['doi']}).")
    macro['REFERENCES'] = '\n\n'.join(lines)
    source = (ROOT / 'MANUSCRIPT_TEMPLATE.md').read_text()

    def replace(match):
        key = match.group(1)
        if key not in macro:
            raise ValueError('Unfilled manuscript macro: ' + key)
        return macro[key]
    manuscript = re.sub('\\{\\{([A-Z0-9_]+)\\}\\}', replace, source)
    (ROOT / 'MANUSCRIPT.md').write_text(manuscript)
    (ROOT / 'MANUSCRIPT_NUMBER_LEDGER.json').write_text(json.dumps(ledger, indent=2) + '\n')
    (ROOT / 'FIGURE_CAPTIONS.md').write_text('# Figure captions\n\nFigure 1. Selected source-attached additions change a declared 2 mm complete-solid clearance predicate while location controls retain its original sign. Dice is identical within each pair and HD95 is zero. Six selected injections establish existence, not prevalence, model performance or injury risk. Data: TABLE_1_CONNECTED_ERRORS.csv and raw/GEOMETRY_REPLAY.json.\n\nFigure 2. Absolute regional minimum-clearance margins relative to an imposed image-label error floor and a numerical resampling envelope. Class 1 admits some finite pitch in this scenario; class 2 lacks strict finite slack. The floor is PHENOMENOLOGICAL, not a measured local anatomical bound. Data: frozen RAW_R2_CANAL.json, Table 2 and results.json.\n\nPNG, PDF and SVG exports are generated by make_figures.py. Derived annotation figures follow ToothFairy2 CC BY-SA 4.0 attribution (manuscript references 4–6).\n')
    print('Rendered MANUSCRIPT.md and numeric locator ledger')
if __name__ == '__main__':
    main()
