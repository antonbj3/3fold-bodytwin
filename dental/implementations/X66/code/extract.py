import re, json, xml.etree.ElementTree as ET
from lxml import etree
from common import *
SOURCES = {'L02': 'PMC10779861', 'L03': 'PMC10485502', 'L04': 'PMC10449641', 'L05': 'PMC10057599'}
TABLES = {'L02': 'materials-17-00263-t003', 'L03': 'jpbs_283_23-t001'}

def norm(s):
    return ' '.join(s.split())

def get_text(e):
    return norm(''.join(e.itertext()))

def load(lid):
    p = DENTAL / 'data/corpus/europepmc/fulltext' / (SOURCES[lid] + '.xml')
    return (p, etree.parse(str(p)))

def observation(lid, t, path, xpath, value, unit, quantity, **kw):
    els = t.xpath(xpath)
    assert len(els) == 1
    return dict(link_id=lid, source_id=SOURCES[lid], source_path=str(path.resolve()), source_sha256=sha(path), xpath=xpath, source_locator=str(path.resolve()) + '#' + xpath, original_text=get_text(els[0]), value=value, unit=unit, quantity=quantity, resolution_level='POPULATION', time_scale='SIMULTANEOUS', **kw)

def main():
    registry = []
    manifest = {}
    contracts = {}
    parsed = {}
    for (lid, pid) in SOURCES.items():
        (p, t) = load(lid)
        parsed[lid] = (p, t)
        lic = ' '.join(t.xpath('/article/front/article-meta/permissions//license//text()'))
        uri = re.search('https?://creativecommons.org/[^\\s)]+', lic)
        contracts[lid] = {'source_id': pid, 'source_path': str(p.resolve()), 'sha256': sha(p), 'doi': t.xpath('/article/front/article-meta/article-id[@pub-id-type="doi"]/text()')[0], 'title': norm(''.join(t.xpath('/article/front/article-meta/title-group/article-title//text()'))), 'license_text': norm(lic), 'license_uri': uri.group(0) if uri else 'UNKNOWN', 'lineage': 'One study; repeated/group summaries are not independent physical replicates.', 'external_referent': {'kind': 'independent_measurement', 'locator': str(p.resolve()), 'compared_quantity': {'L02': 'Aged cFDP fracture load and aging failure count', 'L03': 'Post-cementation site mean vertical marginal gap', 'L04': 'Normalized ISO fatigue survival and printed endpoint relationship', 'L05': 'Transmitted optical dose and no-luting intrapulpal temperature rise'}[lid], 'refutes_us': True}}
        manifest[str(p.resolve())] = {'sha256': sha(p), 'bytes': p.stat().st_size}
    second = {}
    for lid in ['L02', 'L03']:
        (p, t) = parsed[lid]
        tid = TABLES[lid]
        root = ET.parse(p).getroot()
        tab = next((e for e in root.iter('table-wrap') if e.get('id') == tid))
        rows2 = [[get_text(c) for c in tr] for tr in tab.iter('tr')]
        rows = [[get_text(c) for c in tr] for tr in t.xpath('//table-wrap[@id="' + tid + '"]//tr')]
        assert rows == rows2
        second[lid] = {'parser': 'xml.etree.ElementTree', 'exact_cells_match': True, 'rows': len(rows), 'cell_count': sum(map(len, rows))}
    (p, t) = parsed['L02']
    rows = t.xpath('//table-wrap[@id="' + TABLES['L02'] + '"]//tr')[2:]
    crown = []
    for (i, tr) in enumerate(rows):
        cells = [get_text(c) for c in tr]
        grp = cells[0]
        fail = 0 if cells[1] == '-' else int(cells[1])
        mean = float(cells[2])
        sd = float(cells[3])
        for (col, val, q, u) in [(2, fail, 'observed_aging_failures', 'count'), (3, mean, 'survivor_fracture_load_mean', 'N'), (4, sd, 'survivor_fracture_load_SD', 'N')]:
            registry.append(observation('L02', t, p, f'''(//table-wrap[@id="{TABLES['L02']}"]//tr)[{i + 3}]/*[{col}]''', val, u, q, group=grp, initial_n=8, statistic='reported' if q.endswith('failures') else 'reported_mean' if 'mean' in q else 'reported_SD'))
        crown.append({'group': grp, 'early_failures': fail, 'mean_N': mean, 'SD_N': sd, 'initial_n': 8, 'survivor_n': 8 - fail, 'failure_horizon_cycles': 1200000, 'wall_mm': float(grp.rsplit('-', 1)[-1]), 'quantity': 'fracture load of surviving cantilever FDPs, not flexural strength', 'source_locator': registry[-2]['source_locator']})
    contracts['L02'].update({'scope': 'Five n8 cantilever FDP groups, resilient CoCr dies, >=12 mm2 connector, 10,000 thermocycles6.5/60C and1.2e6 cycles108N, survivors destructive fracture; repeated mean does not include aged failures.', 'time_scale': 'HANDOVER', 'unit_locator': '//table-wrap[@id="materials-17-00263-t003"]//th[contains(.,"[N]")]', 'source_conflict': 'Table3 3Y mean571 N vs Results narrative577 N. Keep table571 and flag conflict; do not average.', 'methods_locator': '/article/body/sec[2]'})
    (p, t) = parsed['L03']
    rows = t.xpath('//table-wrap[@id="' + TABLES['L03'] + '"]//tr')[1:]
    fit = []
    for (i, tr) in enumerate(rows):
        c = [get_text(x) for x in tr]
        grp = 'ABC'[i // 4]
        site = c[1]
        n = int(c[2])
        mean = float(c[3])
        sd = float(c[4])
        assert n == 8
        for (col, v, q) in [(4, mean, 'marginal_gap_mean'), (5, sd, 'marginal_gap_SD')]:
            registry.append(observation('L03', t, p, f'''(//table-wrap[@id="{TABLES['L03']}"]//tr)[{i + 2}]/*[{col}]''', v, 'um', q, group=grp, surface_region=site, n_copings=8, statistic='reported_mean' if q.endswith('mean') else 'reported_SD', finest_common_resolution='PER_SURFACE_REGION'))
        fit.append({'group': grp, 'site': site, 'n_copings': 8, 'mean_um': mean, 'SD_um': sd, 'source_locator': registry[-2]['source_locator'], 'resolution_level': 'PER_SURFACE_REGION'})
    contracts['L03'].update({'scope': '24 angled15deg Ti abutments;8 copings/group,0.6mm Zr;A0.75mm vent,Bopen vent,Cprecementation;ZOnE cement,20N60s;stereomicroscope30x;four repeated sites D/M/B/L per coping.', 'time_scale': 'HANDOVER', 'unit_locator': '/article/body/sec[3]/fig/caption/p', 'methods_locator': '/article/body/sec[2]', 'missing_inputs': ['CAD spacer', 'pre-cementation paired gaps', 'per-coping site covariance', 'source mesh'], 'source_conflict': 'Article inferential labels compare within-group sites and call inter-group; retained p-values are not reinterpreted.'})
    (p, t) = parsed['L04']
    xp = '/article/body/sec[3]/sec[2]/sec[2]/p'
    text = get_text(t.xpath(xp)[0])
    a = float(re.search('a = (\\d+(?:\\.\\d+)?)', text).group(1))
    b = -float(re.search('b = −(\\d+(?:\\.\\d+)?)', text).group(1))
    frac = float(re.search('equal to ([0-9.]+)%', text).group(1)) / 100
    ci = [float(v) / 100 for v in re.search('\\[([0-9.]+)%, ([0-9.]+)%\\]', text).groups()]
    iso = {'a': a, 'b': b, 'reported_fraction': frac, 'reported_fraction_CI': ci, 'cycles': 5000000, 'n_at_load_per_design': 4, 'designs': ['50SP', '60SP', '70SP', '50SW', '60SW', '70SW'], 'observed_all_survive_fraction': 0.1, 'observed_at_least_one_failure_fraction': 0.15, 'load_rounding_N': 10, 'mean_absolute_ultimate_load_N': None}
    for k in ['a', 'b', 'reported_fraction']:
        registry.append(observation('L04', t, p, xp, iso[k], '1', k, statistic='reported_fit_coefficient' if k in ['a', 'b'] else 'reported_fit_endpoint'))
    for (i, v) in enumerate(ci):
        registry.append(observation('L04', t, p, xp, v, '1', 'reported_fit_CI_' + ['low', 'high'][i], statistic='reported_95pct_CI'))
    registry.append(observation('L04', t, p, xp, 0.1, '1', 'nominal_fraction_all_four_runout_per_design', n_per_design=4, statistic='reported_nominal_load_level'))
    registry.append(observation('L04', t, p, xp, 0.15, '1', 'nominal_fraction_at_least_one_failure_per_design', n_per_design=4, statistic='reported_nominal_load_level'))
    contracts['L04'].update({'scope': '168 Ti64 SLM samples,4.5mm cylinder and4.5mm porous region;SP/SW50/60/70%;11mm arm30deg R0.1 15Hz5e6cycles;load levels rounded10N;4 specimens per design/load.', 'methods_locator': '/article/body/sec[2]/sec[3]/sec[2]', 'equation_locator': '/article/body/sec[2]/sec[3]/sec[2]/p[2]/disp-formula', 'missing_inputs': ['Digitized individual fatigue/monotonic forces', 'unrounded coefficients', 'same-design measured absolute Fult', 'target as-built mesh'], 'unit_locator': xp})
    (p, t) = parsed['L05']
    xp = '/article/body/sec[2]/p[1]'
    s = get_text(t.xpath(xp)[0])
    xpd = '/article/body/sec[3]/p[6]'
    disc = get_text(t.xpath(xpd)[0])
    patterns = ['20\\.5 ± 0\\.4', '18\\.6 J/cm2 ± 0\\.3', '16\\.6 ± 0\\.3', '14\\.9 ± 0\\.3']
    dose = []
    for (thick, pat) in zip([2.0, 2.5, 3.0, 3.5], patterns):
        match = re.search(pat, s)
        assert match
        nums = [float(x) for x in re.findall('(?<![A-Za-z])\\d+(?:\\.\\d+)?', match.group())]
        val = nums[0]
        sd = nums[-1]
        registry.append(observation('L05', t, p, xp, val, 'J/cm2', 'transmitted_radiant_exposure_mean', ceramic_thickness_mm=thick, reported_uncertainty=sd, uncertainty_kind='as_printed_not_assumed_SD', statistic='reported_mean'))
        dose.append({'ceramic_mm': thick, 'E_J_cm2': val, 'uncertainty_J_cm2': sd, 'duration_s': 40.0, 'source_locator': registry[-1]['source_locator']})
    delta = float(re.search('intrapulpal temperature rise was ([0-9.]+) °C without', disc).group(1))
    assert delta == 4.1
    registry.append(observation('L05', t, p, xpd, delta, 'degC', 'D1_C3p5_no_luting_pulp_temperature_rise', statistic='reported_mean', dentin_mm=1.0, ceramic_mm=3.5))
    pulp = {'doses': dose, 'incident_E_J_cm2': 62.0, 'incident_uncertainty_J_cm2': 0.6, 'duration_s': 40.0, 'measured_no_luting_deltaT_C': delta, 'dentin_mm': 1.0, 'ceramic_mm': 3.5, 'exo_heat_W_m2': None}
    contracts['L05'].update({'scope': 'In vitro representative third molar, varied dentin/ceramic thickness,40s light curing;4.1 C no-luting D1C3.5;no perfusion, repeated tests on same tooth;optical transmission does not measure absorption.', 'methods_locator': '/article/body/sec[4]', 'unit_locator': xp, 'time_scale': 'SIMULTANEOUS', 'exotherm_missing': 'Full cementation exceeds5.5C in text, numerical per-combination values in figures not digitized here.'})
    support = DENTAL / 'notes/occlusion_support_intervals/SUPPORT_INTERVALS.json'
    data = json.loads(support.read_text())
    pdl = next((r for r in data if r['id'] == 'NAT-INC-TIME'))
    manifest[str(support)] = {'sha256': sha(support), 'bytes': support.stat().st_size}
    contracts['L01'] = {'source_id': 'SUPPORT-INTERVALS', 'source_path': str(support), 'sha256': sha(support), 'primary_exists': False, 'required_original': str(DENTAL / 'data/corpus/europepmc/fulltext/PMC9889448.xml'), 'scope': 'Local curator note only, mixed text/table alternatives; no independent primary admission.', 'external_referent': {'kind': 'our_own_fixture', 'locator': str(support) + '#NAT-INC-TIME', 'compared_quantity': 'Curated conflicting healthy10s forces', 'refutes_us': True}, 'license_uri': 'Local research note; no primary redistribution'}
    for rec in registry:
        if rec['link_id'] in ['L02', 'L03']:
            rec['time_scale'] = 'HANDOVER'
        if rec['link_id'] == 'L03':
            rec['resolution_level'] = 'PER_SURFACE_REGION'
        if rec['link_id'] == 'L05':
            rec['resolution_level'] = 'PER_TOOTH' if rec['quantity'] == 'D1_C3p5_no_luting_pulp_temperature_rise' else 'PER_SURFACE_REGION'
    contracts['L03']['license_uri'] = 'https://creativecommons.org/licenses/by-nc-sa/4.0/'
    contracts['L03']['license_basis'] = 'Explicit original permissions text names Attribution-NonCommercial-ShareAlike4.0; URI reconstructed from that name, not a downloaded license claim.'
    expected = {'crown_means': [571.0, 330.0, 363.0, 290.0, 291.0], 'crown_SD': [85.0, 82.0, 134.0, 128.0, 82.0], 'crown_failures': [0, 0, 1, 2, 1], 'fit_means': [158.688, 150.8, 155.138, 169.175, 158.575, 150.525, 148.188, 150.075, 150.288, 155.6, 187.15, 161.288], 'doses': [20.5, 18.6, 16.6, 14.9], 'a': 3.22, 'b': -0.2, 'reported_fraction': 0.1286, 'deltaT': 4.1}
    checks = {'crown_means': [r['mean_N'] for r in crown] == expected['crown_means'], 'crown_SD': [r['SD_N'] for r in crown] == expected['crown_SD'], 'crown_failures': [r['early_failures'] for r in crown] == expected['crown_failures'], 'fit_means': [r['mean_um'] for r in fit] == expected['fit_means'], 'doses': [r['E_J_cm2'] for r in dose] == expected['doses'], 'a': a == expected['a'], 'b': b == expected['b'], 'reported_fraction': abs(frac - expected['reported_fraction']) < 1e-16, 'deltaT': delta == expected['deltaT']}
    assert all(checks.values())
    put('raw/observations.json', {'L01': pdl, 'L02': crown, 'L03': fit, 'L04': iso, 'L05': pulp})
    put('OBSERVATION_REGISTRY.json', registry)
    put('SOURCE_CONTRACTS.json', contracts)
    put('SOURCE_MANIFEST.json', manifest)
    put('FACIT.json', {'source_witnesses': expected, 'checks': checks, 'all_pass': True, 'second_parser': second, 'notice': 'Manually read printed source values, not held-out experiment predictions.'})
    print('Extracted', len(registry), 'quantity records from4 primary originals; PDL primary missing.')
if __name__ == '__main__':
    main()
