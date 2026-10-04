"""Curated primary-source observation selectors. No fitted physical parameters."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, re, xml.etree.ElementTree as ET, hashlib
P = Path(__file__).resolve().parent
S = Path(__import__('os').environ.get('XFACIT_SOURCE_DIR', _release_expand('@DENTAL_WORK_ROOT@/XFACIT-harvest/sources')))

def clean(el):
    if el is None:
        return ''
    return re.sub('\\s+', ' ', ''.join(el.itertext())).strip()

def num(s):
    s = s.replace('−', '-').replace('–', '-')
    return re.findall('(?<![\\w.])[-+]?\\d+(?:\\.\\d+)?', s)

def metadata(file):
    p = S / file
    if file.endswith('.xml'):
        root = ET.parse(p).getroot()
        doi = next((clean(x) for x in root.findall('.//article-id') if x.get('pub-id-type') == 'doi'), None)
        title = clean(root.find('.//article-title'))
        licenses = [clean(x) for x in root.findall('.//license')]
    elif file.endswith('.json'):
        d = json.loads(p.read_text())
        r = d['resultList']['result'][0]
        doi = r.get('doi')
        title = r.get('title')
        licenses = ['Abstract rights not granted: private source verification only']
    else:
        doi = {'frontiers.html': '10.3389/fmech.2026.1788600', 'ferrato.html': '10.11138/ads/2017.8.2.089', 'hattori.pdf': '10.7144/sgf.2.111'}[file]
        title = {'frontiers.html': 'Experimental characterisation of coefficient of friction at bone-implant interface', 'ferrato.html': 'Digital evaluation of occlusal forces: comparison between healthy subjects and TMD patients', 'hattori.pdf': 'Bite force distribution on dental arch during clenching'}[file]
        licenses = ['CC BY: publisher statement'] if file == 'frontiers.html' else ['Reuse license UNKNOWN; private verification only']
    return {'source_file': str(p), 'source_sha256': hashlib.sha256(p.read_bytes()).hexdigest(), 'doi': doi, 'url': 'https://doi.org/' + doi if doi else None, 'title': title, 'license': licenses}
OBS = []

def table(id, file, table_id, row, col, value, unit, quantity, population, protocol, region, edges, variable=None, sd=None, sd_col=None, kind='independent_measurement', stat='reported_group_mean', context=None):
    root = ET.parse(S / file).getroot()
    t = root.find('.//table-wrap[@id="' + table_id + '"]')
    rr = t.findall('.//tr')[row]
    cs = list(rr)
    cell = clean(cs[col])
    printed = num(cell)
    value = str(value)
    sd = str(sd) if sd is not None else None
    if value not in printed:
        raise ValueError((id, 'value not in selected cell', value, cell))
    sdcell = cell if sd_col is None else clean(cs[sd_col])
    if sd is not None and sd not in num(sdcell):
        raise ValueError((id, 'sd not in selected cell', sd, sdcell))
    OBS.append({'id': id, 'claim_type': 'information_link', 'quantity': quantity, 'value': value, 'unit': unit, 'dispersion': {'kind': 'reported_SD', 'value': sd} if sd else {'kind': 'UNKNOWN', 'value': None}, 'statistic': stat, 'resolution_level': 'POPULATION', 'measurement_support_resolution': region, 'population': population, 'protocol': protocol, 'context': context or {}, 'source_kind': kind, 'reference': {**metadata(file), 'locator': f"{clean(t.find('label'))} [{table_id}], raw tr[{row}], cell[{col}]; {cell}", 'selector': {'kind': 'xml_table_cell', 'table_id': table_id, 'row': row, 'col': col}, 'source_cell': cell, 'sd_selector_col': sd_col}, 'edges': edges, 'target_variable': variable, 'transfer_status': 'UNKNOWN_TARGET_PROTOCOL', 'measurement_truth_note': 'Printed observation verified; apparatus error and target transfer not identified.'})

def abstract(id, file, anchor, value, unit, quantity, population, protocol, region, edges, variable=None, sd=None, dispersion='UNKNOWN', stat='reported_group_mean'):
    d = json.loads((S / file).read_text())
    recs = d['resultList']['result']
    r = next((x for x in recs if anchor in (x.get('abstractText') or '')))
    if str(value) not in anchor or (sd is not None and str(sd) not in anchor):
        raise ValueError(id)
    OBS.append({'id': id, 'claim_type': 'information_link', 'quantity': quantity, 'value': str(value), 'unit': unit, 'dispersion': {'kind': dispersion, 'value': str(sd) if sd is not None else None}, 'statistic': stat, 'resolution_level': 'POPULATION', 'measurement_support_resolution': region, 'population': population, 'protocol': protocol, 'context': {}, 'source_kind': 'independent_measurement', 'reference': {**metadata(file), 'doi': r.get('doi'), 'url': 'https://doi.org/' + r['doi'] if r.get('doi') else 'https://pubmed.ncbi.nlm.nih.gov/' + r['id'] + '/', 'locator': 'Primary abstract, exact numerical phrase; PMID ' + r['id'], 'selector': {'kind': 'abstract_phrase', 'record_id': r['id'], 'anchor': anchor}, 'source_cell': anchor}, 'edges': edges, 'target_variable': variable, 'transfer_status': 'UNKNOWN_TARGET_PROTOCOL', 'measurement_truth_note': 'Abstract numerical results only; full table and apparatus uncertainty not available.'})
for (k, base, values, sds, q, u, var, edges) in [('IT', 1, ['34.50', '46.50', '43.50', '46.50'], ['3.59', '3.28', '2.96', '3.28'], 'insertion_torque', 'Ncm', 'insertion_torque', ['D-E-K26']), ('ISQ', 5, ['59.85', '75.85', '71.90', '75.85'], ['1.42', '2.32', '2.40', '2.32'], 'isq', '1', 'isq', ['D-E-K31']), ('RT', 9, ['35.50', '47.25', '43.00', '47.25'], ['3.94', '3.02', '2.51', '3.02'], 'removal_breakaway_torque', 'Ncm', None, ['D-E-K26'])]:
    for (j, (v, sd)) in enumerate(zip(values, sds)):
        g = ['1a', '2a', '1b', '2b'][j]
        table(f'GANDHI_{k}_{g}', 'PMC9930143.xml', 'tbl1', base + j, 3 if j == 0 else 2, v, u, q, '80 implants in ex vivo pig tibia; n=20/group', 'Two implant designs; conventional versus osseodensification; immediate insertion/ISQ then removal, breakaway motor setting', 'PER_TOOTH', edges, var, sd=sd, sd_col=4 if j == 0 else 3, context={'group': g})
from bs4 import BeautifulSoup
soup = BeautifulSoup((S / 'frontiers.html').read_text(), 'html.parser')
ft = next((t for t in soup.find_all('table') if 'Sawbones' in t.get_text()))
for (j, (label, v, sd, n)) in enumerate([('Sawbones_Smooth', '0.27', '0.03', 6), ('Sawbones_Porous', '0.99', '0.09', 6), ('Porcine_Smooth', '0.54', '0.06', 4), ('Porcine_Porous', '0.93', '0.18', 7), ('Bovine1_Smooth', '0.46', '0.02', 2), ('Bovine1_Porous', '1.00', '0.03', 3), ('Bovine2_Porous', '0.68', '0.09', 3), ('Bovine3_Porous', '0.20', '0.01', 3), ('Cadaveric_Smooth', '0.41', '0.08', 3), ('Cadaveric_Porous', '0.99', '0.14', 5)]):
    cells = [[c.get_text(' ', strip=True) for c in row.find_all(['td', 'th'])] for row in ft.find_all('tr')]
    row = next((i for (i, cs) in enumerate(cells) if len(cs) >= 5 and cs[-2:] == [v, sd] and (cs[-3] == str(n))))
    cell = cells[row][-2]
    OBS.append({'id': 'MIN_' + label, 'claim_type': 'information_link', 'quantity': 'friction', 'value': v, 'unit': '1', 'dispersion': {'kind': 'reported_SD', 'value': sd}, 'statistic': 'reported_group_mean', 'resolution_level': 'POPULATION', 'measurement_support_resolution': 'PER_SURFACE_REGION', 'population': f'{label}; n={n}; total42 tests, five bone/surface contexts', 'protocol': 'OUKR tibial keel; smooth Cp-Ti versus porous; fixed horizontal 2kg mass, variable normal load to initiation of slip; Sawbones dry', 'context': {'bone_surface': label}, 'source_kind': 'independent_measurement', 'reference': {**metadata('frontiers.html'), 'locator': f'Table1 HTML row {row}, mean/SD; {label}', 'selector': {'kind': 'html_table_cell', 'contains': 'Sawbones', 'row': row, 'col': len(cells[row]) - 2}, 'source_cell': cell}, 'edges': ['D-E-K19', 'D-E-K26', 'D-E-K27'], 'target_variable': 'friction', 'transfer_status': 'UNKNOWN_TARGET_PROTOCOL', 'measurement_truth_note': 'Knee keel on tibial bone, not dental screw traction. No universal dental friction band.'})
for (mat, v, sd) in [('Katana HT', '1194', '111'), ('Katana UTML', '688', '159'), ('e.max LT', '450', '53')]:
    anchor = v + ' ±' + sd + ' (' + mat + ')'
    abstract('KWON_' + mat.replace(' ', '_'), 'PMID29310875.json', anchor, v, 'MPa', 'flexural_strength_mean', 'n=10/material, 25x4x2mm polished bars', 'Sintered/crystallized bars, fracture at1mm/min', 'PER_TOOTH', ['D-E-K32', 'D-E-K36'], sd=sd, dispersion='reported_plus_minus_unspecified_in_abstract')
for (j, (m, v)) in enumerate([('SA', '58.604'), ('DP', '48.759')], 2):
    table('ZHAO_' + m, 'PMC10067253.xml', 'Tab1', j, 1, v, 'mm2', 'arch_occlusal_contact_area', '22 participants; repeated measures', 'Scanned articulating marks vs Dental Prescale; matched same participant occlusal analyses', 'PER_ARCH', ['D-E-K16', 'D-E-K21'], context={'method': m})
for (direction, v, sd) in [('longitudinal', '0.58', '0.018'), ('circumferential', '0.53', '0.030'), ('radial', '0.54', '0.020')]:
    abstract('DAVIDSON_' + direction, 'PMID11334760.json', v + '+/-' + sd + ' W/mK', v, 'W/(m K)', 'thermal_conductivity', 'Bovine femora mid-diaphysis cortical specimens, n UNKNOWN in abstract', 'Three orthogonal heat-transfer directions: ' + direction, 'PER_SURFACE_REGION', ['D-E-K24', 'D-E-K23'], sd=sd, dispersion='reported_plus_minus_unspecified_in_abstract')
for (j, (name, v, sd)) in enumerate([('none', '6.8', '0.8'), ('water10', '-16.3', '2.1'), ('water23', '-8.6', '1.7'), ('water35', '1.6', '1.2')], 5):
    table('FARAH_' + name, 'PMC6387888.xml', 'T1', j, 2 if j == 5 else 1, v, 'K', 'pulp_temperature_change', '28 molars;7/group', '60s continuous occlusal preparation,50ml/min water;K-type thermocouple2s sampling;baseline-to-extreme', 'PER_POINT', ['D-E-K23', 'D-E-PULP-HEAT', 'D-E-K49'], sd=sd, context={'coolant': name})
for (j, (condition, wear, vals, sds)) in enumerate([('pulp', 'new', ['27.45', '21.57', '20.71'], ['0.66', '0.27', '0.20']), ('pulp', 'fifth', ['29.15', '21.96', '20.74'], ['0.67', '0.19', '0.31']), ('pulp', 'eighth', ['28.67', '21.28', '20.35'], ['1.18', '0.37', '0.21']), ('paste', 'new', ['28.99', '21.14', '20.29'], ['1.50', '0.31', '0.25']), ('paste', 'fifth', ['29.77', '21.25', '20.25'], ['0.47', '0.17', '0.17']), ('paste', 'eighth', ['29.76', '21.53', '20.21'], ['0.41', '0.40', '0.17'])], 2):
    for (k, (v, sd)) in enumerate(zip(vals, sds)):
        table(f'VENEER_{condition}_{wear}_{k}', 'PMC10453045.xml', 'dentistry-11-00197-t002', j, 3 + 2 * k if j in [2, 5] else 2 + 2 * k, v, 'degC', 'pulp_temperature', '30 teeth,3 bur groups;n=5 per pulp-condition/usage group,90 repeated time observations', 'Water-cooled veneer preparation;temperature initially/1min/3min;retained pulp or thermoconductive paste', 'PER_POINT', ['D-E-K23', 'D-E-PULP-HEAT', 'D-E-K49'], 'pulp_temperature', sd=sd, context={'condition': condition, 'wear': wear, 'minute': [0, 1, 3][k]})
for (j, (group, v, sd)) in enumerate([('D0.5', '1393', '301'), ('D1.0', '1789', '233'), ('GIC0.5', '428', '93'), ('GIC1.0', '541', '173'), ('RC0.5', '1043', '339'), ('RC1.0', '1187', '317')], 1):
    table('SUBSTRATE_' + group, 'PMC10756807.xml', 'TB2272278-1', j, 3 if j in [1, 3, 5] else 2, v, 'N', 'fracture_mean', 'Disc/substrate restoration groups;Table1 means;independent cell n UNKNOWN (Methods40 ceramic discs vs6 groups)', '10mm lithium-disilicate discs,0.5/1mm;G10,dentin-analogue with GIC/RC core;1e6cycles100N4Hz37C+1e4thermocycles;load-to-failure0.1mm/min', 'PER_TOOTH', ['D-E-K36', 'D-E-K49'], 'fracture_mean', sd=sd, context={'group': group, 'thickness_mm': group[-3:], 'foundation': group[:-3]})
for (j, (code, it, it_sd, rt, rt_sd, isq, isq_sd)) in enumerate([('12S', '13.67', '1.88', '11.68', '1.06', '51.40', '2.95'), ('06D', '20.19', '1.61', '17.45', '1.28', '53.77', '2.73'), ('06S', '22.30', '1.78', '20.25', '2.47', '55.66', '1.62')], 1):
    for (q, col, v, sd, u, var, edges) in [('IT', 1, it, it_sd, 'Ncm', 'insertion_torque', ['D-E-K26']), ('RT', 2, rt, rt_sd, 'Ncm', None, ['D-E-K26']), ('ISQ', 3, isq, isq_sd, '1', 'isq', ['D-E-K31'])]:
        table('YAMAGUCHI_' + code + '_' + q, 'PMC7438427.xml', 'Tab2', j, col, v, u, {'IT': 'insertion_torque', 'RT': 'removal_torque', 'ISQ': 'isq'}[q], 'Artificial bone implant design groups;sample size seeMethods', 'Matched custom single/double thread implants;immediate insertion,removal and RFA', 'PER_TOOTH', edges, var, sd=sd, context={'design': code})
for (j, (bone, v, sd)) in enumerate([('I', '36.50', '3.37'), ('III', '31.20', '5.92')], 3):
    root = ET.parse(S / 'PMC8630976.xml').getroot()
    t = root.find('.//table-wrap[@id="TB_4"]')
    rows = t.findall('.//tr')
    row = next((i for (i, r) in enumerate(rows) if clean(list(r)[0]) == bone and '±' in clean(r)))
    if bone == 'III':
        x = list(rows[row])[1]
        x.text = (x.text or '').replace('31,20', '31.20')
        orig = (S / 'PMC8630976.xml').read_text()
        assert '31,20' in orig
        expected = '31'
    else:
        expected = v
    table('CLINICAL_IT_' + bone, 'PMC8630976.xml', 'TB_4', row, 1, expected, 'Ncm', 'insertion_torque', '25 clinical implants in9people;typeI n10,typeIII n15', '11.5x3.75mm DSP Biomedical implant;maximum insertion motor torque', 'PER_TOOTH', ['D-E-K26'], 'insertion_torque', sd=sd, context={'bone_type': bone})
    if bone == 'III':
        OBS[-1]['value'] = '31.20'
        OBS[-1]['reference']['numeric_conversion'] = 'decimal comma in cell:31,20->31.20'
vals = ['38.08', '43.18', '45.28', '38.81', '45.97', '50.06', '50.31', '54.39', '66.89', '54.00', '59.94', '66.78']
sds = ['1.01', '2.92', '3.56', '1.93', '3.93', '4.38', '3.45', '4.69', '26.57', '4.46', '6.32', '15.33']
for (j, (v, sd)) in enumerate(zip(vals, sds), 1):
    table('DRILL_' + str(j), 'PMC12423576.xml', 'clr14462-tbl-0001', j, 5, v, 'degC', 'drill_temperature_peak', '120osteotomies in cortical bovine bone;12groups n10', 'Direct thermography at exposed dynamic bone-drill interface;2.4pilot->3.2/3.8cutting diameter;3.5mmcortex;200/600/1000rpm;37C baseline', 'PER_POINT', ['D-E-K24', 'D-E-K28', 'D-E-K49'], None, sd=sd, sd_col=6, context={'irrigation': 'on' if j <= 6 else 'off', 'rpm': [200, 600, 1000][(j - 1) % 3], 'expansion_mm': '0.8' if (j - 1) % 6 < 3 else '1.4'})
for (j, (mat, vals, sds)) in enumerate([('ESS', ['14.65', '48.65', '90.92', '154.02', '246.45'], ['2.25', '6.17', '11.32', '13.50', '31.35']), ('EMX', ['41.29', '124.95', '209.29', '355.65', '709.00'], ['7.03', '15.60', '36.39', '100.76', '87.37']), ('LP', ['104.97', '348.76', '768.17', '1244.02', '1844.81'], ['19.71', '68.40', '140.31', '129.01', '216.77'])], 3):
    for (k, (v, sd)) in enumerate(zip(vals, sds), 1):
        table(f'BIAXIAL_{mat}_{k}', 'PMC10004144.xml', 'materials-16-01997-t002', j, k, v, 'N', 'fracture_mean', '12discs/material-thickness;180total', 'ISO6872biaxial bending,polished discs,3ball10mm support diameter,1.4mmpiston,0.5mm/min;not supported crown', 'PER_TOOTH', ['D-E-K32', 'D-E-K36', 'D-E-K49'], 'fracture_mean', sd=sd, context={'material': mat, 'thickness_mm': ['0.4', '0.7', '1.0', '1.3', '1.6'][k - 1]})
for (k, (mat, vals, sds)) in enumerate([('Zirconia', ['14.24', '27.51', '20.68', '23.41', '21.45'], ['9.7', '12.81', '10.24', '13.85', '12.58']), ('Lithium_disilicate', ['111.64', '23.96', '77.3', '36.22', '62.28'], ['48.94', '9.43', '55.81', '21.33', '51.8']), ('Composite', ['69.72', '35.61', '47.77', '25.71', '44.7'], ['25.8', '8.92', '23.67', '11.79', '24.96'])], 1):
    for (j, (v, sd)) in enumerate(zip(vals, sds), 1):
        table(f'SEM_{mat}_{j}', 'PMC10058296.xml', 'materials-16-02413-t002', j, k, v, 'um', 'marginal_gap', '9crowns/material;27points/region,108points/material;not108independent teeth', 'CAD/CAMcrowns onprinted resin maxillary firstmolar preparation,30umCADspace,SEMaftercementation,deepchamfer;vestibular/palatal/mesial/distal/total', 'PER_SURFACE_REGION', ['D-E-K43', 'D-E-K44', 'D-E-K48'], 'marginal_gap', sd=sd, context={'material': mat, 'region': ['vestibular', 'palatal', 'mesial', 'distal', 'total'][j - 1]})
for (id, anchor, v, sd) in [('THIN', '1.61 +/- 0.24 mm', '1.61', '0.24'), ('THICK', '0.26 +/- 0.08 mm', '0.26', '0.08')]:
    abstract('MBL_' + id, 'PMID19885413.json', anchor, v, 'mm', 'mbl', '46implants(23test/23control)in19patients', 'Prospective1yr;test2mmsupracrestal;initialmucosalthickness withprobe;mesialradiographicbonechange;' + id, 'PER_SURFACE_REGION', ['D-E-K30'], 'mbl', sd=sd, dispersion='reported_SE')
for (t, triplet) in [(4, ['21.2', '72.1', '33.8']), (8, ['21.8', '70.9', '30.0']), (16, ['16.3', '70.2', '29.9'])]:
    anchor = ', '.join((x + '%' for x in triplet[:2])) + ', and ' + triplet[2] + '% at ' + str(t) + ' weeks'
    for (treatment, v) in zip(['untreated', 'HA', 'alkali_heat'], triplet):
        abstract(f'KAWAI_BIC_{t}_{treatment}', 'PMID25242248.json', anchor, v, 'percent', 'bic', 'Rabbit model;n8perreportedtimepoint', 'Rougharc-sprayedTi;' + treatment + 'surface;histologicalbonecontact at' + str(t) + 'weeks', 'PER_SURFACE_REGION', ['D-E-K19', 'D-E-K31'], 'bic')
        OBS[-1]['context'] = {'week': t, 'treatment': treatment}
for (mat, v, sd) in [('I1', '15.0', '25.7'), ('I2', '12.6', '22.1'), ('C', '19.7', '28.6'), ('P1', '23.1', '25.0'), ('P2', '30.5', '31.8'), ('M1', '130.4', '67.0'), ('M2', '211.8', '111.3')]:
    OBS.append({'id': 'HATTORI_' + mat, 'claim_type': 'information_link', 'quantity': 'tooth_force', 'value': v, 'unit': 'N', 'dispersion': {'kind': 'reported_SD', 'value': sd}, 'statistic': 'reported_group_mean', 'resolution_level': 'POPULATION', 'measurement_support_resolution': 'PER_TOOTH', 'population': '42normaldentatesubjects;65mandibularsideswithoutM3contacts;pairedsidesnotindependent', 'protocol': 'Maximalvoluntaryclenchintercuspal;DentalPrescale50HtypeR;MANDIBULARarch(MethodsJapanesetext下顎)', 'context': {'tooth_type': mat, 'arch': 'mandibular'}, 'source_kind': 'independent_measurement', 'reference': {**metadata('hattori.pdf'), 'locator': 'PDFprintedp113Table1 wo/M3mean/SD;visuallycheckedrenderhattori_p113.png', 'selector': {'kind': 'visual_pdf_table', 'page': 113, 'row': mat}, 'source_cell': v + ' (' + sd + ')'}, 'edges': ['D-E-K16', 'D-E-K15'], 'target_variable': 'tooth_force', 'transfer_status': 'UNKNOWN_TARGET_PROTOCOL', 'measurement_truth_note': 'TableNvaluesverifiedvisually;plotsharesnotre-digitized.Maxillaryconsumerassignmentrejected.'})
soup = BeautifulSoup((S / 'ferrato.html').read_text(), 'html.parser')
sec = soup.find(id='t3-79-88')
tt = sec.find_all('table')[0]
rows = tt.find_all('tr')
cells = [[c.get_text(' ', strip=True) for c in r.find_all(['th', 'td'])] for r in rows]
means = ['14.5', '14.7', '7.6', '6.8', '2.2', '1.0', '1.5', '2.2', '0.9', '3.2', '5.4', '6.3', '14.0', '12.7']
sds = ['4.7', '5.2', '3.1', '3.9', '1.5', '1.1', '1.5', '1.6', '0.8', '2.3', '2.3', '3.0', '3.2', '3.8']
meanrow = next((i for (i, r) in enumerate(cells) if r and r[0] == 'Occlusal force'))
sdrow = next((i for (i, r) in enumerate(cells) if r and r[0] == 'Standard dev.'))
for (j, (tooth, v, sd)) in enumerate(zip([17, 16, 15, 14, 13, 12, 11, 21, 22, 23, 24, 25, 26, 27], means, sds), 2):
    assert cells[meanrow][j].replace(',', '.') == v
    OBS.append({'id': 'FERRATO_' + str(tooth), 'claim_type': 'information_link', 'quantity': 'relative_tooth_force', 'value': v, 'unit': 'percent', 'dispersion': {'kind': 'reported_SD', 'value': sd}, 'statistic': 'reported_group_mean', 'resolution_level': 'POPULATION', 'measurement_support_resolution': 'PER_TOOTH', 'population': '28healthycontrolsubjects;upperteeth17..27 excludingM3', 'protocol': 'T-ScanIIImaximumintercuspation;Table3valuespreservedwithoutrenormalisation;sum93percentconflictsTable5sidetotal99.9', 'context': {'FDI': tooth, 'arch': 'maxillary'}, 'source_kind': 'independent_measurement', 'reference': {**metadata('ferrato.html'), 'locator': f'Table3[t3-79-88]controlgroup tooth{tooth},mean/SD', 'selector': {'kind': 'html_section_cell', 'section_id': 't3-79-88', 'table_index': 0, 'row': meanrow, 'col': j}, 'source_cell': cells[meanrow][j]}, 'edges': ['D-E-TANDLAST-PARTIAL-REFERENCE', 'D-E-K16'], 'target_variable': None, 'transfer_status': 'UNKNOWN_TARGET_PROTOCOL', 'measurement_truth_note': 'Numeric table verified;normalisationclosureisUNKNOWN,notperpatientjointforce.'})
for (load, v, anchor) in [('0-0.5', '0.12', '0.12 MPa under load ranging from 0 to 0.5 N'), ('0.5-1', '0.25', '0.25 MPa within the range of 0.5-1.0 N'), ('1-1.5', '0.44', '0.44 MPa under load 1.0-1.5 N')]:
    abstract('YOSHIDA_' + load, 'YoshidaTitle.json', anchor, v, 'MPa', 'pdl_E', 'Humanmaxillarycentralincisor;subjectcountUNKNOWNinabstract', 'Forcepassingthroughcenterofresistance;magneticsensingdisplacements;loadbin' + load + 'N;modulusderivedfrommeasuredforce-displacement', 'PER_TOOTH', ['D-E-K13', 'D-E-K17'], 'pdl_E')
    OBS[-1]['reference']['doi'] = '10.1016/s1350-4533(01)00073-x'
    OBS[-1]['reference']['url'] = 'https://doi.org/10.1016/s1350-4533(01)00073-x'
root = ET.parse(S / 'PMC3607726.xml').getroot()
para = next((p for p in root.findall('.//body//p') if 'S1 (cone angle, 12°) showed a median deviation' in clean(p)))
for (system, value, angle) in [('S1', '3', '12'), ('S2', '2', '15'), ('S3', '2', '5.7'), ('S4', '4', '7.5')]:
    OBS.append({'id': 'RACK_' + system, 'claim_type': 'information_link', 'quantity': 'abutment_vertical_repositioning_median', 'value': value, 'unit': 'um', 'dispersion': {'kind': 'quartiles in primary paragraph', 'value': None}, 'statistic': 'reported_group_median', 'resolution_level': 'POPULATION', 'measurement_support_resolution': 'PER_TOOTH', 'population': 'Four implant systems, six assemblies/system, two test persons', 'protocol': 'Repeated abutment repositioning and hand tightening; cone angle ' + angle + ' degrees; bidirectional vertical deviations', 'context': {'system': system, 'cone_angle_deg': angle}, 'source_kind': 'independent_measurement', 'reference': {**metadata('PMC3607726.xml'), 'locator': 'Results, Vertical displacement of the abutment, paragraph starting S1 (cone angle, 12°); medians, not Table3 extrema', 'selector': {'kind': 'xml_paragraph_value', 'anchor': 'S1 (cone angle, 12°) showed a median deviation', 'system': system}, 'source_cell': clean(para)}, 'edges': ['D-E-K20', 'D-E-K38'], 'target_variable': None, 'transfer_status': 'UNKNOWN_TARGET_PROTOCOL', 'measurement_truth_note': 'Hand-tightening repositioning median is distinct from clamping preload and manufactured assembly gap.'})
if __name__ == '__main__':
    (P / 'VERIFIED_OBSERVATIONS.jsonl').write_text(''.join((json.dumps(x, ensure_ascii=False) + '\n' for x in OBS)))
    print('primary numeric observations', len(OBS), 'sources', len(set((x['reference']['doi'] for x in OBS))))
