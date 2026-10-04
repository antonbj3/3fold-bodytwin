"""Extract specified source cells; all values come from local JATS, never constants."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import hashlib, json, re, xml.etree.ElementTree as ET
LANE = Path(__file__).resolve().parents[1]
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
CACHE = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/datasets/dental_3fold_extra/scratch_MICRO/fulltext'))

def text(e):
    return ' '.join(''.join(e.itertext()).split()) if e is not None else ''

def grid(table):
    """Expand JATS rowspan/colspan, keeping exact source-cell indices."""
    (out, spans, loc) = ([], {}, [])
    for (ri, tr) in enumerate(table.findall('.//tr')):
        (row, locations, col) = ([], [], 0)
        cells = list(tr)
        for (ci, cell) in enumerate(cells):
            while col in spans:
                (val, remaining, source) = spans.pop(col)
                row.append(val)
                locations.append(source)
                if remaining > 1:
                    spans[col] = (val, remaining - 1, source)
                col += 1
            value = text(cell)
            (rs, cs) = (int(cell.get('rowspan', '1')), int(cell.get('colspan', '1')))
            for _ in range(cs):
                source = {'tr': ri + 1, 'cell': ci + 1}
                row.append(value)
                locations.append(source)
                if rs > 1:
                    spans[col] = (value, rs - 1, source)
                col += 1
        while col in spans:
            (val, remaining, source) = spans.pop(col)
            row.append(val)
            locations.append(source)
            if remaining > 1:
                spans[col] = (val, remaining - 1, source)
            col += 1
        out.append(row)
        loc.append(locations)
    return (out, loc)

def numbers(s):
    return [float(v.replace('−', '-')) for v in re.findall('[−-]?\\d+(?:\\.\\d+)?', s)]

def study(pmc):
    path = CORPUS / (pmc + '.xml')
    if not path.exists():
        path = CACHE / (pmc + '.xml')
    r = ET.parse(path).getroot()
    doi = next((a.text for a in r.findall('.//article-id') if a.get('pub-id-type') == 'doi'), None)
    return ({'pmcid': pmc, 'doi': doi, 'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'title': text(r.find('.//article-title')), 'license': text(r.find('.//permissions'))}, r)

def main():
    (rows, sources, tables, notes) = ([], {}, {}, [])

    def prepare(pmc, table_id):
        (src, r) = study(pmc)
        sources[pmc] = src
        t = next((t for t in r.findall('.//table-wrap') if t.get('id') == table_id))
        (g, l) = grid(t.find('table'))
        tables[pmc + '#' + table_id] = {'grid': g, 'cell_indices': l}
        return (r, g, l)

    def add(pmc, tid, g, l, ri, ci, condition, force, obs, **kw):
        v = numbers(g[ri][ci])
        source = sources[pmc]
        d = {'id': pmc + ':' + condition, 'study': pmc, 'doi': source['doi'], 'source_path': source['path'], 'source_sha256': source['sha256'], 'table_id': tid, 'table_locator': pmc + '#' + tid, 'source_cell': l[ri][ci], 'source_text': g[ri][ci], 'condition': condition, 'motion_um': v[0], 'sd_um': v[1] if len(v) > 1 else None, 'unit': 'um', 'force_N': force, 'force_provenance': pmc + ' body/Methods (specified load, NOT table motion cell)', 'observation_operator': obs, 'quantity_kind': 'independent_measurement', 'resolution': 'POPULATION', 'spatial_observation_level': 'PER_POINT', 'interface_relative_direct': False, 'time_scale': 'SIMULTANEOUS', 'bone_class_is_human': False, 'interface_threshold_verdict': 'UNKNOWN: observation is not paired local interface motion', **kw}
        rows.append(d)
        return d
    (pmc, tid) = ('PMC5577443', 'T2')
    (r, g, l) = prepare(pmc, tid)
    ri = next((i for (i, x) in enumerate(g) if 'Mean' in x))
    si = next((i for (i, x) in enumerate(g) if 'SD' in x))
    for (j, (angle, tilt)) in enumerate([(0, 0), (0, 30), (45, 0), (45, 30)]):
        ci = len(g[ri]) - 4 + j
        d = add(pmc, tid, g, l, ri, ci, f'angle{angle}_tilt{tilt}', 200, 'abutment_crosshead', rho_g_cm3=0.32, cortex_mm=2, D_mm=4.3, L_mm=10, e_mm=7 if tilt == 0 else 5, angle_deg=angle, implant_tilt_deg=tilt, n=3, implant='NobelReplace Tapered Groovy', substrate='PU foam + short-fibre epoxy', IT_Ncm=55.8 if tilt == 0 else 48.5, IT_provenance=pmc + ' Results: mean maximum ITV; NOT independent specimen pairing')
        d['sd_um'] = numbers(g[si][-4 + j])[0]
        d['sd_source_cell'] = l[si][-4 + j]
    (pmc, tid) = ('PMC7150554', 'T1')
    (r, g, l) = prepare(pmc, tid)
    for ri in range(3, len(g)):
        cls = g[ri][0]
        if not re.search('Low', cls):
            continue
        rho = 0.32 if 'medium' in cls else 0.16
        cortex = 1 if '1 mm' in cls else 0
        for (ci, shape) in [(4, 'cylindrical'), (5, 'tapered')]:
            d = add(pmc, tid, g, l, ri, ci, f'{shape}_rho{rho}_cortex{cortex}', 15, 'abutment_crosshead', rho_g_cm3=rho, cortex_mm=cortex, D_mm=4.3, L_mm=9.25 if shape == 'cylindrical' else 10.54, e_mm=7, angle_deg=90, implant_tilt_deg=0, n=15, implant='NobelReplace ' + shape + ' Groovy', substrate='PU foam; 0.80 g/cm3 cortical surrogate when present', shape=shape)
            it = numbers(g[ri][1 if shape == 'cylindrical' else 2])
            d['IT_Ncm'] = it[0]
            d['IT_sd_Ncm'] = it[1]
            d['IT_source_cell'] = l[ri][1 if shape == 'cylindrical' else 2]
    (pmc, tid) = ('PMC8705369', 'materials-14-07886-t001')
    (r, g, l) = prepare(pmc, tid)
    ri = next((i for (i, x) in enumerate(g) if x and x[0] == 'Micromov X'))
    tri = next((i for (i, x) in enumerate(g) if x and x[0] == 'Torque'))
    for (ci, mat) in [(3, 'Ti'), (6, 'Zr')]:
        d = add(pmc, tid, g, l, ri, ci, mat, 50, 'optical_collar_X_lab_reference', D_mm=3.8, L_mm=8, angle_deg=None, load_angle_source='6 degrees incidence to occlusal plane; not silently mapped to implant axis', n=21, implant='Klockner SK2 ' + mat, substrate='preserved then rehydrated bovine rib, Lekholm/Zarb type III', biological_class='bovine type III', sensor_resolution_um=2, rho_g_cm3=None, cortex_mm=None)
        it = numbers(g[tri][ci])
        d['IT_Ncm'] = it[0]
        d['IT_sd_Ncm'] = it[1]
    (pmc, tid) = ('PMC10792929', 'Tab1')
    (r, g, l) = prepare(pmc, tid)
    for (ri, x) in enumerate(g):
        if len(x) < 6 or x[3] not in ['B11', 'C11', 'PT', 'PZ', 'MT', 'MZ']:
            continue
        subgroup = x[3]
        for (ci, coordinate) in [(4, 'Z'), (5, 'total')]:
            add(pmc, tid, g, l, ri, ci, subgroup + '_' + coordinate, 50, 'optical_HexMeS_' + coordinate, D_mm={'B11': 4.5, 'C11': 5.5}.get(subgroup), L_mm=11 if subgroup in ['B11', 'C11'] else None, angle_deg=0, implant_tilt_deg=0, n=5, implant='Ankylos ' + subgroup if subgroup in ['B11', 'C11'] else 'root analogue ' + subgroup, rho_g_cm3=0.16, cortex_mm=2, substrate='PU foam + epoxy', subgroup=subgroup, fixation='screwed' if subgroup in ['B11', 'C11'] else 'resin layer between RAI and foam', moment_Nmm=None, moment_note='spoon load generates unreported lever moment; nominal axial force is not pure axial wrench')
    (pmc, tid) = ('PMC9917810', 'jcm-12-00855-t001')
    (r, g, l) = prepare(pmc, tid)
    notes.append({'pmcid': pmc, 'table_id': tid, 'status': 'EXCLUDED_ABSOLUTE_MOTION', 'reason': 'Table equation M has no unit; methods give optical sensitivity in um, not equation ordinate unit. Bone reference locations and regression residuals absent.', 'equations': [{'implant': x[0], 'equation': x[1], 'correlation_coefficient': float(x[2])} for x in g[1:]], 'load_range_N': [20, 200], 'load_angle_deg': 30, 'D_mm': 4, 'L_mm': 12, 'resolution': 'POPULATION', 'bone_quality': 'bovine class 2; hardness MPa is NOT bone modulus'})
    inventory = json.loads((LANE / 'raw/candidate_inventory.json').read_text())
    additional = (LANE / 'raw/additional_candidates.txt').read_text().splitlines()
    paths = {x['path'] for x in inventory} | set(additional) | {str(p) for p in CACHE.glob('*.xml')}
    ledger = []
    for path in sorted(paths):
        pmc = Path(path).stem
        if pmc in sources and pmc != 'PMC9917810':
            reason = 'RETAINED_MEASURED_PROXY_NOT_DIRECT_INTERFACE'
        elif pmc == 'PMC9917810':
            reason = 'REJECT_ABSOLUTE_UNIT_AND_REFERENCE_UNRESOLVED'
        else:
            rr = ET.parse(path).getroot()
            title = text(rr.find('.//article-title'))
            body = text(rr.find('body'))
            if pmc == 'PMC11000944':
                reason = 'REJECT_FE_MICROMOTION; separate displacement table has malformed 6026 FEA cell and independence unresolved (near-duplicate 2017 setup)'
            elif pmc in ['PMC10225166', 'PMC10261134']:
                reason = 'REJECT_NON_DENTAL_HIP_OR_BONE_CANCER_IMPLANT; eligible-looking motion table has wrong domain'
            elif pmc in ['PMC10517556', 'PMC10925053']:
                reason = 'REJECT_FE_OUTPUT_NOT_INDEPENDENT_MEASUREMENT; explicit three-dimensional finite-element study'
            elif rr.get('article-type') == 'review-article' or re.search('systematic review|scoping review|meta-analysis', title, re.I):
                reason = 'REJECT_SECONDARY_SOURCE'
            elif not re.search('dental|oral|implant', title, re.I):
                reason = 'REJECT_WRONG_DOMAIN_OR_NO_DENTAL_MOTION_TABLE'
            elif re.search('abutment|connection', title, re.I) and (not re.search('bone|primary stability|micromotion', title, re.I)):
                reason = 'REJECT_WRONG_INTERFACE'
            elif not any((re.search('(displacement|micromotion|micromovement).*(µm|μm|micron)', text(t), re.I) for t in rr.findall('.//table-wrap'))):
                reason = 'REJECT_NO_ELIGIBLE_ABSOLUTE_LOADED_MOTION_TABLE'
            else:
                reason = 'REJECT_WRONG_QUANTITY_LOAD_OR_DOMAIN; see local source, not positive interface evidence'
        ledger.append({'path': path, 'pmcid': pmc, 'screening_decision': reason})
    for (name, data) in [('measurements.json', rows), ('source_manifest.json', list(sources.values())), ('tables.json', tables), ('source_exclusions.json', notes), ('screening_ledger.json', ledger)]:
        (LANE / 'raw' / name).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'measured_rows': len(rows), 'primary_motion_studies': len({d['study'] for d in rows}), 'direct_local_interface_studies': 0, 'screened_paths': len(ledger)}))
if __name__ == '__main__':
    main()
