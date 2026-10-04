"""Rebuild compact, locator-bound published arms from local sources; no network."""
from dental_release.paths import expand as _release_expand
import json, hashlib, pathlib, re, xml.etree.ElementTree as ET
from collections import Counter
ROOT = pathlib.Path(__file__).resolve().parents[1]
CORPUS = pathlib.Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
LEGACY = ROOT / 'raw/legacy_source_verification'

def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(pmc):
    p = CORPUS / (pmc + '.xml')
    if not p.exists():
        p = LEGACY / (pmc + '.xml')
    return (ET.parse(p), p)

def tables(x):
    return {t.get('id'): [[' '.join(v.itertext()).strip() for v in tr] for tr in t.findall('.//tr')] for t in x.findall('.//table-wrap')}

def scalar(s):
    return float(re.search('[-+]?\\d+(?:\\.\\d+)?', s).group())

def numbers(s):
    return [float(v) for v in re.findall('[-+]?\\d+(?:\\.\\d+)?', s)]

def main():
    rows = []
    sources = {}

    def source(pmc):
        (x, p) = read(pmc)
        sources[pmc] = {'path': str(p), 'sha256': sha(p), 'doi': next((a.text for a in x.findall('.//article-id') if a.get('pub-id-type') == 'doi'), None), 'license': ' '.join(x.find('.//license').itertext()) if x.find('.//license') is not None else 'See publisher', 'tables': list(tables(x))}
        return x

    def add(id, pmc, family, design, D, L, Df, rho, t, rho_c, T, SD, locator, **extra):
        rows.append(dict(id=id, pmc=pmc, study=pmc, family=family, design=design, D_mm=D, L_mm=L, Df_mm=Df, rho_gcc=rho, cortex_mm=t, cortex_rho_gcc=rho_c, torque_Ncm=T, sd_Ncm=SD, density_kind='manufacturer_nominal_PU_density', observation_level='POPULATION', geometry_level='PER_TOOTH', cortex_level='PER_SURFACE_REGION', model_level='PHENOMENOLOGICAL', locator=locator, primary=True, pitch_mm=None, depth_mm=None, **extra))
    source('PMC3596632')
    b = json.loads((LEGACY / 'Bayarchimeg_table_transcription.json').read_text())
    for (tb, values) in b['rows'].items():
        for (k, (z, T, sd)) in enumerate(values):
            rho = 0.016 * z if tb in ('T1', 'T4') else 0.32
            Df = z if tb in ('T2', 'T5') else 3.6
            t = z if tb == 'T3' else 1.5 if tb in ('T4', 'T5') else 0.0
            add(f'B_{tb}_{z:g}', 'PMC3596632', 'Seoul_Osstem', 'SSII', 4.1, 11.5, Df, rho, t, 0.8 if t else 0.0, T, sd, f"PMC3596632#{tb}; publisher PDF p.{b['pages'][tb]}, row {k + 1}, insertion torque average/SD", n=b['n_table'], cortex_material='short_fiber_epoxy_reported_grade_50', cortex_density_note='Nominal #50 converted at 0.016 g/cm3 per pcf. Epoxy identity vs PCF label ambiguity retained; not measured density.', bore_profile='cylindrical', verification='prior K26 numbers independently checked in publisher PDF')
    x = source('PMC8431239')
    tb = tables(x)['ijerph-18-09238-t001']
    mean = tb[3][1:]
    sd = tb[4][1:]
    for j in range(8):
        density = 0.16 if j < 4 else 0.32
        t = 1.0 if j // 2 % 2 else 0.0
        cyl = j % 2 == 0
        add(f'T_{j}', 'PMC8431239', 'Chieti_AoN', 'Cyroth' if cyl else 'IsFour', 4 if cyl else 4.2, 10.0, 3.4 if cyl else None, density, t, 0.48 if t else 0.0, scalar(mean[j]), scalar(sd[j]), f'PMC8431239#ijerph-18-09238-t001; tr[4]/td[{j + 2}] mean; tr[5] SD', n=10, bore_profile='cylindrical' if cyl else 'tapered_conflicting_4.0_vs_4.2', verification='XML cell', rejection=None if cyl else 'Conical bore: abstract underpreparation 0.2 mm vs methods final drill 4.2 mm; no unique scalar Df')
        if not cyl:
            rows[-1]['primary'] = False
    x = source('PMC7216137')
    tb = tables(x)['ijerph-17-02617-t001']
    for j in range(6):
        cyl = j % 3 == 2
        add(f'C_{j}', 'PMC7216137', 'Chieti_AoN', 'RBM' if cyl else 'UNII' if j % 3 == 0 else 'UNIII', 4.1 if cyl else 4.0, 10.0, 3.0 if cyl else None, 0.16 if j < 3 else 0.32, 1.0, 0.48, scalar(tb[2][j + 1]), scalar(tb[3][j + 1]), f'PMC7216137#ijerph-17-02617-t001; tr[3]/td[{j + 2}] mean; tr[4] SD', n=10, bore_profile='cylindrical' if cyl else 'conical_3.5', verification='XML cell', rejection=None if cyl else 'Conical drill nominal diameter is not full bore profile')
        if not cyl:
            rows[-1]['primary'] = False
    x = source('PMC7071536')
    body = ' '.join(x.getroot().itertext())
    for (design, obs, sd) in [('EC', [12.38, 33.52], [2.2, 6.29]), ('Vega', [14.01, 20.47], [2.47, 3.74])]:
        for j in range(2):
            assert str(obs[j]) in body
            add(f'V_{design}_{j}', 'PMC7071536', 'Oviedo_Klockner', design, 4.0, 10.0, 3.6, 0.32, float(j), 1.64 if j else 0.0, obs[j], sd[j], f"PMC7071536#T{(1 if design == 'EC' else 2)}; TORQUE with/without cortical; same values verified in Results paragraph", n=20, bore_profile='cylindrical_with_profile_bur', verification='XML Results numeric paragraph and table caption; table is image', contact_note='EC has 1.5 mm suprabony polished neck; nominal L retained for frozen area predictor; actual contact differs')
    x = source('PMC7150554')
    tb = tables(x)['T1']
    for (j, row) in enumerate(tb[3:7]):
        v = numbers(row[1])
        t = float(j % 2)
        add(f'G_{j}', 'PMC7150554', 'Nagasaki_Nobel', 'NobelReplaceStraight', 4.3, 9.25, 2.8, 0.16 if j < 2 else 0.32, t, 0.8 if t else 0.0, v[0], v[1], f'PMC7150554#T1; tr[{j + 4}]/td[2], cylindrical ITV mean/SD', n=15, bore_profile='step_2.4_2.8_with_cortical_counterbore' if t else 'step_2.4_2.8', verification='XML cell', cortex_counterbore_D_mm=4.3 if t else None, displacement_um=numbers(row[4])[0], displacement_sd_um=numbers(row[4])[1])
    x = source('PMC10384017')
    tb = tables(x)['materials-16-05064-t001']
    mean = next((i for (i, r) in enumerate(tb) if r and r[0] == 'Mean'))
    for j in range(12):
        design = ['TAC', 'INTRALOCK', 'CYROTH'][j % 3]
        t = 1.0 if j // 3 % 2 else 0.0
        add(f'P23_{j}', 'PMC10384017', 'Chieti_AoN', design, 4.0, 15.0, 3.2 if j % 3 != 1 else None, 0.16 if j < 6 else 0.32, t, 0.48 if t else 0.0, scalar(tb[mean][j + 1]), scalar(tb[mean + 1][j + 1]), f'PMC10384017#materials-16-05064-t001; tr[{mean + 1}]/td[{j + 2}] IT mean; next tr SD', n=10, bore_profile='cylindrical' if j % 3 != 1 else 'conical_4.0', verification='XML cell', contact_note='4-mm block, +/-1-mm sheet; 15-mm implant inserted 2 mm below surface at 15-20 degrees. Nominal L is not contact length.', rejection=None if j % 3 != 1 else 'Conical bore profile not recoverable from nominal 4-mm bur')
        if j % 3 == 1:
            rows[-1]['primary'] = False
    seen = {}
    duplicates = []
    for r in rows:
        key = tuple((r[k] for k in ['study', 'design', 'D_mm', 'L_mm', 'Df_mm', 'rho_gcc', 'cortex_mm', 'torque_Ncm']))
        if r['primary'] and key in seen:
            r['primary'] = False
            r['rejection'] = 'Repeated same experimental condition/mean across tables'
            r['duplicate_of'] = seen[key]
            duplicates.append(r['id'])
        elif r['primary']:
            seen[key] = r['id']
    dump('raw/arms.json', rows)
    dump('raw/source_manifest.json', sources)
    human = []
    x = source('PMC10563580')
    for row in tables(x)['T2'][1:]:
        hu = numbers(row[2])
        tor = numbers(row[3])
        human.append(dict(pmc='PMC10563580', locator='T2; ' + row[0], site=row[0], n=int(row[1]), gray_reported_as_HU=hu[0], gray_sd=hu[1], torque_Ncm=tor[0], sd_Ncm=tor[1], observation_level='POPULATION', primary=False, rejection='Regional marginals; implant D/L, drill and cortex not paired to each torque; CBCT HU not independently calibrated'))
    x = source('PMC10847140')
    for ti in ['Tab3', 'Tab4']:
        for row in tables(x)[ti][2:]:
            human.append(dict(pmc='PMC10847140', locator=ti + '; ' + row[0], regional_BMDC=[numbers(row[1]), numbers(row[2])], torque_categories=['<35', '>=35'], observation_level='POPULATION', primary=False, rejection='Ordinal AI density coefficient and categorized torque; no continuous paired torque or cortex/drill table'))
    x = source('PMC10829676')
    tab = tables(x)['tbl1']
    human.append(dict(pmc='PMC10829676', locator='tbl1; MITV mean, FITV mean', kind='PU_anatomical_mandible', peak_torque_Ncm=[scalar(s) for s in tab[2][2:]], final_torque_Ncm=[scalar(s) for s in tab[4][2:]], primary=False, rejection='Pose changes regional cortical contact; cortex thickness and final bore profile not tabulated', observation_level='POPULATION'))
    dump('raw/other_evidence.json', human)
    dump('raw/source_manifest.json', sources)
    summary = {'candidate_arms': len(rows), 'primary_arms': sum((r['primary'] for r in rows)), 'rejected_arms': sum((not r['primary'] for r in rows)), 'rejection_reasons': dict(Counter((r['rejection'] for r in rows if not r['primary']))), 'duplicates': duplicates, 'studies': sorted({r['study'] for r in rows if r['primary']}), 'families': sorted({r['family'] for r in rows if r['primary']}), 'human_paired_complete_arms': 0, 'thread_pitch_depth_complete_families': 0}
    dump('raw/extraction_summary.json', summary)
    print(json.dumps(summary, indent=2))
if __name__ == '__main__':
    main()
