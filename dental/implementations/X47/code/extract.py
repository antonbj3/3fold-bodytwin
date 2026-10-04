"""Source-bound table extraction. Coordinates are zero-based XML rows/cells.
Values come from the local primary tables, never from a fixture.
"""
from dental_release.paths import expand as _release_expand
import json, re, hashlib, xml.etree.ElementTree as ET
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))

def txt(e):
    return ' '.join(''.join(e.itertext()).split()) if e is not None else ''

def plain(e):
    if e.tag == 'sup':
        return ''
    return (e.text or '') + ''.join((plain(c) + (c.tail or '') for c in e))

def numeric(s):
    return [float(x) for x in re.findall('(?<![\\w])[-+]?\\d+(?:\\.\\d+)?', s.replace(',', ''))]

def cells(t, tab):
    e = t.find(f"./body//*[@id='{tab}']")
    if e is None:
        raise ValueError(tab)
    return [[x for x in r if x.tag in ('td', 'th')] for r in e.findall('.//tr')]

def extract():
    manifest = json.loads((ROOT / 'raw/SOURCE_MANIFEST.json').read_text())
    sources = {x['pmc']: x for x in manifest}
    rows = []

    def add(pmc, tab, r, c, n, meta, sd_c=None, skip_prefix=False, sd_r=None):
        t = ET.parse(CORPUS / (pmc + '.xml'))
        cc = cells(t, tab)
        raw = plain(cc[r][c])
        s = raw.split(':', 1)[-1] if skip_prefix else raw
        nums = numeric(s)
        if sd_c is None:
            (mean, sd) = nums[:2]
        else:
            mean = nums[0]
            sd = numeric(plain(cc[r if sd_r is None else sd_r][sd_c]))[0]
        row = {'row_id': f'{pmc}:{tab}:r{r}:c{c}', 'study': pmc, 'doi': sources[pmc]['doi'], 'source_path': sources[pmc]['path'], 'source_sha256': sources[pmc]['sha256'], 'table_id': tab, 'xml_row': r, 'xml_cell': c, 'xml_sd_cell': sd_c, 'xml_sd_row': sd_r, 'strip_label_prefix': skip_prefix, 'table_label': txt(t.find(f"./body//*[@id='{tab}']/label")), 'table_cell_text': txt(cc[r][c]), 'endpoint': 'whole_crown_axial_pull_off', 'unit': 'N', 'mean_force_N': mean, 'sd_force_N': sd, 'n': n, 'height_mm': None, 'total_convergence_deg': None, 'reported_angle_deg': None, 'angle_kind': 'UNKNOWN', 'thermal_cycles': None, 'thermal_low_C': None, 'thermal_high_C': None, 'dwell_s': None, 'transfer_s': None, 'mechanical_cycles': 0, 'mechanical_force_N': None, 'storage_days': None, 'storage_C': None, 'cement': None, 'permanent_resin_cement': 0, 'substrate': 'titanium', 'crown_material': None, 'surface_feature': 'standard', 'surface_treatment': None, 'group': None, 'independent_first_state': True, 'resolution': 'POPULATION', 'specimen_resolution': 'PER_TOOTH', 'timescale': 'HANDOVER', 'covariate_provenance': 'Local body Materials/Methods; explicit unknowns preserved', 'notes': []}
        row.update(meta)
        rows.append(row)
    for (i, group) in enumerate(['Temp', 'TempNE', 'C', 'Temp-S', 'Temp-SU', 'TempNE-S', 'TempNE-SU']):
        add('PMC10028587', 'T2', i + 3, 1, 20 if i < 2 else 10, dict(group=group, height_mm=7, crown_material='zirconia', cement=['TempBond', 'TempBond NE', 'Panavia SA Luting plus', 'Panavia SA Luting plus', 'Panavia SA Luting plus', 'Panavia SA Luting plus', 'Panavia SA Luting plus'][i], permanent_resin_cement=int(i >= 2), thermal_cycles=5000, thermal_low_C=5, thermal_high_C=55, dwell_s=30, storage_days=1, storage_C=37, independent_first_state=i < 3, notes=['Re-cemented groups share prior specimens with initial temporary groups; excluded from independent-state model.'] if i >= 3 else []))
    for (k, dimple) in enumerate([0, 2, 4]):
        for (j, (h, cement, code)) in enumerate([(3, 'TempBond', 'T'), (3, 'PANAVIA F 2.0', 'P'), (6, 'TempBond', 'T'), (6, 'PANAVIA F 2.0', 'P')]):
            add('PMC10154144', 'T1', k + 2, j + 1, 15, dict(group=f'H{h}-{code}{dimple}', height_mm=h, total_convergence_deg=6, reported_angle_deg=6, angle_kind='total_convergence', cement=cement, permanent_resin_cement=int(code == 'P'), crown_material='CoCr', surface_feature=f'{dimple}_dimples', thermal_cycles=10000, thermal_low_C=5, thermal_high_C=55, dwell_s=30, transfer_s=5, notes=['Source says 30 s bath exposure and 5 s dwell; last field is the inter-bath timing, wording ambiguous.']), skip_prefix=True)
    for (c, g) in [(1, 'straight_access'), (2, 'angled_access_24deg')]:
        add('PMC10326962', 'Tab1', 1, c, 7, dict(group=g, crown_material='zirconia', cement='Multilink speed', permanent_resin_cement=1, thermal_low_C=5, thermal_high_C=55, dwell_s=60, mechanical_cycles=250000, mechanical_force_N=100, notes=['Ti-base dimensions 4.5 x 6 x 2.5 mm not explicitly mapped to bonded height; thermal cycle count missing; screw-access angle 24 deg is NOT preparation convergence.']))
    for (c, g) in enumerate(['Control', 'IPS', 'Ceramic', 'Composite'], 1):
        add('PMC10560107', 'tab1', 3, c, 5, dict(group=g, substrate='mixed_post_core', cement='Choice2 Veneer Cement', permanent_resin_cement=1, crown_material='zirconia', reported_angle_deg=6, angle_kind='taper_tool_not_verified_TOC', thermal_cycles=5000, thermal_low_C=5, thermal_high_C=55, dwell_s=30, storage_days=1, storage_C=37, notes=['Separate n=5 fracture and n=5 retention specimens per group; 6 degree tool is not enough to identify total convergence; prepared height unreported.']))
    for (i, g) in enumerate(['ZS', 'ZL', 'LS', 'LL', 'VS', 'VL']):
        add('PMC10594673', 'Tab5', i + 1, 8, 7, dict(group=g, height_mm=4 if i % 2 == 0 else 7, crown_material=['zirconia', 'zirconia', 'lithium_disilicate', 'lithium_disilicate', 'Vita_Enamic', 'Vita_Enamic'][i], cement='Multilink hybrid-abutment', permanent_resin_cement=1, thermal_cycles=5000, thermal_low_C=5, thermal_high_C=55, dwell_s=30, mechanical_cycles=50000, mechanical_force_N=49, storage_days=30, storage_C=37, notes=['Corrected author name only; local correction PMC10848489 read; offset implant placement fixed at 5 mm.']))
    for (c, g) in enumerate(['Control', 'Air abraded', 'Oxidized', 'Oxidized-air abraded'], 1):
        add('PMC10599539', 'pone.0287108.t001', 2, c, 10, dict(group=g, height_mm=5.5, cement='resin cement (final-test brand UNKNOWN)', permanent_resin_cement=1, crown_material='zirconia', surface_treatment=g, thermal_cycles=5000, thermal_low_C=5, thermal_high_C=55, dwell_s=120, storage_days=1, notes=['RelyX Ultimate is named for pre-oxidation bonding in group3; brand of final cementation not explicitly given. Primary article DOI excludes nested peer-review suffix .r005; methods mention 60 crowns but numerical retention table has 4 x 10 specimens.']))
    for (r, g) in [(1, 'Nonanodized'), (2, 'Anodized')]:
        add('PMC10630546', 'T1', r, 2, 12, dict(group=g, crown_material='lithium_disilicate', cement='TempBond NE', permanent_resin_cement=0, thermal_cycles=None, surface_treatment=g, notes=['6 mm total abutment height with 2 mm gingival height: bonded height UNKNOWN. Pilot 2 specimens debonded at 10000 mechanical cycles; test then omitted cyclic loading and re-cemented pilots. Thermal protocol is not specified despite discussion saying thermocycling. Reported n=12 per final group; no aged paired force.']), sd_c=3)
    for (c, g) in [(1, 'Resin nanoceramic'), (2, 'Lithium disilicate')]:
        add('PMC10696722', 'Tab3', 1, c, 10, dict(group=g, crown_material=g, cement='Duo-Link', permanent_resin_cement=1, thermal_cycles=0, storage_days=1, storage_C=37, notes=['6 mm reported analog height is NOT prepared bonded abutment height; no aging performed.']))
    for (r, (g, cement, perm)) in enumerate([('1', 'SpeedCEM Plus', 1), ('2', 'Panavia SA Cement Universal', 1), ('3', 'Panavia V5', 1), ('4', 'RelyX Unicem 2 Automix', 1), ('5', 'VITA ADIVA IA-Cem', 1), ('6', 'Ketac CEM', 0), ('7', 'Hoffmann phosphate', 0)], 1):
        add('PMC10707355', 'materials-16-07314-t005', r, 1, 7 if g == '4' else 8, dict(group=g, crown_material='zirconia', cement=cement, permanent_resin_cement=perm, thermal_low_C=5, thermal_high_C=55, mechanical_cycles=5000000, notes=['TC count and height missing; group 7 includes two premature decementations as zero force; group 4 excludes one abutment/screw fracture after 2778775 cycles.']), sd_c=2)
    for (r, cement) in enumerate(['Harvard zinc phosphate', 'TempBond', 'Dycal', 'Dentotemp', 'Implantlink'], 5):
        for (c, material) in [(1, 'CoCr'), (4, 'zirconia')]:
            add('PMC10758385', 'T2', r, c, 10, dict(group=cement + '_' + material, crown_material=material, cement=cement, thermal_cycles=5000, thermal_low_C=5, thermal_high_C=55, dwell_s=30, storage_days=1, storage_C=37, notes=['8.5 mm specific abutment height with 1.5 mm gingival height: bonded height UNKNOWN; dye immersion adds 24 h.']))
    t = ET.parse(CORPUS / 'PMC10870690.xml')
    cc = cells(t, 't2-medscimonit-30-e943237')
    for (r, crow) in enumerate(cc):
        for (c, e) in enumerate(crow):
            code = txt(e).replace('Gp ', '')
            if re.fullmatch('[MP][CR][12][GB]?', code):
                add('PMC10870690', 't2-medscimonit-30-e943237', r, c + 2, 10, dict(group=code, height_mm=3, total_convergence_deg=10 if code[2] == '1' else 20, reported_angle_deg=10 if code[2] == '1' else 20, angle_kind='total_convergence', substrate='printed_resin_die', crown_material='milled_PMMA' if code[0] == 'M' else 'printed_provisional_resin', cement='TempBond NE' if code[1] == 'C' else 'ProviTemp NE', permanent_resin_cement=0, surface_feature='groove' if code.endswith('G') else 'box' if code.endswith('B') else 'standard', thermal_cycles=10000, thermal_low_C=5, thermal_high_C=55, dwell_s=30, transfer_s=5, notes=['Title calls it stress, but methods explicitly define tensile force in N and Table2 names tensile forces; provisional resin is NOT permanent resin class. Methods mention two crosshead speeds 1 and 5 mm/min without clear group mapping.']), sd_c=c + 3)
    for c in range(1, 9):
        group = ('A' if c <= 4 else 'B') + str((c - 1) % 4 + 1)
        j = (c - 1) % 4
        add('PMC10926708', 't1-medscimonit-30-e943401', 2, c, 10, dict(group=group, height_mm=3, total_convergence_deg=10 if j == 0 else 20, reported_angle_deg=10 if j == 0 else 20, angle_kind='total_convergence', substrate='printed_resin_die', crown_material='zirconia' if c <= 4 else 'NiCr', cement='RIVA glass ionomer', surface_feature=['standard', 'standard', 'groove', 'box'][j], thermal_cycles=10000, thermal_low_C=5, thermal_high_C=55, dwell_s=15, transfer_s=5, notes=['Numeric descriptives in Table1; Tables2/3 explicitly identify N. Text cross-references Table2 for means; value magnitude kept as published.']), sd_c=c, sd_r=3)
    return rows

def ccell(t, tab, r, c):
    return cells(t, tab)[r][c]

def main():
    rows = extract()
    (ROOT / 'raw/RETENTION_GROUPS.json').write_text(json.dumps(rows, indent=2, ensure_ascii=False) + '\n')
    import csv
    with (ROOT / 'raw/RETENTION_GROUPS.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('extracted', len(rows), 'groups', len(set((r['study'] for r in rows))), 'studies')
if __name__ == '__main__':
    main()
