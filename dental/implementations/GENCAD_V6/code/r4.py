from dental_release.paths import expand as _release_expand
from common import *
import xml.etree.ElementTree as ET, itertools, collections, copy
SOURCE = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext/PMC10958416.xml'))

def table():
    root = ET.parse(SOURCE).getroot()
    tab = root.findall('.//table-wrap')[0]
    grid = []
    pending = {}
    for tr in tab.findall('.//tr'):
        row = {}
        next_pending = {}
        for (k, (left, text)) in pending.items():
            row[k] = text
            if left > 1:
                next_pending[k] = (left - 1, text)
        col = 0
        for cell in tr:
            while col in row:
                col += 1
            txt = ' '.join(''.join(cell.itertext()).split())
            rs = int(cell.get('rowspan', '1'))
            cs = int(cell.get('colspan', '1'))
            for k in range(col, col + cs):
                row[k] = txt
                if rs > 1:
                    next_pending[k] = (rs - 1, txt)
            col += cs
        grid.append([row.get(k, '') for k in range(20)])
        pending = next_pending
    teeth = list(map(int, grid[1][4:]))
    rows = []
    for (i, g) in enumerate(grid[2:]):
        rows.append(dict(source_row=i + 2, day=g[0], scale=g[1], method=g[2], time=g[3], values={str(t): v for (t, v) in zip(teeth, g[4:])}))
    return rows

def run():
    st = time.perf_counter()
    pr = read(ROOT / 'PREREG_R4.json')
    assert sha(SOURCE) == pr['source_sha256']
    tt = table()
    assert len(tt) == 30 and all((len(r['values']) == 16 for r in tt))
    forces = {(r['day'], r['time']): r for r in tt if r['method'] == 'T-Scan'}
    paired = []
    for r in tt:
        if not r['method'].startswith(('Silicone', 'Wax', 'Occlusal Foil')):
            continue
        force = forces[r['day'], r['time']]
        for (tooth, value) in r['values'].items():
            paired.append(dict(tooth=int(tooth), method=r['method'], day=r['day'], time=r['time'], count=int(value), force_pp=float(force['values'][tooth]), count_source_row=r['source_row'], force_source_row=force['source_row']))
    collisions = []
    eligible = 0
    for ((tooth, method), group) in itertools.groupby(sorted(paired, key=lambda r: (r['tooth'], r['method'])), key=lambda r: (r['tooth'], r['method'])):
        for (a, b) in itertools.combinations(list(group), 2):
            eligible += 1
            if a['count'] == b['count']:
                collisions.append(dict(tooth=tooth, method=method, count=a['count'], identity_error=0, force_difference_pp=abs(a['force_pp'] - b['force_pp']), A=a, B=b))
    maximum = max(collisions, key=lambda r: r['force_difference_pp'])
    source_probe = next((r for r in paired if r['tooth'] == 16 and r['method'].startswith('Silicone') and (r['day'] == 'Initial') and (r['time'] == '9 a.m.')))

    def valid(records):
        for r in records:
            src = tt[r['count_source_row'] - 2]
            f = tt[r['force_source_row'] - 2]
            if r['count'] != int(src['values'][str(r['tooth'])]) or r['force_pp'] != float(f['values'][str(r['tooth'])]):
                return False
        return True
    mutant = copy.deepcopy(paired)
    mutant[0]['force_pp'] += 1
    controls = dict(source_cell_nominal=source_probe['count'] == 8 and source_probe['force_pp'] == 1.3, all_output_cells_match=valid(paired), changed_force_rejected=not valid(mutant))
    refs = read(ROOT.parent / 'LANE_X2_OCCLUSION_B2B/raw/EXTERNAL_REFERENTS.json')['studies']
    (ferr, hatt) = refs[:2]
    ref_upper = {str(k): v for (k, v) in zip(ferr['fdi'], ferr['mean_percent'])}
    ref_lower = {str(10 * q + i + 1): v for q in (3, 4) for (i, v) in enumerate(hatt['single_side_mean_percent'])}
    r3 = read(ROOT / 'raw/R3_ROWS.json')
    external = []
    for (jaw, ref) in [('upper', ref_upper), ('lower', ref_lower)]:
        good = [r for r in r3 if 'classes' in r]
        hits = []
        width = []
        for r in good:
            for (t, v) in ref.items():
                (lo, hi) = r['classes'][t]['share_hull_pp']
                hits.append(lo <= v <= hi)
                width.append(hi - lo)
        external.append(dict(jaw=jaw, source=ferr['doi'] if jaw == 'upper' else hatt['doi'], quantity='per-tooth relative force share', comparisons=len(hits), covered=sum(hits), mean_hull_width_pp=float(np.mean(width)), reference_sum_pp=sum(ref.values()), joint_normalization='FAIL_93_PERCENT_NOT_RENORMALIZED' if jaw == 'upper' else '100_PERCENT_MARGINAL_COVERAGE_NOT_JOINT_FEASIBILITY', patient_validation='UNKNOWN'))
    native = read(ROOT / 'raw/R3_NATIVE.json')
    area = []
    for (jaw, ref) in [('upper', ref_upper), ('lower', ref_lower)]:
        for r in native:
            shares = r['area_control']['per_tooth_share_pp']
            vals = [abs(shares[t] - v) for (t, v) in ref.items() if shares[t] is not None]
            area.append(dict(case_key=r['case_key'], jaw=jaw, MAE_pp=float(np.mean(vals)) if vals else None, quantity='uniform-pressure area-share practice versus unmatched population profile', resolution='PER_TOOTH -> POPULATION'))
    r1 = read(ROOT / 'raw/R1_ROWS.json')
    norm = []
    for m in FRONTIER:
        candidates = [r for r in r1 if r['participant'] == m and r['family'] == 'lattice_onlay' and (r['level'] == 'normal')]
        rr = [r for r in candidates if r.get('contact', {}).get('status') == 'SCORED']
        rate = sum((r['contact']['predicted']['count'] == 0 for r in rr)) / len(rr) if rr else None
        norm.append(dict(participant=m, tooth_fdi=37, tooth_type='mandibular_second_molar', requested=len(candidates), scored=len(rr), zero_geometric_contact_fraction=rate, published_zero_contact_fraction=0.2, descriptive_difference=rate - 0.2 if rate is not None else None, locator='https://doi.org/10.1016/j.aanat.2023.152112', admission='DESCRIPTIVE_ONLY_GEDAS_SILICONE_VS_PL_GEOMETRY_AND_SELECTED_SUBSET'))
    counts = []
    for r in read(ROOT / 'raw/R2_ROWS.json'):
        if r['status'] == 'SCORED' and r['contact']['status'] == 'SCORED':
            counts.append(dict(uid=r['uid'], fdi=r['source_fdi'], contact_regions=r['contact']['predicted']['count'], area_mm2=r['contact']['predicted']['area_mm2'], resolution='PER_TOOTH'))
    result = dict(source=str(SOURCE), source_sha256=sha(SOURCE), external_referent=pr['external_referent'], table_rows=len(tt), paired_records=len(paired), eligible_session_pairs=eligible, equal_count_pairs=len(collisions), equal_count_unequal_force_pairs=sum((r['force_difference_pp'] > 0 for r in collisions)), largest_reported_difference=maximum, controls=controls, force_reference_comparison=external, area_practice=area, posterior_norm_comparison=norm, source_admission=dict(candidates=5, individual_threshold_admitted=0, rejected=5, rejection_fraction=1.0, reasons=['silicone/foil/PL method mismatch', 'population versus individual criterion', 'pooled tooth types', 'wide natural variation', 'same-tooth repeat table is not a population norm']), seconds=time.perf_counter() - st)
    result['extraction_flow'] = dict(source_tooth_entries=480, contact_count_entries=288, force_lookup_entries=96, qualitative_shimstock_excluded=96, exclusion_fraction=96 / 480, reason='Ordinal +/minus observations are not numeric contact counts', paired_records=288, unpaired_contact_entries=0, missing_unit_or_locator=0)
    dump(ROOT / 'raw/EXTERNAL_PAIRED_OBSERVATIONS.json', paired)
    dump(ROOT / 'raw/EXTERNAL_TABLE_ROWS.json', tt)
    dump(ROOT / 'raw/EXTERNAL_COLLISIONS.json', collisions)
    dump(ROOT / 'raw/CONTACTS_PER_TOOTH.json', counts)
    dump(ROOT / 'rounds/R4.json', result)
    state('R4_DECIDED', str(result['equal_count_unequal_force_pairs']) + ' equal-count pairs differ in reported force', 'Bind external observations, validate outputs and freeze lab-facing predictions')
    return result
if __name__ == '__main__':
    print(clean(run()))
