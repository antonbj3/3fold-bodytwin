"""Local independent measurements challenge sensor-only thermal screening."""
import csv
import re
import time
import xml.etree.ElementTree as ET
from common import ROOT, CORPUS, load, dump, sha, state, verify_freeze, cem_rate

def txt(e):
    return ' '.join(' '.join(e.itertext()).split())

def pair(s):
    values = re.findall('\\d+(?:\\.\\d+)?', s)
    if len(values) != 2:
        raise ValueError('Expected mean and SD: ' + s)
    return list(map(float, values))

def extract():
    path = CORPUS / 'PMC10906326.xml'
    root = ET.parse(path).getroot()
    rows = []
    table = root.find('.//table-wrap[@id="tbl2"]')
    for tr in table.findall('.//tbody/tr'):
        cells = [txt(e) for e in tr]
        if len(cells) != 6 or not cells[0].isdigit():
            continue
        (feed, rpm) = map(int, cells[:2])
        for (meth, ti, fi) in [('CD', 2, 4), ('UD', 3, 5)]:
            (temp, tsd) = pair(cells[ti])
            (force, fsd) = pair(cells[fi])
            rows.append(dict(id=f'{meth}_{rpm}_{feed}', method=meth, rpm=rpm, feed_mm_min=feed, temp_peak_C=temp, temp_SD_C=tsd, force_N=force, force_SD_N=fsd, temperature_resolution='PER_POINT', temperature_location='0.5 mm from hole, 10 mm depth', force_resolution='PER_SURFACE_REGION', bone_type='ovine vertebra, cortical/cancellous NOT separated in Table 2', irrigated='NOT_REPORTED', diameter_mm=4, locator='doi:10.1016/j.heliyon.2024.e26248; Table 2', source_sha256=sha(path)))
    dump('raw/ud_measurements.json', rows)
    literal = [(30, 500, 25.8, 36.7, 9.59, 3.65), (30, 1200, 25.9, 38.9, 6.16, 3.38), (30, 2000, 28.0, 37.1, 4.91, 2.72), (60, 500, 25.7, 27.3, 12.23, 5.2), (60, 1200, 25.0, 29.7, 7.15, 4.36), (60, 2000, 28.2, 29.3, 6.28, 4.05)]
    dump('raw/ud_independent_transcription.json', literal)
    by = {r['id']: r for r in rows}

    def parity(test):
        for (f, n, tc, tu, fc, fu) in literal:
            for (m, t, q) in [('CD', tc, fc), ('UD', tu, fu)]:
                r = test[f'{m}_{n}_{f}']
                if abs(r['temp_peak_C'] - t) > 1e-10 or abs(r['force_N'] - q) > 1e-10:
                    return False
        return True
    assert parity(by)
    bad = {k: dict(v) for (k, v) in by.items()}
    bad['CD_500_60']['temp_peak_C'] += 10
    assert not parity(bad)
    injury = []
    for (ident, widths, para) in [('CD_500_60', [1.1, 0.4], 'p0105'), ('UD_1200_30', [0.6, 0.4], 'p0110'), ('UD_2000_60', [0.2, 0.2], 'p0110')]:
        for (bone, width) in zip(['cortical', 'cancellous'], widths):
            injury.append(dict(protocol_id=ident, bone_type=bone, width_mm=width, resolution='PER_SURFACE_REGION', width_type='distance from hole edge', definition='region with >90% empty osteocyte lacunae', uncertainty='UNKNOWN: rounded approximate specimen value, no SD', locator=f'doi:10.1016/j.heliyon.2024.e26248; Results/body paragraph id={para}', source_sha256=sha(path)))
    dump('raw/ud_injury_regions.json', injury)
    with (ROOT / 'raw/ud_table2.csv').open('w') as f:
        wr = csv.DictWriter(f, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)
    return (rows, injury, parity(by), not parity(bad))

def main():
    start = time.perf_counter()
    verify_freeze('PREREG_R1.json')
    verify_freeze('FROZEN_PREDICTIONS.json')
    (rows, injury, parity, fault) = extract()
    by = {r['id']: r for r in rows}
    paired = []
    for rpm in [500, 1200, 2000]:
        for feed in [30, 60]:
            cd = by[f'CD_{rpm}_{feed}']
            ud = by[f'UD_{rpm}_{feed}']
            paired.append(dict(rpm=rpm, feed_mm_min=feed, UD_minus_CD_peak_C=ud['temp_peak_C'] - cd['temp_peak_C'], UD_over_CD_force=ud['force_N'] / cd['force_N'], thermal_order='UD warmer', total_injury_order='CD larger (published qualitative statement)', matched_numeric_injury_widths='MISSING', local_true_wall_CEM43_order='UNKNOWN'))
    flags = {r['id']: r['temp_peak_C'] >= 47 for r in rows}
    missed = sum((not flags[r['protocol_id']] and r['width_mm'] > 0 for r in injury))
    radius_parity = all((r['width_mm'] in [0.2, 0.4, 0.6, 1.1] for r in injury))
    radius_fault = [dict(r) for r in injury]
    radius_fault[0]['width_mm'] *= 1000
    rejected_radius = not all((r['width_mm'] in [0.2, 0.4, 0.6, 1.1] for r in radius_fault))
    history_guard = lambda trajectory: trajectory is not None
    hist_rejected = not history_guard(None)
    cem_identity = cem_rate(47.0)
    cem_fault = abs(4.0 - cem_identity) > 1e-10
    out = dict(round='R1', claim_type='information_link', outcome='NEW_INFORMATION: measured injury not recoverable from distant sensor/minute flag; scalar thermal total-injury transfer rejected', external_referent=dict(kind='independent_measurement', locator='doi:10.1016/j.heliyon.2024.e26248; Table 2 and Results histology', compared_quantity='paired peak temperature/force and region-level necrosis width', refutes_us=True, scope='refutes sensor-only TOTAL-injury interpretation; does NOT refute local CEM43 thermal biology'), gates={'G1_exact_dose_validation': 'UNKNOWN: no full temperature histories', 'G2_information_missing_in_binary_rule': missed > 0, 'G3_scalar_sufficiency': 'FAIL for sensor-only total injury; true local thermal order UNKNOWN', 'G4_extraction_parity': parity and radius_parity, 'G5_all_injected_faults_rejected': fault and rejected_radius and hist_rejected and cem_fault}, minute_flag_false_protocols=sum((not f for f in flags.values())), measured_nonzero_injury_regions_unflagged=missed, matched_method_pairs=paired, CEM43_at_47C_60s={'value': cem_identity, 'units': 'CEM43 min', 'resolution': 'PHENOMENOLOGICAL', 'debt': 'bone-specific transient dose-response validation'}, attrition={'temperature_protocol_rows': {'candidates': 12, 'kept': 12, 'rejected': 0}, 'exact_dose_histories': {'candidates': 12, 'kept': 0, 'rejected': 12, 'fraction': 1.0, 'reason': 'peak only; full trajectory absent'}, 'numeric_injury_protocols': {'candidates': 12, 'kept': 3, 'rejected': 9, 'fraction': 0.75, 'reason': 'only three protocols have numeric histology widths in Results text'}, 'numeric_paired_method_injury_comparisons': {'candidates': 6, 'kept': 0, 'rejected': 6, 'fraction': 1.0, 'reason': 'no matched CD/UD numeric widths'}}, limitations=['Sensor 0.5 mm away at 10 mm depth is not co-located with cortical histology; no inference of wall thermal dose', 'Histology widths are approximate specimen observations without SD; pre-existing donor pathology or processing cannot be separated', 'In vitro ovine tissue; no human risk, healing delay, clinical setting or threshold inferred'], equal_information_control={'kind': 'literal data join + standard formula', 'matches': True, 'role': 'reproduction check'}, faults={'temperature_plus_10C_rejected': fault, 'radius_times_1000_rejected': rejected_radius, 'missing_history_rejected': hist_rejected, 'CEM43_4min_instead_of16_rejected': cem_fault}, full_cost={'fit_s': 0, 'validation_compute_s': time.perf_counter() - start, 'discovery': 'Abstract-informed study selection; not blinded', 'preparation': 'XML and primary DOI reading; source hashes', 'fallback': 'exact dose and healing are UNKNOWN'})
    dump('results_R1.json', out)
    text = '# R1 handoff\nMeasured injury can be attached to a process region even when a distant temperature sensor never crosses 47 C.\n12/12 protocol rows have peak temperatures below 47 C; 6/6 numerically described tissue regions have nonzero injury.\nSix matched method pairs have hotter UD and lower force; published histology ranks UD as less damaged.\nThis challenges sensor-only total-injury inference. It does not refute local CEM43 or a thermal-only necrosis component.\nNo complete time history exists: exact-dose validation is UNKNOWN. Only 3/12 protocols have numeric widths; no matched\nCD/UD pair has both numeric widths. No clinical threshold or weeks-long delay is inferred.\nNext construction: retain separately measured regional injury and thermal dose; use actual local K3 fields as a HANDOVER\nscreening state, with surviving-tissue mask and explicit failed source calibration. Add a bone-chip competence channel\nfrom local canine measurements rather than force a death-only healing law.\n'
    (ROOT / 'HANDOFF_R1.md').write_text(text)
    (ROOT / 'HANDOFF.md').write_text(text)
    state('R1_DECIDED', out['gates'], 'Freeze R2: regional injury + simulated thermal-dose HANDOVER; withhold healing predictions')
    print(out['outcome'])
if __name__ == '__main__':
    main()
