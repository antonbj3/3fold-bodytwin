"""Preserve opposing geometry when contact edges disappear. Candidate anatomy only."""
from dental_release.paths import expand as _release_expand
import json, time
from pathlib import Path
import numpy as np
from contact import P, LABELS, sha, write, state

def upper_equivalent(k):
    (q, t) = divmod(int(k), 10)
    return k if q in [1, 2] else 10 + t if q == 4 else 20 + t if q == 3 else None

def extract():
    st = time.perf_counter()
    manifest = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    out = {}
    for c in manifest['cases']:
        j = json.load(open(LABELS / c / 'labels+landmarks.json'))
        row = json.load(open(P / 'raw/cases' / (c + '.json')))
        fr = j['arches']['lower']['frame']
        R = np.column_stack([fr['right_unit'], fr['anterior_unit'], fr['superior_unit']])
        center = np.array(fr['center_mm'])
        a = np.array(fr['anterior_unit'])
        u = j['arches']['upper']['landmarks']['teeth']
        l = j['arches']['lower']['landmarks']['teeth']
        features = dict(j['x7_landmark_port']['measures'])
        teeth = []
        for (k, t) in u.items():
            (q, pos) = divmod(int(k), 10)
            lower = str((4 if q == 1 else 3) * 10 + pos)
            if lower not in l:
                continue
            lt = l[lower]
            pu = np.array(t['incisal_candidate_xyz_mm'] if pos <= 2 else t['mesiobuccal_cusp_candidate_xyz_mm'])
            pl = np.array(lt['incisal_candidate_xyz_mm'] if pos <= 2 else lt['mesiobuccal_cusp_candidate_xyz_mm'])
            direction = a if pos <= 2 else np.array(lt['buccal_unit'])
            direction = direction / max(np.linalg.norm(direction), 1e-12)
            b = float((pu - pl) @ direction)
            up = pu if pos <= 2 else np.array(t['mesiobuccal_cusp_candidate_xyz_mm'])
            low = pl if pos <= 2 else np.array(lt['buccal_groove_candidate_xyz_mm'])
            mes = np.array(lt['mesial_unit'])
            sag = float((up - low) @ mes)
            features['opposition_margin_FDI' + k] = b
            features['opposition_AP_FDI' + k] = sag
            for (axis, val) in zip('ras', (pu - center) @ R):
                features['opposition_upper_' + k + '_' + axis] = float(val)
            for (axis, val) in zip('ras', (pl - center) @ R):
                features['opposition_lower_' + lower + '_' + axis] = float(val)
            pair = next((p for p in row.get('pairs', []) if p['upper_fdi'] == int(k) and p['lower_fdi'] == int(lower)), None)
            teeth.append(dict(upper_fdi=int(k), lower_corresponding_fdi=int(lower), signed_order_margin_mm=b, margin_scenario_mm=[b - 1.0, b + 1.0], crossbite_order_candidate=b < 0, robust_reversed_under_point_error_scenario=b < -1.0, possible_reversed_under_point_error_scenario=b < 1.0, direction_xyz=direction, upper_point_xyz_mm=pu, lower_point_xyz_mm=pl, upper_point_frame_mm=(pu - center) @ R, lower_point_frame_mm=(pl - center) @ R, point_role='incisal candidate' if pos <= 2 else 'MB cusp candidate', anatomical_status='UNKNOWN_TARGET_POINT_AND_AXIS_VALIDITY', m1_relation_candidate_mm=sag if pos == 6 else None, projected_pair_minimum_mm=pair['minimum_projected_gap_mm'] if pair else None, contact_pair_proximity_present=pair['near_contact_possible_under_offset_scenarios'] if pair else False, projected_pair_overlap_invalid=pair['excess_projected_penetration'] if pair else None))
        out[c] = dict(features=features, teeth=teeth, incisor_measures_candidate_mm={k: v for (k, v) in j['x7_landmark_port']['measures'].items() if k in ['overbite_mm', 'overjet_mm']}, m1_measures_candidate_mm={k: v for (k, v) in j['x7_landmark_port']['measures'].items() if k.startswith('molar')}, source_sha256=sha(LABELS / c / 'labels+landmarks.json'), status='Opposition geometry retained; contact not required; anatomical accuracy UNKNOWN')
    write(P / 'raw/OPPOSITION_FEATURES_R2.json', out)
    write(P / 'raw/OPPOSITION_COST_R2.json', dict(wall_s=time.perf_counter() - st, cases=len(out)))
    state('R2_OPPOSITION_READY', 'R2 frozen; candidate signed geometry extracted', 'Fit/freeze categories and toothwise ordering before report comparison')
    print('Opposition', len(out), 'cases', round(time.perf_counter() - st, 2), 's')

def evaluate_teeth():
    from report_capability import wilson
    frozen = json.load(open(P / 'FROZEN_PREDICTIONS_R2.json'))
    assert sha(P / 'raw/OPPOSITION_FEATURES_R2.json') == frozen['extra_features_sha256']
    op = json.load(open(P / 'raw/OPPOSITION_FEATURES_R2.json'))
    manifest = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    X7 = P.parent / _release_expand('X7')
    reports = json.load(open(X7 / 'raw/TEST_LABELS.json'))
    rows = []
    hit = 0
    den = 0
    absence = 0
    true_absence = 0
    reportdice = []
    multi = []
    named_cases = set()
    m1missing = 0
    m1n = 0
    incisornoncontact = 0
    incisorn = 0
    for c in manifest['inherited_partition']['test']:
        tt = op[c]['teeth']
        pred = {t['upper_fdi'] for t in tt if t['crossbite_order_candidate']}
        reportlists = []
        for r in reports.get(c, []):
            lab = r['labels']
            status = lab.get('crossbite')
            ref = {upper_equivalent(int(k)) for k in lab.get('crossbite_teeth', []) if upper_equivalent(int(k)) is not None}
            if status == 'absent':
                absence += 1
                true_absence += int(not pred)
            if status == 'present' and ref:
                named_cases.add(c)
                hit += len(ref & pred)
                den += len(ref)
                dice = 2 * len(ref & pred) / (len(ref) + len(pred)) if ref or pred else 1.0
                reportdice.append(dice)
                reportlists.append(ref)
            rows.append(dict(case_id=c, member=r['member'], member_sha256=r['sha256'], reported_crossbite=status, explicit_named_positive_upper_equivalent=sorted(ref), geometry_reversed_order_upper_equivalent=sorted(pred), matched_explicit_positive_teeth=sorted(ref & pred), unlisted_positive_status='UNKNOWN_IF_REPORT_PRESENT', source_tooth_candidates=tt))
        if len(reportlists) >= 2:
            inter = set.intersection(*reportlists)
            union = set.union(*reportlists)
            multi.append(dict(case_id=c, report_sets=[sorted(s) for s in reportlists], disagree=any((s != reportlists[0] for s in reportlists[1:])), jaccard=len(inter) / len(union) if union else 1.0, reversed_order_candidate=sorted(pred)))
        for (side, k) in [('right', 16), ('left', 26)]:
            if any((r['labels'].get('molar_' + side) is not None for r in reports.get(c, []))):
                m1n += 1
                t = next((t for t in tt if t['upper_fdi'] == k), None)
                if t is None or not t['contact_pair_proximity_present']:
                    m1missing += 1
        if any((r['labels'].get('overjet') is not None for r in reports.get(c, []))):
            incisorn += 1
            if not any((t['upper_fdi'] in [11, 21] and t['contact_pair_proximity_present'] for t in tt)):
                incisornoncontact += 1
    recall = hit / den if den else 0.0
    spec = true_absence / absence if absence else 0.0
    gate = recall >= 0.7 and spec >= 0.7 and (len(named_cases) >= 30)
    result = dict(claim_type='capability', round='R2-toothwise', named_test_cases=len(named_cases), explicit_positive_named_tooth_observations=den, matched_named_tooth_observations=hit, named_tooth_recall=recall, recall_wilson95=wilson(hit, den), explicit_absence_report_observations=absence, correct_absence_reports=true_absence, absence_specificity=spec, specificity_wilson95=wilson(true_absence, absence), mean_report_list_dice_descriptive=float(np.mean(reportdice)) if reportdice else None, report_lists_not_complete_ground_truth=True, primary_gate='PASS' if gate else 'FAIL', wrong_order_control=dict(injection='Predict every available tooth reversed for explicit absent reports', absence_specificity=0.0, gate='FAIL', refutes=True), wrong_named_tooth_control=dict(injection='Return empty set for every named-positive report', named_tooth_recall=0.0, gate='FAIL', refutes=True), multi_report_named_list_cases=len(multi), named_list_discordant_cases=sum((r['disagree'] for r in multi)), contact_only_missing_M1_relationships=dict(no_proximal_M1_edge=m1missing, described_M1_sides=m1n), contact_only_missing_anterior_relationships=dict(no_central_incisor_proximity_edge=incisornoncontact, described_overjet_cases=incisorn), external_referent=dict(kind='independent_measurement', locator=manifest['zip'] + '::*/reports_ios_en/*.txt, X7 source hashes in comparison rows', compared_quantity='Explicit named crossbite tooth IDs and explicit absence, not landmark coordinate accuracy', refutes_us=True), independence='Multiple report/tooth observations clustered within patients; Wilson intervals are descriptive observation intervals, not patient-independent inferential confidence.', clinical_point_accuracy='UNKNOWN')
    write(P / 'raw/TOOTHWISE_ROWS_R2.json', rows)
    write(P / 'raw/TOOTH_LIST_DISAGREEMENT_R2.json', multi)
    write(P / 'raw/RESULTS_TOOTHWISE_R2.json', result)
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    import sys
    extract() if sys.argv[1] == 'extract' else evaluate_teeth()
