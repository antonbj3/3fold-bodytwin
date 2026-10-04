import argparse
import collections
import datetime
import itertools
import json
import resource
import time
import zipfile
from common import ROOT, OUTPUT, ARCHIVE, sha, write, read, check_frozen, checkpoint
from extract import source_rows, extract_row, case_id, header
from clinical_rules import periodontal, caries, control_periodontal, control_caries

def run(mode):
    start = time.perf_counter()
    cpu = time.process_time()
    check_frozen('PREREG_R1_REPORT_LINKS.json')
    if mode == 'R2':
        check_frozen('PREREG_R2_LOCALIZED_OBSERVATIONS.json')
    (rows, volumes, csv_hash) = source_rows()
    snapshots = []
    observations = []
    rejections = collections.Counter()
    unique_volumes = {}
    record_cases = set()
    paired_cases = set()
    with zipfile.ZipFile(ARCHIVE) as z:
        for (n, row) in enumerate(rows, 2):
            original = row['Filename']
            pid = case_id(original)
            record_cases.add(pid)
            member = f'MMDental/{original}/{original}.nii.gz'
            paired = member in volumes
            if not paired:
                rejections['NO_CBCT_CASE_MATCH'] += 1
                continue
            if pid not in unique_volumes:
                try:
                    h = header(z, member)
                except Exception as e:
                    rejections['BAD_HEADER_' + type(e).__name__] += 1
                    continue
                info = volumes[member]
                unique_volumes[pid] = {'member': member, 'uncompressed_bytes': info.file_size, 'zip_crc32': info.CRC, **h, 'header_not_full_volume_integrity': True, 'visit_contemporaneity': 'UNKNOWN'}
            paired_cases.add(pid)
            obs = extract_row(row, n, csv_hash, mode)
            observations.extend(obs)
            snapshots.append({'patient_id': pid, 'snapshot_row': n, 'volume': unique_volumes[pid], 'observation_ids': [o['observation_id'] for o in obs], 'BD10': periodontal(obs), 'BD08': caries(obs), 'no_clinical_certification': True, 'acquisition_date': 'UNKNOWN'})
    pred = {'schema': 'X57-frozen-predictions-v1', 'round': mode, 'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'claim_type': 'information_link', 'csv_sha256': csv_hash, 'code_hashes': {p.name: sha(p.read_bytes()) for p in (ROOT / 'code').glob('*.py')}, 'predictions': [{k: s[k] for k in ['patient_id', 'snapshot_row', 'BD10', 'BD08']} for s in snapshots], 'diagnosis_used_for_predictions': False, 'endpoint': 'reported advanced destruction / reported caries-pulp penetration; no clinical stage or treatment truth'}
    write(f'FROZEN_PREDICTIONS_{mode}.json', pred)
    (OUTPUT / f'FROZEN_PREDICTIONS_{mode}.sha256').write_text(sha((OUTPUT / f'FROZEN_PREDICTIONS_{mode}.json').read_bytes()) + '\n')
    concordance = collections.defaultdict(list)
    for s in snapshots:
        diag = rows[s['snapshot_row'] - 2]['Diagnosis']
        import re
        pdiag = bool(re.search('(?<!apical )\\b(?:chronic\\s+)?periodontitis\\b|K05\\.3', diag, re.I))
        clean = re.sub('(?:chronic\\s+)?apical\\s+periodontitis', '', diag, flags=re.I)
        pdiag = bool(re.search('\\bperiodontitis\\b|K05\\.3', clean, re.I))
        endo = bool(re.search('pulpitis|K04\\.00|pulp\\s+necros', diag, re.I))
        for (chain, flag, label) in [('BD10', s['BD10']['advanced_threshold_definitely_reported'], pdiag), ('BD08', s['BD08']['pulp_penetration_definitely_reported'], endo)]:
            if flag:
                concordance[chain].append({'patient_id': s['patient_id'], 'snapshot_row': s['snapshot_row'], 'diagnosis_category_present': label, 'diagnosis_sha256': sha(diag.encode()), 'reference_column': 'Diagnosis', 'interpretation': 'concordance only; label absence not disease-negative'})
    accepted = [o for o in observations if o['accepted']]
    sets = {c: set((o['patient_id'] for o in accepted if o['chain'] == c)) for c in ['BD10', 'BD08']}
    counts = {c: len(sets[c]) for c in sets}
    controls = []
    for s in snapshots:
        obs = [o for o in observations if o['snapshot_row'] == s['snapshot_row']]
        p = control_periodontal(obs)
        e = control_caries(obs)
        control_ok = p == s['BD10']['conditional_stages'] and all((x['classes'] == y['conditional_lesion_classes'] for (x, y) in zip(e, s['BD08']['lesions'])))
        controls.append({'snapshot_row': s['snapshot_row'], 'ok': control_ok})
    n_joint = len(sets['BD10'] & sets['BD08'])
    joint_localized = []
    for s in snapshots:
        local = [o for o in accepted if o['snapshot_row'] == s['snapshot_row'] and o['teeth']]
        if all((any((o['chain'] == chain for o in local)) for chain in ['BD10', 'BD08'])):
            joint_localized.append(s['snapshot_row'])
    summary = {'round': mode, 'claim_type': 'information_link', 'resolution': 'POPULATION', 'timescale': 'SIMULTANEOUS', 'record_rows': len(rows), 'record_cases': len(record_cases), 'archive_volumes': len(volumes), 'paired_snapshots': len(snapshots), 'paired_cases': len(paired_cases), 'rejected_rows': dict(rejections), 'rejected_row_fraction': sum(rejections.values()) / len(rows), 'observations_candidates': len(observations), 'accepted_observations': len(accepted), 'rejected_observations': len(observations) - len(accepted), 'observation_rejections': dict(collections.Counter((o['rejection'] for o in observations if not o['accepted']))), 'by_chain': dict(collections.Counter((o['chain'] for o in accepted))), 'by_quantity': dict(collections.Counter((o['quantity'] for o in accepted))), 'paired_cases_by_chain': counts, 'joint_case_count': n_joint, 'joint_localized_snapshot_rows': joint_localized, 'joint_snapshot_rows': [s['snapshot_row'] for s in snapshots if s['BD10']['supports'] and s['BD08']['lesions']], 'report_link_gate': 'PASS' if min(counts.values()) >= 10 and n_joint >= 1 else 'FAIL', 'R2_joint_information_gate': 'PASS' if mode == 'R2' and joint_localized else 'FAIL' if mode == 'R2' else 'NOT_APPLICABLE', 'clinical_answer_gate': 'UNKNOWN_MISSING_MEASUREMENTS_AND_INDEPENDENT_ADJUDICATION', 'complete_clinical_answers': 0, 'same_information_control_disagreements': sum((not x['ok'] for x in controls)), 'concordance': {c: {'flagged_snapshots': len(rs), 'same_record_diagnosis_category': sum((x['diagnosis_category_present'] for x in rs)), 'meaning': 'Not stage/lesion calibration; duplicated visits not independent'} for (c, rs) in concordance.items()}, 'full_cost': {'wall_s': time.perf_counter() - start, 'cpu_s': time.process_time() - cpu, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'fit_s': 0, 'questions': 0, 'gpu_s': 0, 'original_acquisition': 'UNKNOWN', 'clinical_validation': 'NOT_PERFORMED', 'agent_tokens': 'UNKNOWN'}, 'csv_sha256': csv_hash, 'archive_sha256': 'UNKNOWN_NOT_FULLY_REREAD', 'external_referent': read('PREREG_R1_REPORT_LINKS.json')['external_referent'], 'strict_patient_anatomy_model': 'CBCT header and exact case locator only; no unlabeled anatomy inferred'}
    write(f'raw/OBSERVATIONS_{mode}.json', observations)
    write(f'raw/PATIENT_MODELS_{mode}.json', snapshots)
    write(f'raw/DIAGNOSIS_CONCORDANCE_{mode}.json', dict(concordance))
    write(f'raw/CONTROLS_{mode}.json', controls)
    write(f'raw/RESULTS_{mode}.json', summary)
    ranked = sorted(observations, key=lambda o: sha(('producer-audit-' + o['observation_id']).encode()))[:30]
    write(f'raw/REVIEW_SAMPLE_{mode}.json', ranked)
    checkpoint(mode + '_COMPUTED', summary['report_link_gate'] + ' observation-link; clinical answer UNKNOWN', 'Review source scope; freeze changed construction without changing R1 gates')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--round', choices=['R1', 'R2'], default='R1')
    run(p.parse_args().round)
