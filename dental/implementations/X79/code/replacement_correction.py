"""R5: changed observation representation for absent native + present prosthetic site.
Preserves R1/R2/R3 and their parser. Only a new source correction artifact is emitted.
"""
from dental_release.paths import expand as _release_expand
from common import *
from data import reports, load_geometry
import copy, re

def main():
    freeze('PREREG_R5_NATIVE_VS_SITE.json', dict(frozen_utc=now(), claim_type='information_link', capability='Keep native anatomical absence and presence of a prosthetic replacement as simultaneous compatible per-tooth-site observations', obstacle='R1 grammar propagates absence of native46 onto a reduced-size prosthetic element present at46. A later full negative-assertion review found the source error.', changed_operation='Preserve the frozen grammar and every original result. New representation distinguishes native tooth from tooth-site replacement; source clauses with named absent tooth and prosthetic element explicitly present receive positive prosthetic-site material state.', consumer='R2 native-only material applicability; R3 absence remains replacement-excluded', metric=dict(original_wrong_polarity_must_be_rejected=True, corrected_observation_source_identity_exact=True, no_refit_no_threshold_change=True), strongest_equal_information_control='Read literal explicit present/absent predicates with separate native/prosthetic objects; same source bytes', external_referent=dict(kind='published_dataset', locator=str(ZIP) + _release_expand('::@DENTAL_CASE_ID@/reports_intraoral-photo_en/18_3735.txt + @DENTAL_CASE_ID@/reports_ios_en/18_1842.txt'), compared_quantity='Native46 is absent while reduced-size prosthetic replacement is present; material presence is positive at tooth site46', refutes_us=True), resolution='PER_TOOTH', timescale='SIMULTANEOUS', full_cost=dict(preparation='Existing text/geometry reused with hashes', fit='none', discovery='full emitted named-negative review', validation='source predicate guard + injected polarity', queries='explicit replacement update', fallback='UNKNOWN exact material, boundary, anatomical replacement mapping', threads=2, max_intermediate_bytes=3000000000)))
    freeze('FROZEN_PREDICTIONS_R5.json', dict(frozen_utc=now(), phase='Semantic source prediction before correction join; no new physical measurement', prediction='Native tooth46 absence is compatible with positive prosthetic replacement state at site46; a negative restorative-state parser output must be rejected', prereg_sha256=sha(ROOT / 'PREREG_R5_NATIVE_VS_SITE.json')))
    st = time.perf_counter()
    cpu = time.process_time()
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    ids = m['train'] + m['calibration'] + m['test']
    corrections = []
    for modality in ['ios', 'intraoral-photo']:
        rr = reports(ids, modality)
        for (c, reps) in rr.items():
            for r in reps:
                for a in r['parsed']['assertions']:
                    if not a['polarity'] and a['teeth'] and re.search('absence of tooth\\s+\\d+.*prosthetic element is present', a['clause'], re.I):
                        corrections.append(dict(case_id=c, fdi=a['teeth'][0], source_member=r['member'], source_sha256=r['sha256'], clause=a['clause'], original_restoration_polarity=0, corrected_prosthetic_site_polarity=1, native_tooth_state='ABSENT_REPORTED', site_replacement_state='PRESENT_REPORTED', resolution='PER_TOOTH', timescale='SIMULTANEOUS'))
    ports = read('raw/MATERIAL_APPLICABILITY_PORT.json')
    corrected = copy.deepcopy(ports)
    photo = {(r['case_id'], r['fdi']): r for r in corrections if '/reports_intraoral-photo_en/' in r['source_member']}
    for p in corrected:
        r = photo.get((p['case_id'], p['fdi']))
        if r:
            p.update(restoration_state='PROSTHETIC_SITE_REPLACEMENT_REPORTED', reported_class='prosthetic_replacement', natural_only_closure='REFUTED_BY_REPORT', native_tooth_state='ABSENT_REPORTED', source={k: r[k] for k in ['source_member', 'source_sha256', 'clause']})
    write('raw/R5_SOURCE_CORRECTIONS.json', corrections)
    write('raw/MATERIAL_APPLICABILITY_PORT_R5.json', corrected)
    from reviewed_source_guard import source_guard
    checks = dict(wrong_FDI_rejected=not source_guard(dict(corrections[0], fdi=11), 1), wrong_hash_rejected=not source_guard(dict(corrections[0], source_sha256='0' * 64), 1), original_injected_wrong_polarity_rejected=bool(corrections) and all((not source_guard(r, 0) for r in corrections)), corrected_polarity_accepted=bool(corrections) and all((source_guard(r, 1) for r in corrections)), wrong_source_case_rejected=bool(corrections) and (not source_guard(dict(corrections[0], case_id=_release_expand('@DENTAL_CASE_ID@')), 1)))
    write('raw/R5_RESULTS.json', dict(claim_type='information_link', outcome='PASS_NATIVE_VS_SITE_SOURCE_CORRECTION' if all(checks.values()) else 'FAIL', external_referent=read('PREREG_R5_NATIVE_VS_SITE.json')['external_referent'], corrections=corrections, original_source_scope_gate='FAIL_ON_F4970: do not erase or relabel original R1/R2 artifacts', corrected_positive_ports=sum((r['natural_only_closure'] == 'REFUTED_BY_REPORT' for r in corrected)), new_positive_patient_count=len({r['case_id'] for r in corrected if r['natural_only_closure'] == 'REFUTED_BY_REPORT'}), injections=checks, cost=cost(st, cpu)))
    (ROOT / 'HANDOFF_R5.md').write_text(_release_expand('R5 source refutes original parser on @DENTAL_CASE_ID@: native46 absent, prosthetic replacement present. New native-vs-site state corrects the material port, preserving original metrics/labels/thresholds. Next: named-rater tooth/site labels and physical material/interface measurements.\n'))
    state('R5_DECIDED', _release_expand('Source correction; original source scope FAIL on @DENTAL_CASE_ID@'), 'Package original and corrected ports, complete-negative measurement panel, one-command replay')
    print('R5 source correction', len(corrections), 'clauses', checks, flush=True)
if __name__ == '__main__':
    main()
