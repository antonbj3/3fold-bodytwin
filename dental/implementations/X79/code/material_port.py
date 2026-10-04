from dental_release.paths import expand as _release_expand
from common import *
from data import load_geometry, treatment_rows, reports
import collections, numpy as np, zipfile

def main():
    st = time.perf_counter()
    cpu = time.process_time()
    g = load_geometry()
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    freeze('FROZEN_PREDICTIONS_R2.json', dict(frozen_utc=now(), claim_type='information_link', phase='No-report consumer recipe frozen before new observation join; source labels were evaluated in R1', prediction='Natural-only closure would be unrefuted at every represented X11 tooth; every source-linked treatment will contradict it', domain_sha256=sha(ROOT / 'raw/TOOTH_FEATURES.json'), parser_sha256=sha(ROOT / 'code/report_parser.py'), prereg_sha256=sha(ROOT / 'PREREG_R2_MATERIAL_INFORMATION_PORT.json')))
    observations = {}
    source_rows = []
    excluded = []
    for split in ['train', 'calibration', 'test']:
        rr = read('raw/R1_' + ('TEST' if split == 'test' else split) + '_OBSERVATIONS.json')
        observations.update(rr)
        (rows, ex) = treatment_rows(rr, g)
        source_rows.extend((dict(**r, split=split) for r in rows))
        excluded.extend(ex)
    bykey = {(r['case_id'], r['fdi']): r for r in source_rows}
    ports = []
    manifest_by_case = {r['case_id']: r for r in read('raw/GEOMETRY_INPUT_MANIFEST.json')}
    for (c, teeth) in sorted(g.items()):
        for (f, geom) in sorted(teeth.items()):
            if geom['label_count'] < 100:
                continue
            r = bykey.get((c, f))
            positive = r is not None and r['y'] == 1
            ports.append(dict(case_id=c, fdi=f, restoration_state='RESTORATION_REPORTED' if positive else 'NO_RESTORATION_REPORTED' if r else 'UNKNOWN', reported_class=r['kind'] if positive else None, material_identity='UNKNOWN', restoration_boundary_volume='UNKNOWN', natural_only_closure='REFUTED_BY_REPORT' if positive else 'NOT_VALIDATED', geometry_resolution='PER_TOOTH', report_resolution='PER_TOOTH', edge_resolution='PER_TOOTH', timescale='SIMULTANEOUS', time_alignment='Same case; paired acquisition time UNKNOWN', source={k: r[k] for k in ['source_member', 'source_sha256', 'original_it_member', 'original_it_sha256', 'clause', 'scope']} if r else None, geometry=dict(source_input_sha256=manifest_by_case[c]['source_STL_hashes'], label_count=geom['label_count'])))
    write('raw/MATERIAL_APPLICABILITY_PORT.json', ports)
    positives = {(r['case_id'], r['fdi']) for r in source_rows if r['y'] == 1}
    control = [(r['case_id'], r['fdi']) in positives for r in ports]
    candidate = [r['natural_only_closure'] == 'REFUTED_BY_REPORT' for r in ports]
    parity = sum((a != b for (a, b) in zip(candidate, control)))

    def admissible(r):
        source = observations.get(r['case_id'], [])
        if not source:
            return False
        s = source[0]
        return r['source_member'] == s['member'] and r['source_sha256'] == s['sha256'] and any((a['polarity'] == r['y'] and (r['fdi'] in a['teeth'] or (a['scope'] == 'ALL_REPORTED_TEETH' and r['y'] == 0)) for a in s['parsed']['assertions']))
    good = next((r for r in source_rows if r['y'] == 1))
    wrongf = next((f for f in range(11, 49) if f not in {a for r in observations[good['case_id']][0]['parsed']['assertions'] for a in r['teeth']}))
    inj = dict(correct_rows_accepted=all((admissible(r) for r in source_rows)), wrong_FDI_rejected=not admissible(dict(good, fdi=wrongf)), wrong_case_rejected=not admissible(dict(good, case_id=_release_expand('@DENTAL_CASE_ID@'))), wrong_hash_rejected=not admissible(dict(good, source_sha256='0' * 64)), wrong_polarity_rejected=not admissible(dict(good, y=0)), same_information_control_injected_state_rejected=not all((a == b for (a, b) in zip([not candidate[0]] + candidate[1:], control))))
    restored = [r for r in source_rows if r['y'] == 1]
    test = [r for r in restored if r['split'] == 'test']
    port_keys = {(p['case_id'], p['fdi']) for p in ports}
    qa = read('raw/PARSER_GOLD_VALIDATION.json')
    gate = len(test) >= 50 and len({r['case_id'] for r in test}) >= 20 and (qa['precision'] >= 0.95) and (parity == 0) and all(inj.values())
    result = dict(claim_type=['information_link', 'capability'], outcome='PASS_SOURCE_LINKED_APPLICABILITY' if gate else 'FAIL', external_referent=read('PREREG_R2_MATERIAL_INFORMATION_PORT.json')['external_referent'], restored_teeth=len(restored), restored_patients=len({r['case_id'] for r in restored}), test_restored_teeth=len(test), test_restored_patients=len({r['case_id'] for r in test}), exported_ports=len(ports), exported_positive_ports=sum(candidate), reported_restored_without_100_vertex_support=sum(((r['case_id'], r['fdi']) not in port_keys for r in restored)), natural_material_assumptions_now_refuted=len(restored), matched_control_parity_mismatches=parity, injections=inj, clinical_material_calibration='UNKNOWN; report cannot identify material constitutive law or restoration volume', negative_source='Explicit no-treatment means clinician did not observe a treatment; it does not validate enamel-only state', cost=cost(st, cpu))
    write('raw/R2_RESULTS.json', result)
    (ROOT / 'HANDOFF_R2.md').write_text(f"R2 {result['outcome']}: {len(restored)} named treatment teeth/{result['restored_patients']} patients, {len(test)} teeth/{result['test_restored_patients']} test patients. New per-tooth material applicability input; mechanics is not calibrated. Exact source membership control agrees. Next construction: quantify presence detector and determine whether spacing/wear sources identify full binary endpoints.\n")
    state('R2_DECIDED', result['outcome'], 'R3 measure explicit absence vs named presence and refuse unsupported spacing specificity')
    print(json.dumps(result), flush=True)
if __name__ == '__main__':
    main()
