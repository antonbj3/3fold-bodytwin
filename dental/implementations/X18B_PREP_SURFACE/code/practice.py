"""Re-evaluate historical generic anatomy on its unchanged X18 reference grid."""
from dental_release.paths import expand as _release_expand
from common import *

def run():
    path = H / 'PREREG_R6.json'
    if not path.exists():
        dump(path, dict(round='R6', claim_type='information_link', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), capability='Quantify preoperative source-surface information against the existing generic-donor practice proxy on a common reference grid', obstacle='R1 preserves source exteriors but its new0.15mm sampling cannot be numerically compared to historical X18 summary grids', changed_operation='Query exported source carrier and historical frozen generic roof at the same original X18 points, with the same recorded antagonist ceiling', consumer='Measured-anatomy information comparison', cohort='Cases127,122,188, mandibular34–37 where original X18 frozen prediction exists; no new fit', metrics=dict(contact_IoU_gain_min=0.15, contact_area_error_ratio_max=0.8, surface_centroid_error_max_mm=1e-05, patch_count_error_max=0), strongest_equally_informed_control='Untouched original source roof at identical coordinates', practice='Original X18 generic donor without registered antagonist; an explicit software proxy, no commercial CAD evaluation', falsifiers=['Different grids called same measurement', 'Injected exterior shift passes source check', 'Preoperative source used but result called withheld anatomy prediction'], external_referent=dict(kind='published_dataset', locator=_release_expand('@DENTAL_DATA_ROOT@/geometry/Bits2Bites/Bits2Bites_v01.zip'), compared_quantity='Original scan-derived proximity masks on original X18 sampling support', refutes_us=True), prospective_validity='This is a preregistered descriptive historical replay after R1 results were seen; no new held-out prediction claim', resolution='PER_SURFACE_REGION', full_cost=dict(preparation='Read existing frozen priors and source carriers', fit='none; inherited donor fit charged as existing cost', discovery='R1 known; historical grid mismatch disclosed', validation='Same-grid independent source queries', questions='12 source-envelope queries', fallback='UNKNOWN if prior or source absent')))
    pr = json.loads(path.read_text())
    old = json.loads((X18 / 'raw/PREDICTIONS_R1.json').read_text())
    rows = []
    files = []
    cur = None
    for case in [127, 122, 188]:
        (data, _) = pair(case)
        for fdi in [34, 35, 36, 37]:
            r = next((r for r in old if r['case'] == case and r['fdi'] == fdi and (r['status'] == 'DESIGNED')), None)
            if r is None:
                rows.append(dict(case=case, fdi=fdi, status='REJECTED', reason='Historical original prediction absent'))
                continue
            files.append(artifact(r['file']))
            z = np.load(r['file'])
            xy = z['xy']
            carrier = np.load(D / f'{case}_{fdi}_carrier.npz')
            external = carrier['external_vertices'][carrier['external_faces']]
            original = data['lower']['tri'][data['lower']['owner'] == fdi]
            (source_z, _) = query_height(original, xy, False)
            (retained_z, _) = query_height(external, xy, False)
            truth = contact(z['ceiling'] - source_z, source_z, xy, z['index'], z['weights'], 0.1, 0.1)
            candidate = contact(z['ceiling'] - retained_z, retained_z, xy, z['index'], z['weights'], 0.1, 0.1)
            practice = contact(z['ceiling'] - z['z_practice'], z['z_practice'], xy, z['index'], z['weights'], 0.1, 0.1)
            ce = score(candidate, truth, xy)
            pe = score(practice, truth, xy)
            truth.pop('mask')
            candidate.pop('mask')
            practice.pop('mask')
            rows.append(dict(case=case, fdi=fdi, status='SCORED', resolution='PER_SURFACE_REGION', original=truth, retained_source=candidate, generic_practice=practice, candidate_errors=ce, practice_errors=pe, IoU_gain=ce['IoU'] - pe['IoU'], prior_file=r['file'], scope='Target surface was preoperative measured input, not held out from construction'))
    p = H / 'raw/PREDICTIONS_R6.json'
    dump(p, rows)
    files.append(artifact(p))
    dump(H / 'FROZEN_PREDICTIONS_R6.json', dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(path), files=files, code_sha256={str(H / 'code/practice.py'): sha(H / 'code/practice.py')}, future_physical_measurement_status='NOT_RUN', historical_reference_already_known=True))
    rr = [r for r in rows if r['status'] == 'SCORED']
    gain = float(np.median([r['IoU_gain'] for r in rr]))
    ratio = float(np.median([r['candidate_errors']['absolute_area_error_mm2'] for r in rr]) / max(0.1, np.median([r['practice_errors']['absolute_area_error_mm2'] for r in rr])))
    gates = dict(IoU=gain >= pr['metrics']['contact_IoU_gain_min'], area=ratio <= pr['metrics']['contact_area_error_ratio_max'], count=all((r['candidate_errors']['absolute_patch_count_error'] == 0 for r in rr)), centroid=all((r['candidate_errors']['centroid_error_mm'] <= pr['metrics']['surface_centroid_error_max_mm'] for r in rr)))
    out = dict(round='R6', claim_type='information_link', external_referent=pr['external_referent'], decision='COMMON_GRID_INFORMATION_PASS' if all(gates.values()) else 'COMMON_GRID_INFORMATION_FAIL', gates=gates, attempted=len(rows), retained=len(rr), rejection_fraction=(len(rows) - len(rr)) / len(rows), median_IoU_gain=gain, area_error_ratio=ratio, median_practice_IoU=float(np.median([r['practice_errors']['IoU'] for r in rr])), median_candidate_IoU=float(np.median([r['candidate_errors']['IoU'] for r in rr])), median_practice_count_error=float(np.median([r['practice_errors']['absolute_patch_count_error'] for r in rr])), median_practice_centroid_error_mm=float(np.median([r['practice_errors']['centroid_error_mm'] for r in rr])), rows=rows, scope=pr['prospective_validity'])
    dump(H / 'rounds/R6.json', out)
    (H / 'HANDOFF_R6.md').write_text(f"R6 {out['decision']}: historical generic-donor proxy versus exact measured target exterior on same original sample support; n={len(rr)}, medianIoUgain={gain}. No independent clinical-contact or withheld shape prediction. All previous R1/R2/R4 failures retained.\n")
    print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
if __name__ == '__main__':
    run()
