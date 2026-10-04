"""R2: contralateral information, tested with the target pair excluded."""
import json
from pathlib import Path
import time
import zipfile
import numpy as np
from measure import HERE, ROOT, DATA, POSTERIOR, tf2_io, unit, sha, dump, check_frozen, Local, trace, restoration_near, scalar_decide, write_csv

def pose(j, counts, fdi):
    homolog = fdi + 10 if fdi // 10 == 1 else fdi - 10
    tooth_no = fdi % 10
    pairs = []
    for k in range(1, 9):
        if k == tooth_no:
            continue
        if min(int(counts.get(str(10 + k), 0)), int(counts.get(str(20 + k), 0))) >= 500:
            cr = np.asarray(j['teeth'][str(10 + k)]['centroid'])
            cl = np.asarray(j['teeth'][str(20 + k)]['centroid'])
            pairs.append((cr, cl))
    if len(pairs) < 3 or int(counts.get(str(homolog), 0)) < 500:
        return None
    origin = np.mean([(l + r) / 2 for (r, l) in pairs], axis=0)
    n = unit(np.mean([l - r for (r, l) in pairs], axis=0))
    p = np.asarray(j['teeth'][str(homolog)]['crest_entry'])
    a = -unit(j['teeth'][str(homolog)]['axis'])
    anchor = p - 2 * np.dot(p - origin, n) * n
    axis = unit(a - 2 * np.dot(a, n) * n)
    normal = unit(n - np.dot(n, axis) * axis)
    tangent = unit(np.cross(normal, axis))
    return (anchor, axis, normal, tangent, len(pairs))

def run():
    start = time.monotonic()
    check_frozen('PREREG_R2')
    check_frozen('FROZEN_PREDICTIONS_R2')
    cases = [json.loads(s) for s in (ROOT / 'results/BONE_HU/tf2_census.jsonl').read_text().splitlines()]
    (wanted, validation) = ([], [])
    for c in cases:
        q = c['counts']
        if sum((int(q.get(str(i), 0)) for i in (5, 6))) < 500:
            continue
        jp = ROOT / f"results/NV1_canals/per_case/{c['case']}.json"
        j = json.loads(jp.read_text())
        for t in POSTERIOR:
            try:
                f = pose(j, q, t)
            except KeyError:
                continue
            if f is None:
                continue
            (anchor, axis, normal, tangent, npairs) = f
            rec = {'case': c['case'], 'fdi': t, 'stratum': 'PM' if t % 10 < 6 else 'M', 'anchor': anchor.tolist(), 'axis': axis.tolist(), 'normal': normal.tolist(), 'tangent': tangent.tolist(), 'n_other_pairs': npairs, 'frame_sha256': sha(jp), 'clinical_decision': 'UNKNOWN', 'location_truth': 'UNKNOWN_EDENTULOUS_SITE'}
            if str(t) not in q:
                wanted.append(rec)
            elif int(q[str(t)]) >= 500:
                ptrue = np.asarray(j['teeth'][str(t)]['crest_entry'])
                atrue = -unit(j['teeth'][str(t)]['axis'])
                d = anchor - ptrue
                rec.update(crest_error_3d_mm=float(np.linalg.norm(d)), transverse_error_mm=float(np.linalg.norm(d - np.dot(d, atrue) * atrue)), axis_error_deg=float(np.degrees(np.arccos(np.clip(np.dot(axis, atrue), -1, 1)))), true_crest=ptrue.tolist(), true_axis=atrue.tolist())
                validation.append(rec)
    dump(HERE / 'MIRROR_SELECTED_SITES.json', wanted)
    dump(HERE / 'MIRROR_HELDOUT_PREDICTIONS.json', validation)
    vr = {'n_heldout': len(validation), 'p95_transverse_error_mm': float(np.percentile([r['transverse_error_mm'] for r in validation], 95)) if validation else None, 'median_axis_error_deg': float(np.median([r['axis_error_deg'] for r in validation])) if validation else None, 'median_crest_error_3d_mm': float(np.median([r['crest_error_3d_mm'] for r in validation])) if validation else None}
    vr['location_gate'] = 'UNKNOWN_MINIMUM_N' if len(validation) < 30 else 'PASS' if vr['p95_transverse_error_mm'] <= 2 else 'FAIL'
    vr['axis_gate'] = 'UNKNOWN_MINIMUM_N' if len(validation) < 30 else 'PASS' if vr['median_axis_error_deg'] <= 10 else 'FAIL'
    (rows, manifest) = ([], [])
    with zipfile.ZipFile(tf2_io.ZIP) as z:
        for case in sorted(set((r['case'] for r in wanted))):
            name = f'{tf2_io.ROOT}/labelsTr/{case}.mha'
            raw = z.read(name)
            (lab, sp, header) = tf2_io.read_mha_bytes(raw)
            manifest.append({'member': name, 'bytes': len(raw), 'sha256': __import__('hashlib').sha256(raw).hexdigest()})
            for r in [r for r in wanted if r['case'] == case]:
                (anchor, axis, normal, tangent) = (np.asarray(r[k]) for k in ('anchor', 'axis', 'normal', 'tangent'))
                local = Local(lab, sp, anchor, axis)
                r.update(trace(local, anchor, axis, normal))
                (near, n) = restoration_near(local, anchor, axis)
                r['restoration_voxels_near'] = n
                if near:
                    r.update(valid=False, reason='restoration_at_proposed_gap')
                h = r.get('bone_height_mm')
                w = r.get('widths_mm', {}).get('1')
                r['rule5_flag'] = 'lift_trigger' if h is not None and h < 5 else 'no_lift_trigger' if h is not None else 'unknown'
                r['scalar_nominal'] = scalar_decide(r, half=0)
                r['scalar_interval'] = scalar_decide(r)
                r['height_interval_mm'] = [max(0, h - 0.6), h + 0.6] if h is not None else None
                r['width1_interval_mm'] = [max(0, w - 0.6), w + 0.6] if w is not None else None
                r['shortest_scalar_catalog_length_mm'] = 6 if r['scalar_interval'] == 'short_geometry' else None
                reread = trace(local, anchor, axis, normal, order=0)
                errors = {k: abs(r[k] - reread[k]) for k in ('bone_height_mm', 'crest_to_sinus_label_mm') if k in r and k in reread}
                for d in ('1', '3', '7'):
                    if r.get('widths_mm', {}).get(d) is not None and reread.get('widths_mm', {}).get(d) is not None:
                        errors['width' + d] = abs(r['widths_mm'][d] - reread['widths_mm'][d])
                r['independent_voxel_read'] = reread
                r['cross_read_errors_mm'] = errors
                r['cross_read_pass'] = bool(errors) and all((e <= 0.15 + 1e-08 for e in errors.values()))
                path = DATA / f"{case}_{r['fdi']}_mirror.npz"
                local.save(path, anchor, axis, normal, tangent)
                manifest.append({'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size})
                rows.append(r)
                print(case, r['fdi'], r['reason'], r.get('bone_height_mm'), r['scalar_interval'], flush=True)
            del lab, raw
    with open(HERE / 'RAW_R2.jsonl', 'w') as f:
        for r in rows:
            f.write(json.dumps(r) + '\n')
    write_csv(rows, HERE / 'MEASUREMENTS_R2.csv')
    dump(HERE / 'DATA_MANIFEST_R2.json', manifest)
    r1 = [json.loads(s) for s in (HERE / 'RAW_R1.jsonl').read_text().splitlines()]
    r1keys = {(r['case'], r['fdi']) for r in r1}
    new = [r for r in rows if (r['case'], r['fdi']) not in r1keys]
    valid = [r for r in new if r['valid']]
    summary = dict(vr, n_selected=len(rows), n_new_selected=len(new), n_new_evaluable=len(valid), coverage_gate='PASS' if len(valid) >= 8 else 'FAIL', counts_new={k: sum((r['scalar_interval'] == k for r in new)) for k in ('short_geometry', 'height_augmentation_geometry', 'uncertain')}, elapsed_s=time.monotonic() - start)
    combined = [dict(r, construction='R1') for r in r1] + [dict(r, construction='R2') for r in new]
    with open(HERE / 'COMBINED_SITES.jsonl', 'w') as f:
        for r in combined:
            f.write(json.dumps(r) + '\n')
    dump(HERE / 'SUMMARY_R2.json', summary)
    dump(HERE / 'CURRENT_WORK_STATE.json', {'lane': 'X20-short-implant-sinus', 'milestone': 'R2_completed', 'latest_gate': summary, 'next_operation': 'Finite-body feasibility and local location/edge information required to resolve each site'})
    print(json.dumps(summary, indent=2))
if __name__ == '__main__':
    run()
