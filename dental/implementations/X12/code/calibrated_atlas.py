"""R8: physical CT identity inherited, side semantics calibrated on two anchors."""
import csv, hashlib, json, os, resource, sys, time, zipfile
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from archive_io import index, nifti, ARCHIVE
from pair_atlas_r6 import ROOT, DATA, FDI, transform, tooth_measure, write_csv
from tf2_io import ZIP, ROOT as TFROOT, read_mha_bytes

def semantic_map(pulp, tooth):
    coords = np.argwhere(np.isin(pulp, list(FDI)))
    v = pulp[tuple(coords.T)]
    tv = tooth[tuple(coords.T)]
    present = sorted((int(q) for q in np.unique(v)))
    anchors = present[:2]
    swap = np.where(v < 40, v + 10, v - 10)
    profiles = []
    for (name, canonical) in [('identity', v), ('left_right_swap', swap)]:
        inside = tv == canonical
        fractions = {str(fdi): float(np.mean(inside[v == fdi])) for fdi in present}
        calibration_ok = len(anchors) == 2 and all((fractions[str(fdi)] >= 0.95 for fdi in anchors))
        heldout = ~np.isin(v, anchors)
        profiles.append({'map': name, 'calibration_pass': calibration_ok, 'fractions': fractions, 'heldout_containment': float(np.mean(inside[heldout])) if heldout.any() else None, 'overall_containment': float(np.mean(inside)) if len(v) else 0})
    accepted = [p for p in profiles if p['calibration_pass']]
    if len(accepted) != 1:
        return (None, {'anchors': anchors, 'profiles': profiles, 'reason': 'NONUNIQUE_OR_FAILED_ANCHORS'})
    winner = accepted[0]
    if winner['heldout_containment'] is None or winner['heldout_containment'] < 0.95:
        return (None, {'anchors': anchors, 'profiles': profiles, 'reason': 'HELDOUT_CONTAINMENT_FAIL'})
    return (winner['map'], {'anchors': anchors, 'profiles': profiles, 'selected_map': winner['map'], 'heldout_containment': winner['heldout_containment'], 'overall_containment': winner['overall_containment']})

def main():
    started = time.time()
    cpu = time.process_time()
    rows = index()
    lookup = {r['name']: r for r in rows}
    identity = json.loads((ROOT / 'raw/R6_pairs.json').read_text())
    if not (ROOT / 'raw/R6_cost.json').exists():
        raise RuntimeError('R6 must finish before R8 consumes its complete evidence')
    source_stat = {str(ARCHIVE): {'size': ARCHIVE.stat().st_size, 'mtime_ns': ARCHIVE.stat().st_mtime_ns}, str(ZIP): {'size': Path(ZIP).stat().st_size, 'mtime_ns': Path(ZIP).stat().st_mtime_ns}}
    original_source = json.loads((ROOT / 'raw/R6_source_stat.json').read_text()) if (ROOT / 'raw/R6_source_stat.json').exists() else None
    if original_source is not None and source_stat != original_source:
        raise RuntimeError('Source archive stat changed since full CT equality; rerun R6')
    out = []
    measure = []
    excluded = []
    controls = []
    manifest = []
    saved = 0
    with zipfile.ZipFile(ZIP) as z:
        for record in identity:
            pid = record['case']
            if not record.get('paired') or len(record['image_transforms']) != 1 or (not record.get('label_grid_consistent', False)):
                out.append({'case': pid, 'accepted': False, 'reason': 'NO_UNIQUE_EXACT_CT_IDENTITY'})
                continue
            image_pose = record['image_transforms'][0]
            pname = 'Pulpy3D/' + pid + '/gt_instance.nii.gz'
            (pimg, psha) = nifti(lookup[pname])
            pulp = np.asanyarray(pimg.dataobj)
            pulp = transform(pulp, image_pose[0], image_pose[1])
            lname = f"{TFROOT}/labelsTr/{record['tf2_case']}.mha"
            raw = z.read(lname)
            (tooth, sp, h) = read_mha_bytes(raw)
            if psha != record['pulpy_pulp_sha256'] or hashlib.sha256(raw).hexdigest() != record['TF2_tooth_sha256']:
                raise RuntimeError('Source changed after R6 pairing')
            (mapping, cal) = semantic_map(pulp, tooth)
            rec = {'case': pid, 'paired': True, 'accepted': mapping is not None, 'calibration': cal, 'CT_identity_evidence': 'raw/R6_pairs.json#' + pid, 'pulpy_CT_sha256': record['pulpy_CT_sha256'], 'tf2_CT_sha256': record['tf2_CT_sha256'], 'pulp_sha256': psha, 'tooth_sha256': record['TF2_tooth_sha256']}
            if mapping:
                if mapping == 'left_right_swap':
                    lut = np.arange(65536, dtype=np.uint16)
                    for fdi in FDI:
                        lut[fdi] = fdi + 10 if fdi < 40 else fdi - 10
                    mapped = lut[pulp]
                else:
                    mapped = pulp
                boxes = ndi.find_objects(tooth, max_label=48)
                for fdi in FDI:
                    if boxes[fdi - 1] is None:
                        excluded.append({'case': pid, 'fdi': fdi, 'reason': 'MISSING_TOOTH'})
                        continue
                    (result, exc) = tooth_measure(pid, fdi, tooth, mapped, sp, h, validate=len(controls) < 10, box=boxes[fdi - 1])
                    if exc:
                        excluded.append(exc)
                        continue
                    (m, control, crop) = result
                    source_fdi = fdi + 10 if fdi < 40 else fdi - 10
                    if mapping == 'identity':
                        source_fdi = fdi
                    m['source_pulpy_fdi'] = source_fdi
                    m['semantic_map'] = mapping
                    m['calibration_anchor'] = source_fdi in cal['anchors']
                    measure.append(m)
                    if control:
                        controls.append(control)
                    if saved < 3 and fdi % 10 == 6 and (not m['domain_truncated']):
                        file = DATA / f'R8_{pid}_{fdi}_paired.npz'
                        np.savez_compressed(file, **crop)
                        manifest.append({'path': str(file), 'sha256': hashlib.sha256(file.read_bytes()).hexdigest(), 'bytes': file.stat().st_size, 'case': pid, 'fdi': fdi})
                        saved += 1
            out.append(rec)
            for (name, obj) in [('calibration', out), ('excluded', excluded), ('controls', controls), ('manifest', manifest)]:
                (ROOT / 'raw' / f'R8_{name}.json').write_text(json.dumps(obj, indent=1) + '\n')
            write_csv(ROOT / 'raw/R8_teeth.csv', measure)
            if len(out) % 25 == 0:
                print('cases', len(out), 'accepted', sum((r.get('accepted', False) for r in out)), 'teeth', len(measure), flush=True)
            (ROOT / 'CURRENT_WORK_STATE.json').write_text(json.dumps({'lane': 'X12-pulpy3d', 'milestone': 'R8_SEMANTICALLY_CALIBRATED_ATLAS_RUNNING', 'latest_case': pid, 'accepted_cases': sum((r.get('accepted', False) for r in out)), 'measured_teeth': len(measure), 'last_gate': 'Two-anchor map; held-out95% FDI containment', 'next_operation': 'completeatlas; external anatomy compatibility; consumer ports and exports'}, indent=2) + '\n')
    cost = {'wall_s': time.time() - started, 'cpu_s': time.process_time() - cpu, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'source_stat': source_stat, 'cases_considered': len(out), 'accepted_cases': sum((r.get('accepted', False) for r in out)), 'measured_teeth': len(measure), 'fit': 'two-anchor discrete semantic calibration, no geometric fit', 'upstream_cost': 'raw/R6_cost.json fully charged'}
    (ROOT / 'raw/R8_cost.json').write_text(json.dumps(cost, indent=2) + '\n')
    print(json.dumps(cost), flush=True)
if __name__ == '__main__':
    main()
