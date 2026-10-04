import csv, json, resource, time, hashlib, sys
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from uncertainty import counts_in_region
from atlas import ROOT, DATA, X12, sha, dump, csvwrite
sys.path.insert(0, str(X12 / 'code'))
from archive_io import index, nifti

def main():
    start = time.monotonic()
    cpu = time.process_time()
    manifests = json.loads((ROOT / 'raw/R1_manifest.json').read_text())
    pairs = {r['case']: r for r in json.loads((X12 / 'raw/R6_pairs.json').read_text())}
    teeth = list(csv.DictReader((ROOT / 'raw/R1_teeth.csv').open()))
    lookup = {(r['case'], int(r['fdi'])): r for r in teeth}
    inc = [r for r in manifests if r['fdi'] in [31, 32, 41, 42]]
    rows = []
    for row in inc:
        if sha(row['path']) != row['sha256']:
            raise RuntimeError('Crop hash changed')
        with np.load(row['path']) as a:
            p = a['pulp']
        z = np.argwhere(p)[:, 0]
        (count, profile) = counts_in_region(p, int(z.min()), int(z.max()))
        rows.append({'case': row['case'], 'fdi': row['fdi'], 'resolution_level': 'PER_TOOTH', 'rootward_R1_count': int(lookup[row['case'], row['fdi']]['persistent_section_count']), 'whole_pulp_persistent_count': count, 'all_plane_counts': json.dumps(profile), 'whole_pulp_is_not_root_landmark': 'True'})
    csvwrite(ROOT / 'raw/R4_full_pulp_incissors.csv', rows)
    summary = []
    expected = {31: 0.3, 32: 0.365, 41: 0.335, 42: 0.335}
    for (f, p) in expected.items():
        rs = [r for r in rows if r['fdi'] == f]
        k = sum((r['whole_pulp_persistent_count'] >= 2 for r in rs))
        summary.append({'fdi': f, 'n': len(rs), 'second_count': k, 'fraction': k / len(rs), 'external_fraction': p, 'external15pp_gate_pass': abs(k / len(rs) - p) <= 0.15, 'resolution_level': 'POPULATION'})
    available = {}
    for r in inc:
        available.setdefault(r['case'], set()).add(r['fdi'])
    selected = [c for c in available if len(available[c]) == 4][:8]
    pred = ROOT / 'FROZEN_PREDICTIONS_CT_REVIEW.json'
    if not pred.exists():
        import datetime
        frozen = {'frozen_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'claim_type': 'capability', 'review_state': 'INDEPENDENT_REVIEW_NOT_PERFORMED', 'cases': selected, 'quantity': 'complete-mask axial section multiplicity, not clinically established canal count', 'predictions': [r for r in rows if r['case'] in selected], 'tooth_pulp_hashes': [r for r in inc if r['case'] in selected], 'clinical_working_length': 'UNKNOWN', 'Schneider_Pruett_landmarks': 'UNKNOWN', 'anatomical_apex': 'UNKNOWN', 'blinding': 'Reviewer should receive CT crops without these predictions/overlays; coordinator unblinds afterwards'}
        dump(pred, frozen)
        pred.with_suffix('.sha256').write_text(sha(pred) + '  ' + pred.name + '\n')
    else:
        if sha(pred) != pred.with_suffix('.sha256').read_text().split()[0]:
            raise RuntimeError('Frozen predictions drift')
        if json.loads(pred.read_text())['cases'] != selected:
            raise RuntimeError('Frozen review selection drift')
    archive = {r['name']: r for r in index()}
    out = []
    mutation = []
    for case in selected:
        (img, digest) = nifti(archive[f'Pulpy3D/{case}/data.nii.gz'])
        if digest != pairs[case]['pulpy_CT_sha256']:
            raise RuntimeError('CT SHA drift')
        im = np.asanyarray(img.dataobj)
        (perm, flip) = pairs[case]['image_transforms'][0]
        im = im.transpose(perm)[tuple((slice(None, None, -1) if f else slice(None) for f in flip))]
        for row in [r for r in inc if r['case'] == case]:
            with np.load(row['path']) as a:
                origin = a['origin_zyx']
                shape = a['pulp'].shape
                spacing = a['spacing']
            sl = tuple((slice(int(o), int(o + s)) for (o, s) in zip(origin, shape)))
            crop = im[sl].copy()
            cp = DATA / f"REVIEW_CT_{case}_{row['fdi']}.npz"
            if not cp.exists():
                np.savez_compressed(cp, CT=crop, origin_zyx=origin, spacing=spacing)
            with np.load(cp) as a:
                valid = np.array_equal(a['CT'], im[sl])
                poison = not np.array_equal(np.roll(a['CT'], 1, axis=2), im[sl])
            if not valid:
                raise RuntimeError('CT crop differs from source')
            mutation.append({'case': case, 'fdi': row['fdi'], 'valid_crop_equal': valid, 'one_voxel_shift_rejected': poison, 'wrong_source_sha_rejected': digest != '0' * 64})
            out.append({'case': case, 'fdi': row['fdi'], 'path': str(cp), 'sha256': sha(cp), 'bytes': cp.stat().st_size, 'source_CT_sha256': digest, 'origin_zyx': origin.tolist(), 'spacing_mm': spacing.tolist(), 'coordinate_frame': 'TF2 full-volume index(zyx);world origin/direction in sourceMHA', 'review_performed': False})
    dump(ROOT / 'raw/R4_CT_review_manifest.json', out)
    res = {'claim_type': 'capability', 'full_pulp_population_summary': summary, 'crop_explanation_gate_pass': all((r['external15pp_gate_pass'] for r in summary)), 'raw_CT_review_cases': len(selected), 'raw_CT_review_teeth': len(out), 'review_performed': False, 'crop_identity_controls': mutation, 'all_valid_crop_passed': all((r['valid_crop_equal'] for r in mutation)), 'all_injected_errors_rejected': all((r['one_voxel_shift_rejected'] and r['wrong_source_sha_rejected'] for r in mutation)), 'wall_s': time.monotonic() - start, 'cpu_s': time.process_time() - cpu, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024}
    dump(ROOT / 'raw/R4_result.json', res)
    print(json.dumps(res), flush=True)
if __name__ == '__main__':
    main()
