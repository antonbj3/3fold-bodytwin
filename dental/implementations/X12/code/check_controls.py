from dental_release.paths import expand as _release_expand
import json, sys, time, os
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from archive_io import read_member, index
from consumer_ports import uniform_shell_decision, heat_transfer, compatible_anatomy
from pair_atlas_r6 import exact_match
PACKAGE = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('X12_RUN_ROOT', str(PACKAGE)))

def main():
    t0 = time.time()
    checks = []
    r = next((r for r in index() if r['name'] == 'Pulpy3D/dataset.json'))
    good = read_member(r)
    corrupt = r.copy()
    corrupt['crc32'] ^= 1
    try:
        read_member(corrupt)
        bad = False
    except ValueError:
        bad = True
    checks.append({'control': 'CRC integrity', 'valid_passed': bool(good), 'injected_error': 'CRC xor1', 'injected_rejected': bad})
    q = np.load(_release_expand('@DENTAL_WORK_ROOT@/X12/R6_P1_36_paired.npz'))
    t = q['tooth']
    p = q['pulp']
    sp = q['spacing']
    b = t & ~ndi.binary_erosion(t, structure=ndi.generate_binary_structure(3, 1))
    dist = ndi.distance_transform_edt(~b, sampling=sp)
    pc = np.argwhere(p)
    tree = cKDTree(np.argwhere(b) * sp)
    v = tree.query(pc * sp, workers=1)[0]
    dv = dist[tuple(pc.T)]
    err = float(np.max(np.abs(v - dv)))
    checks.append({'control': 'Independent EDT/KDTree', 'valid_passed': err <= 1e-06, 'max_abs_error_mm': err, 'injected_error': 'distance+1mm', 'injected_rejected': bool(np.max(np.abs(v - (dv + 1))) > 1e-06)})
    d2 = ndi.distance_transform_edt(~b, sampling=sp * 3)
    changed = float(np.max(np.abs(d2[tuple(pc.T)] - dv)))
    checks.append({'control': 'Physical scale', 'valid_passed': bool(np.allclose(sp, 0.3)), 'injected_error': '0.9mm rather than0.3mm spacing, x3', 'distance_change_mm': changed, 'injected_rejected': changed >= 0.5})
    import zipfile
    from tf2_io import ZIP, ROOT as TFROOT, read_mha_bytes
    from archive_io import nifti
    from pair_atlas_r6 import transform
    from calibrated_atlas import semantic_map
    source = next((r for r in index() if r['name'] == 'Pulpy3D/P1/gt_instance.nii.gz'))
    (img, _) = nifti(source)
    source_pulp = transform(np.asanyarray(img.dataobj), (2, 1, 0), (False, True, False))
    with zipfile.ZipFile(ZIP) as z:
        (full_tooth, _, _) = read_mha_bytes(z.read(TFROOT + '/labelsTr/ToothFairy2P_001.mha'))
    (validmap, cal) = semantic_map(source_pulp, full_tooth)
    wrong_tooth = full_tooth.copy()
    wrong_tooth[wrong_tooth == 46] = 47
    (badmap, badcal) = semantic_map(source_pulp, wrong_tooth)
    checks.append({'control': 'Canonical FDI membership', 'valid_passed': validmap == 'left_right_swap', 'injected_error': 'real held-out TF2 tooth46 relabelled47, anchors31/32 retained', 'injected_rejected': badmap is None, 'poison_reason': badcal.get('reason')})
    ct = np.arange(60).reshape(3, 4, 5)
    poison = ct.copy()
    poison[1, 2, 3] += 1
    checks.append({'control': 'Full CT equality', 'valid_passed': bool(exact_match(ct, ct)), 'injected_error': 'single CT voxel+1', 'injected_rejected': not exact_match(ct, poison), 'scope': 'Algorithm poison fixture; actual-data voxel poison also in R6 firstpair record'})
    a = uniform_shell_decision(2.0, 0.52, 0.5, 0.5)
    bq = uniform_shell_decision(0.2, 0.52, 0.5, 0.5)
    checks.append({'control': 'Inverse depth feasibility', 'valid_passed': a['decision'] == 'LOCAL_CONDITION_ALL_ENCLOSURES', 'injected_error': 'replace2mm thickness with0.2mm', 'injected_rejected': bq['decision'] == 'LOCAL_CONDITION_NO_ENCLOSURE'})
    refs = json.loads((PACKAGE / 'sources/anatomical_referents.json').read_text())['referents']
    wrong = next((q for q in refs if q['id'] == 'MAXILLARY2019_NEGATIVE_SCOPE'))
    checks.append({'control': 'External anatomical compatibility', 'valid_passed': compatible_anatomy(refs[0]['compared_quantity'], 'mandibular molars', refs[0]), 'injected_error': 'upper-molar reference assigned to lower FDI36', 'injected_rejected': not compatible_anatomy(wrong['compared_quantity'], 'mandibular first molars', wrong), 'scope': 'Landmark equality also required; all present projected-distance comparisons remain UNMATCHED'})
    checks.append({'control': 'Heat distance domain', 'valid_passed': heat_transfer(1, 0.2, 5) is not None, 'injected_error': 'negative remaining thickness', 'injected_rejected': heat_transfer(-1, 0.2, 5) is None})
    result = {'checks': checks, 'all_valid_passed': all((c['valid_passed'] for c in checks)), 'all_injected_rejected': all((c['injected_rejected'] for c in checks)), 'wall_s': time.time() - t0}
    (ROOT / 'raw/CONTROLS.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
    assert result['all_valid_passed'] and result['all_injected_rejected']
if __name__ == '__main__':
    main()
