"""X58: same-patient published canal revisions. No full-archive extraction."""
from dental_release.paths import expand as _release_expand
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, csv, datetime, hashlib, io, json, resource, subprocess, sys, time, zipfile
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
ROOT = Path(__file__).resolve().parents[1]
DENT = ROOT.parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X58-canal-wall-spread'))
TF1 = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/8838D60F38D5FBDE/ToothFairy_Dataset.zip'))
MAX = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/8838D60F38D5FBDE/IAN_Maxillo_dataset.zip'))
TF2 = Path(_release_expand('@DENTAL_CORPUS_ROOT@/geometry/ToothFairy2/ToothFairy2_Dataset.zip'))
SP = 0.3

def sha_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(4 << 20), b''):
            h.update(b)
    return h.hexdigest()

def dump(p, obj):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    tmp.replace(p)

def state(stage, gate, next_op, **kwargs):
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(lane='X58-canal-wall-spread', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), stage=stage, latest_gate=gate, next_operation=next_op, **kwargs))

def read_member(z, n):
    if z.getinfo(n).compress_type == 9:
        return subprocess.check_output(['7z', 'e', '-so', '-mmt=1', str(z.filename), n], stderr=subprocess.DEVNULL)
    return z.read(n)

def npy(z, n):
    b = read_member(z, n)
    a = np.load(io.BytesIO(b), allow_pickle=False)
    return (a, hashlib.sha256(b).hexdigest())

def mha(z, n):
    b = z.read(n)
    i = b.index(b'ElementDataFile')
    j = b.index(b'\n', i) + 1
    h = dict((s.split(' = ', 1) for s in b[:j].decode().strip().split('\n') if ' = ' in s))
    assert h.get('CompressedData', 'False') == 'False' and h['ElementDataFile'].strip() == 'LOCAL'
    dt = {'MET_UCHAR': 'u1', 'MET_CHAR': 'i1', 'MET_DOUBLE': 'f8', 'MET_FLOAT': 'f4', 'MET_SHORT': 'i2', 'MET_INT': 'i4', 'MET_UINT': 'u4', 'MET_USHORT': 'u2'}[h['ElementType']]
    shape = tuple((int(v) for v in h['DimSize'].split()[::-1]))
    sp = tuple((float(v) for v in h['ElementSpacing'].split()[::-1]))
    a = np.frombuffer(b, dtype='<' + dt, count=int(np.prod(shape)), offset=j).reshape(shape)
    return (a, sp, h, hashlib.sha256(b).hexdigest())

def boundary(m):
    return np.argwhere(m & ~ndi.binary_erosion(m, structure=ndi.generate_binary_structure(3, 1), border_value=0))

def dist(p, q):
    if not len(p) or not len(q):
        raise ValueError('Empty boundary')
    return cKDTree(q.astype(float) * SP).query(p.astype(float) * SP, workers=1)[0]

def stats(a):
    a = np.asarray(a)
    return {'n': len(a), 'median_mm': float(np.median(a)), 'p95_mm': float(np.quantile(a, 0.95)), 'max_mm': float(a.max())} if len(a) else {'n': 0, 'median_mm': None, 'p95_mm': None, 'max_mm': None}

def region_map(case, points):
    """Located NV1 mappings; endpoint windows are marked unvalidated anatomical proxies."""
    p = DENT / 'results/NV1_canals/per_case' / f'{case}.json'
    if not p.exists():
        return (np.full(len(points), -1, dtype=np.int16), np.full(len(points), -1, dtype=np.int8), {'status': 'NO_NV1_MAP'})
    nv = json.loads(p.read_text())
    C = []
    fdi = []
    seg = []
    for (sd, c) in nv.get('canals', {}).items():
        if not c.get('present') or not c.get('C'):
            continue
        cc = np.array(c['C'])
        ss = np.array(c['s'])
        reg = c['region']
        front = c.get('anterior_end')
        ai = int(np.argmin(np.linalg.norm(cc - np.array(front['point']), axis=1))) if front else None
        ds = np.abs(ss - ss[ai]) if ai is not None else None
        for (i, x) in enumerate(cc):
            fd = -1 if reg[i] is None else int(reg[i])
            sg = -1
            if ds is not None and ds[i] <= 5:
                sg = 0
            elif fd in (36, 37, 38, 46, 47, 48):
                sg = 1
            elif ds is not None and ds[i] >= ds.max() - 10:
                sg = 2
            C.append(x)
            fdi.append(fd)
            seg.append(sg)
    if not C:
        return (np.full(len(points), -1, dtype=np.int16), np.full(len(points), -1, dtype=np.int8), {'status': 'NO_CENTERLINE'})
    (d, k) = cKDTree(C).query(points * SP, workers=1)
    fd = np.array(fdi, dtype=np.int16)[k]
    sg = np.array(seg, dtype=np.int8)[k]
    fd[d > 3] = -1
    sg[d > 3] = -1
    return (fd, sg, {'status': 'NEAREST_NV1_CENTERLINE_WITHIN_3MM', 'nv1_sha256': sha_file(p), 'unresolved_points': int((sg < 0).sum()), 'total_points': len(points), 'mental_window_mm': 5, 'posterior_window_mm': 10, 'anatomical_validation': 'UNKNOWN'})

def verify_images(a, c, k):
    (old, oldsha) = npy(a, f'ToothFairy_Dataset/Dataset/{k}/data.npy')
    (new, sp, h, newsha) = mha(c, f'Dataset112_ToothFairy2/imagesTr/ToothFairy2P_{int(k[1:]):03d}_0000.mha')
    eq = old.shape == new.shape and np.array_equal(old, new)
    ans = {'exact_numeric_equal': bool(eq), 'tf1_image_npy_sha256': oldsha, 'tf2_image_mha_sha256': newsha, 'shape_tf1': list(old.shape), 'shape_tf2': list(new.shape), 'spacing_mm': sp, 'transform': h.get('TransformMatrix'), 'origin': h.get('Offset'), 'resolution': 'PER_POINT'}
    del old, new
    return ans

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args()
    for p in ROOT.glob('PREREG_*.json'):
        assert sha_file(p) == p.with_suffix('.sha256').read_text().split()[0]
    DATA.mkdir(parents=True, exist_ok=True)
    assert sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file())) < 3000000000
    start = time.perf_counter()
    rows = []
    failures = []
    regional = []
    manifest = []
    with zipfile.ZipFile(TF1) as a, zipfile.ZipFile(MAX) as b, zipfile.ZipFile(TF2) as c:
        aa = {i.filename.split('/')[-2]: i for i in a.infolist() if i.filename.endswith('gt_alpha.npy')}
        bb = {i.filename.split('/')[-2]: i for i in b.infolist() if i.filename.endswith('gt_alpha.npy')}
        cn = set(c.namelist())
        ai = {i.filename.split('/')[-2]: i for i in a.infolist() if i.filename.endswith('data.npy')}
        bi = {i.filename.split('/')[-2]: i for i in b.infolist() if i.filename.endswith('data.npy')}
        orig_identity = [{'patient': k, 'tf1_member': ai[k].filename, 'maxillo_member': v.filename, 'size_equal': v.file_size == ai[k].file_size, 'crc_equal': v.CRC == ai[k].CRC} for (k, v) in bi.items() if k in ai]
        identity = {'tf1_patients': len(ai), 'maxillo_patients': len(bi), 'dense_overlap': len(bb.keys() & aa.keys()), 'original_images': orig_identity, 'exact_sha_image_checks': [], 'independent_annotator_ids_present': False, 'independent_repeat_masks_documented': False, 'R1_gate': 'FAIL_INDEPENDENCE_NOT_IDENTIFIABLE', 'reason': '40 tool/protocol reannotations and 51 reused labels are published; independent six-patient repeat-annotation masks are not identified in either archive'}
        for k in sorted(bi, key=lambda s: int(s[1:]))[:3]:
            ba = read_member(a, ai[k].filename)
            ha = hashlib.sha256(ba).hexdigest()
            del ba
            ba = read_member(b, bi[k].filename)
            hb = hashlib.sha256(ba).hexdigest()
            del ba
            identity['exact_sha_image_checks'].append({'patient': k, 'tf1_sha256': ha, 'maxillo_sha256': hb, 'exact_byte_equal': ha == hb})
        dump(ROOT / 'raw/IDENTITY.json', identity)
        for (ii, k) in enumerate(sorted(bb, key=lambda s: int(s[1:])), 1):
            cache = ROOT / 'raw/patients' / f'{k}.json'
            if args.resume and cache.exists():
                row = json.loads(cache.read_text())
                ci = row.get('source_central_index')
                assert ci == {'new_crc': aa[k].CRC, 'new_bytes': aa[k].file_size, 'old_crc': bb[k].CRC, 'old_bytes': bb[k].file_size}, 'Archive member metadata changed: run fresh measurement'
                assert sha_file(row['point_artifact']['path']) == row['point_artifact']['sha256'], 'Point artifact hash drift'
                rows.append(row)
                regional.extend(row['regional'])
                manifest.append(row['point_artifact'])
                continue
            t = time.perf_counter()
            try:
                (old, oh) = npy(b, bb[k].filename)
                (new, nh) = npy(a, aa[k].filename)
                assert old.shape == new.shape
                om = old > 0
                nm = new > 0
                del old, new
                equal = bool(np.array_equal(om, nm))
                po = boundary(om)
                pn = boundary(nm)
                do = dist(po, pn)
                dn = dist(pn, po)
                sym = np.r_[do, dn]
                inter = int((om & nm).sum())
                no = int(om.sum())
                nn = int(nm.sum())
                added = nn - inter
                row = {'patient': k, 'case_tf2': f'ToothFairy2P_{int(k[1:]):03d}', 'resolution': 'PER_SURFACE_REGION', 'time_scale': 'SIMULTANEOUS', 'old_member': bb[k].filename, 'new_member': aa[k].filename, 'old_mask_sha256': oh, 'new_mask_sha256': nh, 'numeric_mask_equal': equal, 'shape': list(om.shape), 'spacing_mm': [SP] * 3, 'old_voxels': no, 'new_voxels': nn, 'intersection_voxels': inter, 'added_voxels': added, 'author_added_fraction_percent': 100 * added / no, 'dice': 2 * inter / (no + nn), 'old_to_new': stats(do), 'new_to_old': stats(dn), 'symmetric': stats(sym), 'regional': []}
                case = row['case_tf2']
                tm = f'Dataset112_ToothFairy2/labelsTr/{case}.mha'
                image_ok = False
                if tm in cn:
                    vi = verify_images(a, c, k)
                    row['image_identity'] = vi
                    image_ok = vi['exact_numeric_equal']
                    assert all((abs(x - SP) < 1e-14 for x in vi['spacing_mm']))
                    (label, sp, h, lh) = mha(c, tm)
                    row['tf2_label_sha256'] = lh
                    mm = np.isin(label, [3, 4])
                    del label
                    row['tf2_mask_equal_tf1'] = bool(mm.shape == nm.shape and np.array_equal(mm, nm))
                    row['tf2_mask_equal_maxillo'] = bool(mm.shape == om.shape and np.array_equal(mm, om))
                else:
                    mm = None
                    row['image_identity'] = {'exact_numeric_equal': False, 'status': 'MISSING_TF2'}
                if image_ok:
                    (fd_o, sg_o, map_o) = region_map(case, po)
                    (fd_n, sg_n, map_n) = region_map(case, pn)
                else:
                    fd_o = np.full(len(po), -1, dtype=np.int16)
                    fd_n = np.full(len(pn), -1, dtype=np.int16)
                    sg_o = np.full(len(po), -1, dtype=np.int8)
                    sg_n = np.full(len(pn), -1, dtype=np.int8)
                    map_o = map_n = {'status': 'NO_VERIFIED_IDENTITY'}
                row['region_mapping'] = {'old': map_o, 'new': map_n, 'clinical_segment_validity': 'UNVALIDATED_PROXIES'}
                for (sg, name) in [(0, 'mental_endpoint_5mm_proxy'), (1, 'molar_nearest_FDI'), (2, 'posterior_endpoint_10mm_proxy'), (-1, 'UNRESOLVED')]:
                    ds = np.r_[do[sg_o == sg], dn[sg_n == sg]]
                    rr = {'patient': k, 'case': case, 'segment': name, 'revised': not equal, 'resolution': 'PER_SURFACE_REGION', **stats(ds)}
                    row['regional'].append(rr)
                    regional.append(rr)
                for fd in sorted(set(fd_o) | set(fd_n)):
                    if fd < 0:
                        continue
                    row['regional'].append({'patient': k, 'case': case, 'fdi': int(fd), 'segment': 'nearest_FDI', 'revised': not equal, 'resolution': 'PER_TOOTH', **stats(np.r_[do[fd_o == fd], dn[fd_n == fd]])})
                artifact = DATA / f'{k}_walls.npz'
                data = {'old_boundary_zyx': po.astype(np.int16), 'new_boundary_zyx': pn.astype(np.int16), 'old_to_new_mm': do, 'new_to_old_mm': dn, 'old_region': sg_o, 'new_region': sg_n, 'old_fdi': fd_o, 'new_fdi': fd_n, 'old_voxels_zyx': np.argwhere(om).astype(np.int16), 'new_voxels_zyx': np.argwhere(nm).astype(np.int16)}
                if image_ok and mm is not None:
                    data['tf2_voxels_zyx'] = np.argwhere(mm).astype(np.int16)
                np.savez_compressed(artifact, **data)
                row['point_artifact'] = {'path': str(artifact), 'sha256': sha_file(artifact), 'bytes': artifact.stat().st_size, 'resolution': 'PER_POINT'}
                row['source_central_index'] = {'new_crc': aa[k].CRC, 'new_bytes': aa[k].file_size, 'old_crc': bb[k].CRC, 'old_bytes': bb[k].file_size}
                row['wall_seconds'] = time.perf_counter() - t
                rows.append(row)
                manifest.append(row['point_artifact'])
                dump(cache, row)
                del om, nm, mm, po, pn, do, dn, data
                if ii % 5 == 0 or ii == len(bb):
                    dump(ROOT / 'raw/PAIRED_PATIENTS.json', rows)
                    dump(ROOT / 'raw/POINT_MANIFEST.json', manifest)
                    state('R2_MEASURING', 'R1_INDEPENDENCE_FAIL', 'Continue dense revisions and exact TF2 image identity', patients_complete=len(rows), patients_total=len(bb), elapsed_seconds=time.perf_counter() - start)
                print(f"{ii}/{len(bb)} {k} equal={equal} P95={row['symmetric']['p95_mm']:.4f} image_equal={image_ok} {row['wall_seconds']:.1f}s", flush=True)
            except Exception as e:
                import traceback
                failures.append({'patient': k, 'error': repr(e), 'traceback': traceback.format_exc()})
                dump(ROOT / 'raw/FAILURES.json', failures)
                print(k, repr(e), flush=True)
        dump(ROOT / 'raw/PAIRED_PATIENTS.json', rows)
        dump(ROOT / 'raw/POINT_MANIFEST.json', manifest)
        dump(ROOT / 'raw/FAILURES.json', failures)
        rev = [r for r in rows if not r['numeric_mask_equal']]
        same = [r for r in rows if r['numeric_mask_equal']]
        mean_added = float(np.mean([r['author_added_fraction_percent'] for r in rev])) if rev else None
        out = {'claim_type': 'information_link', 'completed': len(rows), 'total_dense_pairs': len(bb), 'revised_patients': len(rev), 'copied_patients': len(same), 'image_matched_tf2': sum((r['image_identity']['exact_numeric_equal'] for r in rows)), 'tf2_equal_tf1': sum((r.get('tf2_mask_equal_tf1', False) for r in rows)), 'median_of_patient_median_mm': float(np.median([r['symmetric']['median_mm'] for r in rev])) if rev else None, 'median_of_patient_p95_mm': float(np.median([r['symmetric']['p95_mm'] for r in rev])) if rev else None, 'max_patient_p95_mm': max((r['symmetric']['p95_mm'] for r in rev), default=None), 'mean_author_added_percent': mean_added, 'patients_p95_gt_0p3mm': sum((r['symmetric']['p95_mm'] > 0.3 for r in rev)), 'gates': {'published_40_revisions_51_copies': len(rev) == 40 and len(same) == 51, 'published_added_61p9_percent': mean_added is not None and abs(mean_added - 61.9) <= 0.5, 'consequential_wall_revision': len(rev) > 0 and sum((r['symmetric']['p95_mm'] > 0.3 for r in rev)) / len(rev) >= 0.05, 'independent_annotator_floor': False}, 'dropout': {'eligible_pairs': len(bb), 'failed': len(failures), 'failed_fraction': len(failures) / len(bb), 'independent_pairs_admitted': 0, 'independent_pair_rejection_fraction': 1.0, 'rejection_reason': 'missing independent/blinded annotator attribution'}, 'cost': {'wall_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'gpu': False, 'scratch_bytes': sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file())), 'fit_seconds': 0, 'questions': 0}}
        dump(ROOT / 'raw/R2_RESULT.json', out)
        state('R2_COMPLETE', out['gates'], 'R3 located clearance transfer; retain failed independence and published-count gates', patients_complete=len(rows))
        print(json.dumps(out, indent=2), flush=True)
if __name__ == '__main__':
    main()
