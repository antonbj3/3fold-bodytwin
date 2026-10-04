from dental_release.paths import expand as _release_expand
import argparse, csv, hashlib, json, os, resource, sys, time, zipfile
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from operators import measure, H
ROOT = Path(__file__).resolve().parents[1]
X12 = ROOT.parent / _release_expand('X12')
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X30'))
sys.path.insert(0, str(X12 / 'code'))
from archive_io import index, nifti, ARCHIVE
from tf2_io import ZIP, ROOT as TFROOT, read_mha_bytes

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, o):
    Path(p).write_text(json.dumps(o, indent=2, allow_nan=False) + '\n')

def csvwrite(p, rows):
    if not rows:
        return
    with Path(p).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--tag', default='R1')
    args = ap.parse_args()
    start = time.monotonic()
    cpu = time.process_time()
    tag = args.tag
    pairs = {r['case']: r for r in json.loads((X12 / 'raw/R6_pairs.json').read_text())}
    calibration = json.loads((X12 / 'raw/R8_calibration.json').read_text())
    rows = list(csv.DictReader((X12 / 'raw/R8_teeth.csv').open()))
    targets = {}
    for r in rows:
        if r['domain_truncated'] == 'False':
            targets.setdefault(r['case'], []).append(r)
    for (path, stat) in json.loads((X12 / 'raw/R6_source_stat.json').read_text()).items():
        p = Path(path)
        if dict(size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns) != stat:
            raise RuntimeError('Source archive stat changed')
    lookup = {r['name']: r for r in index()}
    out = []
    paths = []
    sections = []
    graphs = []
    controls = []
    manifest = []
    excluded = []
    with zipfile.ZipFile(ZIP) as zf:
        selected = [r for r in calibration if r['accepted'] and r['case'] in targets]
        if args.limit:
            selected = selected[:args.limit]
        for (ri, cal) in enumerate(selected):
            case = cal['case']
            pair = pairs[case]
            (pimg, psha) = nifti(lookup[f'Pulpy3D/{case}/gt_instance.nii.gz'])
            if psha != cal['pulp_sha256']:
                raise RuntimeError('Pulpy SHA drift')
            pulp = np.asanyarray(pimg.dataobj)
            (perm, flip) = pair['image_transforms'][0]
            pulp = pulp.transpose(perm)[tuple((slice(None, None, -1) if f else slice(None) for f in flip))]
            raw = zf.read(f"{TFROOT}/labelsTr/{pair['tf2_case']}.mha")
            if hashlib.sha256(raw).hexdigest() != cal['tooth_sha256']:
                raise RuntimeError('TF2 SHA drift')
            (tooth, spacing, header) = read_mha_bytes(raw)
            if tuple(spacing) != (H, H, H):
                raise RuntimeError('Invalid physical scale')
            if not np.allclose(np.fromstring(header['TransformMatrix'], sep=' ').reshape(3, 3)[2], [0, 0, 1]):
                raise RuntimeError('Superior-axis unsupported')
            for tr in targets[case]:
                fdi = int(tr['fdi'])
                source = int(tr['source_pulpy_fdi'])
                tc = np.argwhere(tooth == fdi)
                lo = np.maximum(0, tc.min(0) - 3)
                hi = np.minimum(tooth.shape, tc.max(0) + 4)
                sl = tuple((slice(int(a), int(b)) for (a, b) in zip(lo, hi)))
                t = tooth[sl] == fdi
                p = pulp[sl] == source
                containment = float((p & t).sum() / max(1, p.sum()))
                if containment < 0.95:
                    raise RuntimeError('Inherited FDI pairing fails')
                (rec, graph, pc, bc, d) = measure(p, t, lo)
                rec = {'case': case, 'fdi': fdi, 'source_fdi': source, 'tooth_type': tr['tooth_type'], 'calibration_anchor': tr['calibration_anchor'], 'containment': containment, 'pulp_sha256': psha, 'tooth_sha256': cal['tooth_sha256'], 'voxel_mm': H, 'resolution_level': 'PER_TOOTH', **rec}
                out.append(rec)
                graphs.append({'case': case, 'fdi': fdi, 'z_range': graph['z_range'], 'counts': graph['counts'], 'edges': graph['edges']})
                for (j, path) in enumerate(graph['paths']):
                    paths.append({'case': case, 'fdi': fdi, 'branch_id': j, 'resolution_level': 'PER_SURFACE_REGION', **path})
                for node in graph['nodes']:
                    sections.append({'case': case, 'fdi': fdi, 'resolution_level': 'PER_SURFACE_REGION', **node})
                if len(controls) < 10 and len(d):
                    kd = cKDTree(bc * H).query(pc * H, workers=1)[0]
                    controls.append({'case': case, 'fdi': fdi, 'max_abs_error_mm': float(np.max(abs(d - kd))), 'gate_pass': bool(np.max(abs(d - kd)) <= 1e-06), 'poison_distance_rejected': bool(np.max(abs(d + 1 - kd)) > 1e-06)})
                cp = DATA / f'{case}_{fdi}.npz'
                if not cp.exists():
                    np.savez_compressed(cp, pulp=p, tooth=t, origin_zyx=lo, spacing=np.array(spacing))
                manifest.append({'case': case, 'fdi': fdi, 'path': str(cp), 'bytes': cp.stat().st_size, 'sha256': sha(cp)})
            mask = np.isin(tooth, [3, 4])
            endpoints = [(r, rr['apical_annotation_endpoint_mm']) for r in paths if r['case'] == case for rr in [r]]
            if mask.any() and endpoints:
                mb = mask & ~ndi.binary_erosion(mask, structure=ndi.generate_binary_structure(3, 1), border_value=0)
                tree = cKDTree(np.argwhere(mb) * H)
                for (r, endpoint) in endpoints:
                    r['annotation_endpoint_IAN_boundary_mm'] = float(tree.query(endpoint)[0])
                    r['endpoint_IAN_digital_halfwidth_mm'] = float(np.sqrt(3) * H)
            for r in paths:
                if r['case'] == case:
                    r.setdefault('annotation_endpoint_IAN_boundary_mm', None)
            csvwrite(ROOT / 'raw' / f'{tag}_teeth.csv', out)
            if (ri + 1) % 10 == 0:
                print('cases', ri + 1, 'teeth', len(out), 'paths', len(paths), flush=True)
            dump(ROOT / 'CURRENT_WORK_STATE.json', {'lane': 'X30-root-canal', 'stage': tag + '_ATLAS_RUNNING', 'last_case': case, 'cases_complete': ri + 1, 'teeth': len(out), 'last_gate': 'source CRC/SHA, containment and 0.3-mm grid checked', 'next_operation': 'complete population; test external prevalence; freeze uncertainty successor'})
    csvwrite(ROOT / 'raw' / f'{tag}_teeth.csv', out)
    with (DATA / f'{tag}_paths.jsonl').open('w') as f:
        for r in paths:
            f.write(json.dumps(r, allow_nan=False) + '\n')
    with (DATA / f'{tag}_sections.jsonl').open('w') as f:
        for r in sections:
            f.write(json.dumps(r, allow_nan=False) + '\n')
    with (DATA / f'{tag}_graphs.jsonl').open('w') as f:
        for r in graphs:
            f.write(json.dumps(r, allow_nan=False) + '\n')
    dump(ROOT / 'raw' / f'{tag}_manifest.json', manifest)
    dump(ROOT / 'raw' / f'{tag}_distance_controls.json', controls)
    cost = {'wall_s': time.monotonic() - start, 'cpu_s': time.process_time() - cpu, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'cases': len(selected), 'teeth': len(out), 'paths': len(paths), 'sections': len(sections), 'GPU_s': 0, 'fit_s': 0, 'parent_R6_cost': json.loads((X12 / 'raw/R6_cost.json').read_text()), 'parent_R8_cost': json.loads((X12 / 'raw/R8_cost.json').read_text())}
    dump(ROOT / 'raw' / f'{tag}_cost.json', cost)
    print(json.dumps(cost), flush=True)
if __name__ == '__main__':
    main()
