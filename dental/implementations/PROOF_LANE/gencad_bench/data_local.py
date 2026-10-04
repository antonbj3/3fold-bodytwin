"""Bounded local-data preparation. No downloads and no cross-subject bite joining."""
from dental_release.paths import expand as _release_expand
import json, time
from pathlib import Path
import hashlib
import numpy as np
from .io import ROOT, DATA, dump, sha, digest, freeze
from .geometry import profile_from_points, dome_profile, rounded
TDS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/geometry/Teeth3DS'))
TYPES = [1, 3, 4, 6]

def files(split):
    names = (TDS / 'Teeth3DS_train_test_split' / f'{split}_upper.txt').read_text().splitlines()
    rows = []
    for name in names:
        name = name.strip()
        pid = name.rsplit('_', 1)[0]
        for part in range(1, 8):
            p = TDS / f'data_part_{part}' / 'upper' / pid / (name + '.obj')
            if p.exists() and p.with_suffix('.json').exists():
                rows.append((pid, p))
                break
    return sorted(rows, key=lambda x: hashlib.sha256(x[0].encode()).hexdigest())

def read_obj(p):
    vv = []
    ff = []
    for line in p.open():
        if line.startswith('v '):
            vv.append([float(x) for x in line.split()[1:4]])
        elif line.startswith('f '):
            ids = [int(x.split('/')[0]) - 1 for x in line.split()[1:]]
            if min(ids) < 0:
                raise ValueError('Unsupported negative OBJ indices')
            for j in range(1, len(ids) - 1):
                ff.append([ids[0], ids[j], ids[j + 1]])
    return (np.asarray(vv), np.asarray(ff, dtype=np.int32))

def arch(p):
    labels = np.asarray(json.loads(p.with_suffix('.json').read_text())['labels'])
    if any((np.sum(labels == q * 10 + t) < 100 for q in [1, 2] for t in TYPES)):
        return None
    (v, f) = read_obj(p)
    if len(v) != len(labels):
        raise ValueError('Label/vertex mismatch')
    z = np.linalg.eigh(np.cov(v.T))[1][:, 0]
    if np.any(labels == 0) and np.dot(z, v[labels > 0].mean(0) - v[labels == 0].mean(0)) < 0:
        z = -z
    x = np.array([1.0, 0.0, 0.0])
    x -= np.dot(x, z) * z
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    frame = np.array([x, y, z])
    v = v @ frame.T
    teeth = {}
    for fdi in [10 + t for t in TYPES] + [20 + t for t in TYPES]:
        ids = np.flatnonzero(labels == fdi)
        vv = v[ids]
        origin = np.r_[np.median(vv[:, :2], axis=0), np.quantile(vv[:, 2], 0.02)]
        vv = vv - origin
        keep = np.all(labels[f] == fdi, axis=1) & np.all(v[f, 2] >= origin[2], axis=1)
        lookup = np.full(len(v), -1, int)
        lookup[ids] = np.arange(len(ids))
        ff = lookup[f[keep]]
        if len(ff) < 20:
            return None
        teeth[fdi] = dict(vertices=vv, faces=ff, origin=origin, profile=profile_from_points(vv))
    return (teeth, frame)

def prepare():
    tic = time.perf_counter()
    DATA.mkdir(parents=True, exist_ok=True)
    paths = {'training': files('training'), 'testing': files('testing')}
    training = []
    evaluation = []
    manifest = []
    for (split, need, out) in [('training', 6, training), ('testing', 2, evaluation)]:
        for (pid, p) in paths[split]:
            a = arch(p)
            if a is None:
                continue
            (teeth, frame) = a
            manifest.append(dict(subject=pid, split=split, path=str(p), sha256=sha(p), labels_sha256=sha(p.with_suffix('.json'))))
            out.append((pid, teeth, frame, p))
            if len(out) == need:
                break
        if len(out) != need:
            raise ValueError('Insufficient locally eligible ' + split + ' subjects')
    assert not {r[0] for r in training} & {r[0] for r in evaluation}
    means = {}
    for t in TYPES:
        ps = [a[1][10 + t]['profile'] for a in training]
        means[t] = dict(z_mm=rounded(np.mean([x['z_mm'] for x in ps], 0)), radii_mm=rounded(np.mean([x['radii_mm'] for x in ps], 0)), apex_mm=round(float(np.mean([x['apex_mm'] for x in ps])), 6))
    tasks = []
    private = {}
    for (idx, (pid, teeth, frame, src)) in enumerate(evaluation):
        for ty in TYPES:
            fdi = 10 + ty
            target = teeth[fdi]
            donor = teeth[20 + ty]
            dv = donor['vertices'].copy()
            dv[:, 0] *= -1
            mirrored = profile_from_points(dv)
            tp = target['profile']
            R = 0.55 * float(np.median(tp['radii_mm'][0]))
            H = 0.62 * tp['apex_mm']
            prep = dome_profile(R, H)
            inner = dome_profile(R + 0.08, H + 0.08)
            gid = f's{idx + 1}_FDI{fdi}'
            private_file = DATA / (gid + '_reference.npz')
            np.savez_compressed(private_file, vertices=target['vertices'], faces=target['faces'])
            private[gid] = dict(path=str(private_file), sha256=sha(private_file), source_subject=pid, source_path=str(src), source_sha256=sha(src), FDI=fdi, frame_rotation=frame.tolist(), origin=target['origin'].tolist())
            mats = [('KATANA_ML', 0.4 if ty <= 3 else 0.5), ('KATANA_UTML', 0.8 if ty <= 3 else 1.0), ('IPS_e_max_CAD_adhesive_1mm', 1.0)]
            for (mat, w) in mats:
                task = dict(schema_version='0.1', task_id=gid + '_' + mat, geometry_id=gid, subject_group=hashlib.sha256(pid.encode()).hexdigest()[:16], split='development' if idx == 0 else 'test', units='mm', frame='local_crown_mm', tooth_FDI=fdi, tooth_type=ty, geometry_reference=dict(dataset='Teeth3DS', kind='withheld_measured_surface', locator='https://osf.io/xctdy/', sha256=private[gid]['sha256']), preparation=dict(kind='synthetic_dome', profile=prep, physical_measurement=False), material=dict(product=mat, requirement_source='data/materials.json#' + mat), requirements=dict(wall_min_mm=w, cement_min_mm=0.03, mill_radius_mm=0.5, mill_allowance_mm=0.05, allowed_exit_directions=[[0, 0, -1]], required_checks=['geometry', 'insertion', 'wall', 'cement', 'milling', 'occlusion']), truth_tiers=['exact_digital_scope', 'withheld_anatomy'], antagonist=dict(status='UNKNOWN', reason='No measured registered opposing arch for this subject'), public=dict(contralateral_profile=mirrored, population_profile=means[ty], intaglio_template=inner, population_training_groups=[hashlib.sha256(a[0].encode()).hexdigest()[:16] for a in training]), information_limit='Synthetic prep discloses target size; no native target surface to plugin')
                tasks.append(task)
    dump(ROOT / 'data/tasks.json', tasks)
    dump(ROOT / 'data/evaluator_private.json', private)
    record = dict(tasks_sha256=sha(ROOT / 'data/tasks.json'), references=private, source_manifest=manifest, n_tasks=len(tasks), split_unit='subject', preparation_wall_s=time.perf_counter() - tic)
    dump(ROOT / 'raw/data_preparation.json', record)
    frozen = {k: v for (k, v) in record.items() if k != 'preparation_wall_s'}
    freeze(ROOT / 'data/TASKSET_FROZEN.json', frozen)
    return tasks
