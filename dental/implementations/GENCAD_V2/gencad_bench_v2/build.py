"""Local-only cohort construction. Target references never enter the public scene."""
from dental_release.paths import expand as _release_expand
import sys, time, zipfile, hashlib
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from .common import *
from .geometry import *
sys.path.insert(0, str(ROOT / 'gencad_bench_v2/vendor'))
from obstacle import constraints
LABEL = Path(_release_expand('@DENTAL_WORK_ROOT@/X11/targets'))
ZIPS = {d: Path(_release_expand('@DENTAL_DATA_ROOT@/geometry')) / d / f for (d, f) in [('Bits2Bites', 'Bits2Bites_v01.zip'), ('Bite2Text', 'Bite2Text.zip')]}
FDI = {'molar_crown': 36, 'premolar_crown': 34, 'anterior_crown': 31, 'bridge3': 35, 'implant_crown': 46, 'inlay': 36, 'onlay': 46, 'veneer': 31, 'lattice_onlay': 37}

def select():
    selected = []
    all_hashes = {}
    excluded = []
    for dataset in ZIPS:
        records = sorted((p for p in (LABEL / dataset).iterdir() if p.is_dir()), key=lambda p: digest('PROOF_LANE-v2|' + dataset + '|' + p.name))
        n = 0
        for p in records:
            meta = read(p / 'labels+landmarks.json')
            if meta.get('input_status') == 'REFUSED_EMPTY_SOURCE':
                continue
            pairhash = tuple((meta['arches'][j]['input_sha256'] for j in ['upper', 'lower']))
            if pairhash in all_hashes:
                excluded.append(dict(dataset=dataset, case=p.name, duplicate_of=all_hashes[pairhash]))
                continue
            split = ['train', 'dev', 'test'][n // 4]
            case_key = digest(dataset + '|' + p.name)[:16]
            row = dict(dataset=dataset, case=p.name, case_key=case_key, split=split, metadata=str(p / 'labels+landmarks.json'), metadata_sha256=sha(p / 'labels+landmarks.json'), geometry_hashes=pairhash)
            selected.append(row)
            all_hashes[pairhash] = case_key
            n += 1
            if n == 12:
                break
        if n != 12:
            raise ValueError('Insufficient local metadata')
    freeze(ROOT / 'data/COHORT_LOCK.json', dict(cases=selected, excluded_exact_duplicates=excluded, latent_cross_dataset_identity='UNKNOWN'))
    return selected

def pair(row):
    meta = read(row['metadata'])
    out = {}
    sources = []
    if sha(row['metadata']) != row['metadata_sha256']:
        raise ValueError('label metadata drift')
    with zipfile.ZipFile(ZIPS[row['dataset']]) as z:
        for jaw in ['upper', 'lower']:
            suffix = '/' + row['case'] + '/' + jaw + '.stl' if row['dataset'] == 'Bits2Bites' else row['case'] + '/ios/ios_' + jaw + '.stl'
            matches = [s for s in z.namelist() if s.endswith(suffix)]
            if len(matches) != 1:
                raise ValueError('Missing/nonunique geometry member ' + suffix)
            member = matches[0]
            blob = z.read(member)
            info = meta['arches'][jaw]
            if hashlib.sha256(blob).hexdigest() != info['input_sha256']:
                raise ValueError('scan hash mismatch')
            (v, f) = parse_stl(blob)
            lp = Path(row['metadata']).parent / (jaw + '_labels.npz')
            if sha(lp) != info['labels']['sha256']:
                raise ValueError('label hash mismatch')
            labs = np.load(lp, allow_pickle=False)['labels']
            if len(labs) != len(v):
                raise ValueError('label indexing')
            fl = labs[f]
            owner = np.where(fl[:, 0] == fl[:, 1], fl[:, 0], np.where(fl[:, 0] == fl[:, 2], fl[:, 0], np.where(fl[:, 1] == fl[:, 2], fl[:, 1], 0)))
            out[jaw] = dict(v=v, f=f, labels=labs, owner=owner, meta=info)
            sources.append(dict(member=member, zip_path=str(ZIPS[row['dataset']]), bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest(), label_sha256=sha(lp)))
    return (out, sources)

def make_site(pair, family):
    arch = pair['lower']
    fdi = FDI[family]
    pp = arch['v'][arch['labels'] == fdi]
    if len(pp) < 100:
        raise ValueError('missing predicted tooth ' + str(fdi))
    fr = arch['meta']['frame']
    occ = unit(fr['occlusal_unit'])
    right = unit(fr['right_unit'])
    x = unit(right - occ * (right @ occ))
    y = unit(np.cross(occ, x))
    R = np.stack([x, y, occ], axis=1)
    base = pp.mean(0)
    if family == 'veneer':
        lm = arch['meta']['landmarks']['teeth'].get(str(fdi), {})
        normal = unit(lm['buccal_unit'])
        x = unit(np.cross(occ, normal))
        y = unit(np.cross(normal, x))
        R = np.stack([x, y, normal], axis=1)
    elif family == 'bridge3':
        a = arch['v'][arch['labels'] == 34]
        b = arch['v'][arch['labels'] == 36]
        if min(len(a), len(b)) < 100:
            raise ValueError('bridge abutments missing')
        x = unit(b.mean(0) - a.mean(0))
        x = unit(x - occ * (x @ occ))
        y = unit(np.cross(occ, x))
        R = np.stack([x, y, occ], axis=1)
        pp = arch['v'][np.isin(arch['labels'], [34, 35, 36])]
        base = pp.mean(0)
    local = (pp - base) @ R
    band = local[local[:, 2] <= np.quantile(local[:, 2], 0.35)] if family != 'veneer' else local
    (lo, hi) = np.quantile(band[:, :2], [0.02, 0.98], axis=0)
    center = (hi + lo) / 2
    half = (hi - lo) / 2
    if np.min(half) < 0.7 or np.max(half) > 18:
        raise ValueError('invalid predicted footprint')
    if family == 'inlay':
        half *= 0.48
    if family == 'onlay':
        half *= 0.82
    if family == 'veneer':
        half *= 0.8
    base = base + R[:, :2] @ center
    (xy, faces, uv) = grid(half, family)
    target_labels = [34, 35, 36] if family == 'bridge3' else [fdi]
    lv = (arch['v'] - base) @ R
    tri = lv[arch['f']]
    target_tri = crop(tri[np.isin(arch['owner'], target_labels)], xy)
    native = height(target_tri, xy)
    if np.isfinite(native).mean() < 0.25:
        raise ValueError('too little native projected reference coverage')
    ni = ~np.isin(arch['labels'], target_labels + [0])
    neighbors = lv[ni]
    if len(neighbors) < 100:
        raise ValueError('missing adjacent geometry')
    dd = np.linalg.norm(neighbors[:, :2], axis=1)
    near = neighbors[dd < max(half) * 2 + 5]
    if len(near) < 30:
        near = neighbors[np.argsort(dd)[:200]]
    ztop = float(np.quantile(near[:, 2], 0.85))
    U = ((pair['upper']['v'] - base) @ R)[pair['upper']['f']]
    U = crop(U, xy)
    if family == 'veneer':
        U = np.empty((0, 3, 3))
    if len(U):
        (A, b) = constraints(U, xy, faces, 0.0)
    else:
        A = csr_matrix((0, len(xy)))
        b = np.empty(0)
    if A.shape[0] and (not np.isfinite(b).all() or np.min(A.data) < -1e-07 or np.max(np.abs(np.asarray(A.sum(1)).ravel() - 1)) > 1e-07):
        raise ValueError('invalid continuous obstacle rows')
    ceiling = height(U, xy, True)
    return dict(xy=xy, faces=faces, uv=uv, A=A, b=b, ceiling=ceiling, reference=native, target_tri=target_tri, antagonist=U, weights=area_weights(xy, faces), neighbor_points=near[np.linspace(0, len(near) - 1, min(len(near), 1500), dtype=int)], ztop=ztop, half=half, base=base, R=R, source_fdi=fdi)

def wall(family, level):
    if family == 'veneer':
        return (0.6, 'IPS_e_max_CAD_veneer')
    if family == 'inlay':
        return (1.0, 'IPS_e_max_CAD_inlay')
    if family == 'onlay':
        return (1.5, 'IPS_e_max_CAD_partial')
    if family == 'lattice_onlay':
        return (0.5, 'engineering_lattice_no_clinical_IFU')
    if family == 'bridge3':
        return (0.5 if level in ['easy', 'normal'] else 1.0, 'KATANA_HTML_PLUS' if level in ['easy', 'normal'] else 'KATANA_STML')
    anterior = family == 'anterior_crown'
    return ((0.4 if anterior else 0.5) if level in ['easy', 'normal'] else 1.0, 'KATANA_HTML_PLUS' if level in ['easy', 'normal'] else 'IPS_e_max_CAD_adhesive_1mm')

def run():
    if (ROOT / 'data/BENCHMARK_LOCK.json').exists():
        lock = read(ROOT / 'data/BENCHMARK_LOCK.json')
        for (p, h) in lock['payload']['files'].items():
            if sha(DATA / p) != h:
                raise ValueError('Locked benchmark data drift ' + p)
        return lock['payload']
    start = time.perf_counter()
    cohort = select()
    donors = {k: [] for k in FAMILIES}
    source_rows = []
    donor_errors = []
    for row in cohort:
        if row['split'] != 'train':
            continue
        (p, src) = pair(row)
        source_rows.append(dict(case_key=row['case_key'], sources=src))
        for family in FAMILIES:
            try:
                s = make_site(p, family)
                n = s['reference'] - s['ztop']
                ok = np.isfinite(n)
                if ok.sum() < 5:
                    raise ValueError('sparse donor')
                (_, j) = cKDTree(s['uv'][ok]).query(s['uv'][~ok])
                n[~ok] = n[ok][j]
                donors[family].append(n)
            except (ValueError, KeyError) as e:
                donor_errors.append(dict(case_key=row['case_key'], family=family, reason=str(e)))
        print('donor', row['case_key'], flush=True)
    template = {k: np.median(v, axis=0) for (k, v) in donors.items() if v}
    fit_s = time.perf_counter() - start
    dump(ROOT / 'raw/DONOR_FIT.json', dict(case_keys=[r['case_key'] for r in cohort if r['split'] == 'train'], counts={k: len(v) for (k, v) in donors.items()}, errors=donor_errors, seconds=fit_s))
    public = DATA / 'public'
    private = DATA / 'private'
    private.mkdir(parents=True, exist_ok=True)
    private.chmod(448)
    (public / 'scenes').mkdir(parents=True, exist_ok=True)
    tasks = []
    manifest = []
    for row in cohort:
        if row['split'] == 'train':
            continue
        state('BUILDING_TASKS', 'Donor-only shape fit complete', 'Native task ' + row['case_key'], tasks_built=len(tasks))
        (p, src) = pair(row)
        source_rows.append(dict(case_key=row['case_key'], sources=src))
        for (fi, family) in enumerate(FAMILIES):
            sid = row['case_key'] + '_' + family
            err = None
            try:
                s = make_site(p, family)
                if family not in template:
                    raise ValueError('no train template')
                prior = s['ztop'] + template[family]
                sc = public / 'scenes' / (sid + '.npz')
                A = s['A']
                np.savez_compressed(sc, xy=s['xy'], faces=s['faces'], uv=s['uv'], prior=prior, ceiling=s['ceiling'], weights=s['weights'], neighbor_points=s['neighbor_points'], antagonist=s['antagonist'], a_data=A.data, a_indices=A.indices, a_indptr=A.indptr, a_rows=A.shape[0], b=s['b'])
                ref = private / (sid + '.npz')
                np.savez_compressed(ref, reference=s['reference'], target_tri=s['target_tri'], base=s['base'], R=s['R'])
                if row['split'] == 'dev':
                    dev = DATA / 'development_reference'
                    dev.mkdir(exist_ok=True)
                    np.savez_compressed(dev / (sid + '.npz'), reference=s['reference'])
                prep0 = float(np.quantile(prior, 0.1))
                reference = dict(file=ref.name, sha256=sha(ref), coverage=float(np.isfinite(s['reference']).mean()))
            except (ValueError, KeyError) as e:
                err = str(e)
            for (li, level) in enumerate(LEVELS):
                tid = sid + '_' + level
                task = dict(task_id=tid, case_key=row['case_key'], dataset=row['dataset'], split=row['split'], family=family, level=level, source_fdi=FDI[family], units='mm', frame='site_' + sid, status='UNKNOWN_SITE' if err else 'READY', site_error=err)
                if not err:
                    (w, material) = wall(family, level)
                    reduction = [2.0, 1.5, 1.0, 1.0][li]
                    prep = prep0 - reduction
                    if level == 'boundary':
                        ceiling_min = float(s['b'].min()) if len(s['b']) else float(np.min(prior))
                        sign = 1 if int(row['case_key'][-1], 16) % 2 else -1
                        prep = ceiling_min - 0.02 - 0.04 - w - sign * 0.02
                    req = dict(wall_mm=w, film_min_mm=0.04, film_max_mm=0.12, clearance_mm=0.02, tool_radius_mm=[0.25, 0.5, 0.8, 1.0][li], connector_mm2=9.0 if material == 'KATANA_HTML_PLUS' else 16.0, channel_max_deg=25.0, strut_mm=0.35)
                    task.update(geometry_file=sc.name, geometry_sha256=sha(sc), material=material, requirements=req, preparation_height_mm=prep, preparation_rule=dict(convergence_total_deg=6.0, occlusal_reduction_scenario_mm=reduction, axial_reduction_nominal_mm=1.0, margin='rounded chamfer 0.5 mm specified for full-crown extension; not instantiated in roof-only geometry', actual_geometry='horizontal preparation patch; no cervical axial walls'), connector_x=[-float(s['half'][0]) * 0.31, float(s['half'][0]) * 0.31], implant_axis=[float(np.sin(np.radians([0, 10, 20, 25.02][li]))), 0, float(np.cos(np.radians([0, 10, 20, 25.02][li])))], minimum_strut_mm=float(min(s['half']) / 6), partial_restoration=True)
                    manifest.append(dict(task_id=tid, reference=reference, source_case=row['case_key']))
                dump(public / (tid + '.json'), task)
                tasks.append(task)
        print('tasks', row['case_key'], len(tasks), flush=True)
        budget()
    dump(private / 'references.json', manifest)
    dump(ROOT / 'raw/SOURCE_MEMBERS.json', source_rows)
    files = {str(p.relative_to(DATA)): sha(p) for d in [public, private, DATA / 'development_reference'] for p in sorted(d.rglob('*')) if p.is_file()}
    payload = dict(files=files, task_count=len(tasks), ready=sum((t['status'] == 'READY' for t in tasks)), unknown=sum((t['status'] != 'READY' for t in tasks)), counts={f: {l: sum((t['family'] == f and t['level'] == l for t in tasks)) for l in LEVELS} for f in FAMILIES}, source_records=len(cohort), case_counts={s: sum((r['split'] == s for r in cohort)) for s in ['train', 'dev', 'test']}, fit_seconds=fit_s, preparation_seconds=time.perf_counter() - start, bytes=budget())
    freeze(ROOT / 'data/BENCHMARK_LOCK.json', payload)
    return payload
