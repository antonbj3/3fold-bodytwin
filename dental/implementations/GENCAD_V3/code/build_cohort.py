"""Author-only local data preparation. Evaluator never invokes this after release."""
from dental_release.paths import expand as _release_expand
import sys, hashlib, time, collections, resource
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'legacy/vendor'))
from util import *
from source_ops import pair, make_site, wall, FDI
from scipy.spatial import cKDTree
from packing import reduce_equal_rows, savez
LABEL = Path(_release_expand('@DENTAL_WORK_ROOT@/X11/targets'))

def inventory():
    lock = PAYLOAD / 'private/COHORT.json'
    if lock.exists():
        return read(lock)['payload']['cases']
    old = read(ROOT.parent / 'PROOF_LANE_GENCAD_V2/data/COHORT_LOCK.json')['payload']['cases']
    legacy = {(r['dataset'], r['case']) for r in old}
    rows = []
    excluded = []
    parents = {}

    def find(x):
        parents.setdefault(x, x)
        if parents[x] != x:
            parents[x] = find(parents[x])
        return parents[x]
    seen = {}
    for ds in ['Bite2Text', 'Bits2Bites']:
        for p in sorted((LABEL / ds).iterdir()):
            if not p.is_dir():
                continue
            m = read(p / 'labels+landmarks.json')
            key = digest(ds + '|' + p.name)[:16]
            if m.get('input_status') == 'REFUSED_EMPTY_SOURCE':
                excluded.append(dict(dataset=ds, source_record=p.name, reason='REFUSED_EMPTY_SOURCE'))
                continue
            hs = [m['arches'][j]['input_sha256'] for j in ['upper', 'lower']]
            r = dict(dataset=ds, case=p.name, case_key=key, metadata=str(p / 'labels+landmarks.json'), metadata_sha256=sha(p / 'labels+landmarks.json'), geometry_hashes=hs, legacy_v2=(ds, p.name) in legacy)
            rows.append(r)
            find(key)
            for h in hs:
                if h in seen:
                    parents[find(key)] = find(seen[h])
                else:
                    seen[h] = key
    for r in rows:
        r['group'] = find(r['case_key'])
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r['group']].append(r)
    old_groups = {g for (g, rs) in groups.items() if any((r['legacy_v2'] for r in rs))}
    primary = sorted([g for (g, rs) in groups.items() if any((r['dataset'] == 'Bite2Text' for r in rs)) and g not in old_groups], key=lambda g: digest('PROOF_LANE-v3-20261002|' + g))
    split = {g: 'train' if i < int(0.2 * len(primary)) else 'dev' if i < int(0.4 * len(primary)) else 'test' for (i, g) in enumerate(primary)}
    split.update({g: 'dev' for g in old_groups})
    keep = []
    primary_seen = set()
    for r in rows:
        r['split'] = split.get(r['group'], 'auxiliary') if r['dataset'] == 'Bite2Text' else 'auxiliary'
        if r['dataset'] == 'Bite2Text' and r['group'] in primary_seen:
            excluded.append(dict(dataset=r['dataset'], source_record=r['case'], reason='duplicate_arch_patient_group', group=r['group']))
            continue
        if r['dataset'] == 'Bite2Text':
            primary_seen.add(r['group'])
        keep.append(r)
    keep.sort(key=lambda r: (['train', 'dev', 'test', 'auxiliary'].index(r['split']), digest('PROOF_LANE-v3-20261002|' + r['group'])))
    out = dict(cases=keep, excluded=excluded, attempted=len(rows) + sum((e['reason'] == 'REFUSED_EMPTY_SOURCE' for e in excluded)), retained=len(keep), patient_identity='Within Bite2Text case/patient per upstream documentation. Cross-dataset and Teeth3DS-pretraining overlap UNKNOWN.', source_locator='https://ditto.ing.unimore.it/bite2text/', counts=dict(collections.Counter((r['dataset'] + '/' + r['split'] for r in keep))))
    freeze(lock, out)
    dump(ROOT / 'raw/COHORT_SUMMARY.json', {k: v for (k, v) in out.items() if k != 'cases'})
    selected = []
    for sp in ['dev', 'test']:
        selected.extend([dict(case_key=r['case_key'], split=sp, dataset=r['dataset'], group=r['group']) for r in keep if r['split'] == sp][:12])
    freeze(PAYLOAD / 'private/PREP_SELECTION.json', selected)
    return keep

def fit(cohort):
    fn = PAYLOAD / 'public/templates.npz'
    if fn.exists():
        with np.load(fn, allow_pickle=False) as a:
            return dict(a)
    donors = {k: [] for k in FAMILIES}
    err = []
    t0 = time.perf_counter()
    for (ii, row) in enumerate((r for r in cohort if r['split'] == 'train')):
        (p, sources) = pair(row)
        for family in FAMILIES:
            try:
                s = make_site(p, family, compute_obstacles=False)
                n = s['reference'] - s['ztop']
                ok = np.isfinite(n)
                if ok.sum() < 5:
                    raise ValueError('sparse template')
                if not ok.all():
                    (_, j) = cKDTree(s['uv'][ok]).query(s['uv'][~ok])
                    n[~ok] = n[ok][j]
                donors[family].append(n)
            except (ValueError, KeyError) as e:
                err.append(dict(case_key=row['case_key'], family=family, reason=str(e)))
        if ii % 10 == 0:
            print('donor', ii + 1, flush=True)
    template = {k: np.median(v, axis=0) for (k, v) in donors.items() if v}
    fn.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(fn, **template)
    dump(PAYLOAD / 'private/TEMPLATE_FIT.json', dict(groups=[r['group'] for r in cohort if r['split'] == 'train'], counts={k: len(v) for (k, v) in donors.items()}, rejections=err, seconds=time.perf_counter() - t0))
    return template

def build_case(row, template):
    dest = PAYLOAD / 'public/tasks' / (row['case_key'] + '.json')
    if dest.exists():
        return
    (p, sources) = pair(row)
    tasks = []
    refs = {}
    errors = []
    for family in FAMILIES:
        sid = row['case_key'] + '_' + family
        err = None
        try:
            s = make_site(p, family)
            if family not in template:
                raise ValueError('missing train template')
            prior = s['ztop'] + template[family]
            (A, b) = reduce_equal_rows(s['A'], s['b'])
            sc = PAYLOAD / 'public/scenes' / (sid + '.npz')
            sc.parent.mkdir(parents=True, exist_ok=True)
            savez(sc, xy=s['xy'], faces=s['faces'], uv=s['uv'], prior=prior, ceiling=s['ceiling'], weights=s['weights'], a_data=A.data, a_indices=A.indices, a_indptr=A.indptr, a_rows=A.shape[0], b=b)
            refs[family] = s['reference']
            prep0 = float(np.quantile(prior, 0.1))
        except (ValueError, KeyError) as e:
            err = str(e)
            errors.append(dict(family=family, reason=err))
        for (li, level) in enumerate(LEVELS):
            tid = sid + '_' + level
            t = dict(task_id=tid, case_key=row['case_key'], split=row['split'], family=family, level=level, units='mm', frame='site_' + sid, status='UNKNOWN_SITE' if err else 'READY', site_error=err, resolution='PER_POINT')
            if not err:
                (w, material) = wall(family, level)
                reduction = [2.0, 1.5, 1.0, 1.0][li]
                prep = prep0 - reduction
                if level == 'boundary':
                    ceiling_min = float(s['b'].min()) if len(s['b']) else float(prior.min())
                    sign = 1 if int(row['case_key'][-1], 16) % 2 else -1
                    prep = ceiling_min - 0.02 - 0.04 - w - sign * 0.02
                req = dict(wall_mm=w, film_min_mm=0.04, film_max_mm=0.12, clearance_mm=0.02, connector_mm2=9.0 if material == 'KATANA_HTML_PLUS' else 16.0, channel_max_deg=25.0, strut_mm=0.35)
                t.update(geometry_file='scenes/' + sc.name, requirements=req, preparation_height_mm=prep, material=material, connector_x=[-float(s['half'][0]) * 0.31, float(s['half'][0]) * 0.31], implant_axis=[float(np.sin(np.radians([0, 10, 20, 25.02][li]))), 0.0, float(np.cos(np.radians([0, 10, 20, 25.02][li])))], minimum_strut_mm=float(min(s['half']) / 6), scope='Restricted affine roof; axial preparation is evaluated separately in anatomical track; milling NOT_SCORED')
            tasks.append(t)
    rdir = PAYLOAD / 'private/references'
    rdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(rdir / (row['case_key'] + '.npz'), **refs)
    dump(PAYLOAD / 'private/provenance' / (row['case_key'] + '.json'), dict(case_key=row['case_key'], sources=sources, errors=errors))
    dump(dest, tasks)

def main():
    start = time.perf_counter()
    cohort = inventory()
    state('BUILDING_COHORT', 'Frozen splits and prep selection', 'Fit primary train only, stream registered scans', records=len(cohort))
    template = fit(cohort)
    for (i, row) in enumerate((r for r in cohort if r['split'] != 'train')):
        build_case(row, template)
        if i % 5 == 0:
            print('case', i + 1, row['split'], row['case_key'], 'seconds', round(time.perf_counter() - start), flush=True)
            state('BUILDING_COHORT', 'Case checkpoint saved', 'Next frozen source record', case_index=i + 1, data_bytes=budget())
    dump(ROOT / 'raw/BUILD_COST.json', dict(wall_seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, data_bytes=budget()))
if __name__ == '__main__':
    main()
