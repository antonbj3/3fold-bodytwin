from common import *
import csv, collections, hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def icc(values):
    x = np.asarray(values, float)
    if x.ndim != 2 or x.shape[0] < 3 or x.shape[1] not in [2, 3] or (not np.isfinite(x).all()) or np.any((x < 1) | (x > 5)):
        raise ValueError('ICC needs complete crossed >=3 objects, 2-3 raters, scores1..5')
    (n, k) = x.shape
    grand = x.mean()
    MSR = k * np.sum((x.mean(1) - grand) ** 2) / (n - 1)
    MSC = n * np.sum((x.mean(0) - grand) ** 2) / (k - 1)
    MSE = np.sum((x - x.mean(1)[:, None] - x.mean(0)[None, :] + grand) ** 2) / ((n - 1) * (k - 1))
    d1 = MSR + (k - 1) * MSE + k * (MSC - MSE) / n
    dk = MSR + (MSC - MSE) / n
    return dict(ICC_2_1=float((MSR - MSE) / d1) if d1 else None, ICC_2_k=float((MSR - MSE) / dk) if dk else None, objects=n, raters=k, MSR=MSR, MSC=MSC, MSE=MSE, scope='Two-way random absolute agreement; ordinal1..5 ratings treated as interval; negative ICC retained')

def analyze(path):
    rr = list(csv.DictReader(open(path)))
    by = collections.defaultdict(dict)
    raters = set()
    for r in rr:
        if not r.get('score'):
            continue
        key = (r['image_id'], r['dimension'])
        rid = r['rater_id']
        if rid in by[key]:
            raise ValueError('Duplicate image/dimension/rater')
        by[key][rid] = float(r['score'])
        raters.add(rid)
    if not by:
        return dict(status='UNKNOWN_NO_RATINGS', ICC=None)
    out = {}
    rids = sorted(raters)
    for dim in sorted({k[1] for k in by}):
        objects = [v for (k, v) in by.items() if k[1] == dim]
        complete = [v for v in objects if set(v) == raters]
        if len(complete) < 3 or len(rids) not in [2, 3]:
            out[dim] = dict(status='UNKNOWN_INCOMPLETE_CROSSED_DESIGN', complete=len(complete), requested=len(objects))
            continue
        arr = np.array([[v[r] for r in rids] for v in complete])
        q = icc(arr)
        rng = np.random.default_rng(6105)
        vals = []
        for _ in range(1000):
            z = icc(arr[rng.integers(0, len(arr), len(arr))])
            vals.append(z['ICC_2_1'])
        finite = [v for v in vals if v is not None]
        q.update(bootstrap_object_CI95=np.quantile(finite, [0.025, 0.975]) if finite else None, interval_scope='Object bootstrap descriptive only; designs share patient; not independent patient uncertainty', incomplete_excluded=len(objects) - len(complete))
        out[dim] = q
    return dict(status='DESCRIPTIVE_RATINGS', dimensions=out)

def make_pack():
    from functional import sample
    dest = ROOT / 'blind_pack'
    dest.mkdir(exist_ok=True)
    frozen = ROOT / 'FROZEN_IMAGE_PACKAGE.json'
    if frozen.exists():
        manifest = read(frozen)
        if manifest.get('renderer_sha256') != sha(Path(__file__)):
            raise ValueError('Frozen renderer changed: preserve previous package before versioning')
        if manifest.get('private_key_sha256') != sha(ROOT / 'raw/BLIND_KEY_PRIVATE.json'):
            raise ValueError('Frozen blind key drift')
        actual = {p.name: sha(p) for p in dest.iterdir() if p.is_file()}
        if actual != manifest['files']:
            raise ValueError('Frozen image package changed: preserve previous package before versioning')
        return dict(images=manifest['objects'], requested_objects=33, unrenderable_generation_failures=33 - manifest['objects'], raters_planned=[2, 3], ICC_status='UNKNOWN_NO_RATINGS', images_mode='FROZEN_BYTES_VERIFIED', analysis=analyze(dest / 'RATINGS_TEMPLATE.csv'))
    keys = [r['key'] for r in read(ROOT / 'EXTERNAL_INPUTS.json')['rows']]
    records = []
    for key in keys:
        for folder in ['whole_predictions', 'prep_predictions']:
            for name in sorted((V4 / 'payload' / folder).iterdir()):
                p = name / key / 'mesh.npz'
                if p.exists():
                    records.append(dict(key=key, participant=name.name, track=folder, path=p))
        p = DATA / 'external_meshes' / key / 'mesh.npz'
        if p.exists():
            records.append(dict(key=key, participant='ToothCraft_normal', track='external', path=p))
    bounds = {}
    for key in keys:
        pp = npz(DATA / 'external_inputs' / f'{key}.npz')
        cloud = []
        for r in records:
            if r['key'] == key:
                cloud.append(npz(r['path'])['vertices'])
        v = np.vstack(cloud)
        lo = v.min(0)
        hi = v.max(0)
        bounds[key] = ((lo + hi) / 2, float(np.max(hi - lo) / 2 * 1.12))
    rng = np.random.default_rng(6105)
    rng.shuffle(records)
    private = []
    template = []
    for (index, r) in enumerate(records):
        ident = 'B' + hashlib.sha256(('v5-blind-seed6105/' + str(index)).encode()).hexdigest()[:10]
        m = npz(r['path'])
        p = npz(DATA / 'external_inputs' / f"{r['key']}.npz")
        tri = m['vertices'][m['faces']]
        ctx = np.concatenate([p['neighbor_triangles'], npz(V4 / 'payload/whole_private' / r['key'] / 'reference.npz')['antagonist_triangles']])
        prep = p['preparation_vertices'][p['preparation_faces']]
        fig = plt.figure(figsize=(10, 4))
        (center, rad) = bounds[r['key']]
        ctx = ctx[np.all((ctx.mean(1) >= center - rad) & (ctx.mean(1) <= center + rad), axis=1)]
        prep = prep[np.all((prep.mean(1) >= center - rad) & (prep.mean(1) <= center + rad), axis=1)]
        for (ii, (el, az)) in enumerate([(90, -90), (15, -90), (15, 0)]):
            ax = fig.add_subplot(1, 3, ii + 1, projection='3d')
            ax.add_collection3d(Poly3DCollection(tri, facecolors='#64b4d8', linewidths=0, alpha=1, shade=True, lightsource=LightSource(315, 45)))
            ax.add_collection3d(Poly3DCollection(ctx, facecolor='#999999', edgecolor='none', alpha=0.13))
            ax.add_collection3d(Poly3DCollection(prep, facecolor='#e7ab65', edgecolor='none', alpha=0.3))
            for (k, fn) in enumerate([ax.set_xlim, ax.set_ylim, ax.set_zlim]):
                fn(center[k] - rad, center[k] + rad)
            ax.set_box_aspect((1, 1, 1))
            ax.set_proj_type('ortho')
            ax.view_init(el, az)
            ax.set_axis_off()
            ax.set_title(['Occlusal / local Z', 'Side A / local -Y', 'Side B / local X'][ii])
            q = center.copy()
            q[2] = center[2] - rad * 0.85
            q[0] = center[0] - rad * 0.65
            q[1] = center[1] - rad * 0.65
            end = q.copy()
            end[1 if ii == 2 else 0] += 2
            ax.plot(*np.stack([q, end]).T, color='black', lw=2)
            ax.text(*(q + end) / 2, '2 mm', fontsize=8)
        fig.suptitle(ident)
        fig.text(0.5, 0.025, 'Blue: candidate | grey: measured opposing/neighbour surfaces | orange: virtual preparation', ha='center', fontsize=8)
        fig.tight_layout(rect=[0, 0.05, 1, 0.95])
        file = dest / (ident + '.png')
        fig.savefig(file, dpi=130)
        plt.close(fig)
        private.append(dict(image_id=ident, key=r['key'], participant=r['participant'], track=r['track'], mesh_sha256=sha(r['path']), image_sha256=sha(file), file=file.name))
        for rid in ['R1', 'R2', 'R3']:
            for dim in ['occlusal_morphology', 'proximal_relation', 'marginal_form', 'overall_visual_plausibility']:
                template.append(dict(image_id=ident, rater_id=rid, dimension=dim, score='', unable_to_assess='', comment=''))
    with (dest / 'RATINGS_TEMPLATE.csv').open('w') as f:
        w = csv.DictWriter(f, list(template[0]))
        w.writeheader()
        w.writerows(template)
    (dest / 'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Blinded crown assessment</title><style>body{font:18px sans-serif;max-width:1100px;margin:auto}img{width:100%}article{border-bottom:1px solid #ccc;padding:2em 0}</style><h1>Blinded crown images</h1><p>Blue: candidate. Grey: measured neighbours and antagonist. Orange: virtual preparation. Use the supplied rating protocol; these images do not establish physical fit.</p>' + ''.join(('<article><h2>' + r['image_id'] + '</h2><img src="' + r['file'] + '" alt="' + r['image_id'] + '"></article>' for r in private)))
    dump(ROOT / 'raw/BLIND_KEY_PRIVATE.json', private)
    manifest = {p.name: sha(p) for p in dest.iterdir() if p.is_file()}
    f = ROOT / 'FROZEN_IMAGE_PACKAGE.json'
    if not f.exists():
        freeze(f, dict(files=manifest, objects=len(records), independent_patient_cases=1, requested_objects=33, unrenderable_generation_failures=33 - len(records), anatomic_direction_independently_verified=False, renderer_sha256=sha(Path(__file__)), private_key_sha256=sha(ROOT / 'raw/BLIND_KEY_PRIVATE.json'), rating_status='NOT_RUN', key='raw/BLIND_KEY_PRIVATE.json: withhold from raters', note='First three tooth types from one case; feasibility pilot, no population agreement claim'))
    elif read(f)['files'] != manifest:
        raise ValueError('Frozen image package drift: keep previous and create new version')
    return dict(images=len(records), requested_objects=33, unrenderable_generation_failures=33 - len(records), images_mode='RENDERED_AND_FROZEN', raters_planned=[2, 3], ICC_status='UNKNOWN_NO_RATINGS', analysis=analyze(dest / 'RATINGS_TEMPLATE.csv'))
if __name__ == '__main__':
    if len(sys.argv) > 1:
        dump(ROOT / 'LAB_RATINGS_ANALYSIS.json', analyze(sys.argv[1]))
    else:
        dump(ROOT / 'raw/LAB_PACKAGE.json', make_pack())
