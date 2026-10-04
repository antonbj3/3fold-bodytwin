"""Executed positive controls and review regressions; no clinical interpretation."""
from dental_release.paths import expand as _release_expand
import copy, tempfile, shutil, os, sys, json, time, subprocess, io, zipfile, warnings
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from util import *
from integrity import Bundle, IntegrityError, strict_json, safe_bytes, relative, h, predictions, load_npz
from legacy.checks import all_checks, feasibility
from legacy.generators import submit, optimize
from legacy.geometry import shell, section_weights
from task_io import load_local
from case_statistics import cluster_summary

def fixture(family='molar_crown'):
    xy = np.array([[-2.0, -2.0], [2.0, -2.0], [-2.0, 2.0], [2.0, 2.0]])
    f = np.array([[0, 1, 2], [1, 3, 2]])
    return dict(task_id='fixture', frame='fixed', family=family, xy=xy, faces=f, A=csr_matrix(np.eye(4)), obstacle_b=np.ones(4) * 2.0, preparation_z=np.zeros(4), requirements=dict(wall_mm=1.0, film_min_mm=0.04, film_max_mm=0.12, clearance_mm=0.02, connector_mm2=4.0, channel_max_deg=25.0, strut_mm=0.35), connector_x=[-0.5, 0.5], implant_axis=[0, 0, 1], minimum_strut_mm=0.5, prior=np.ones(4) * 1.2, weights=np.ones(4))

def geometry_controls():
    rows = []

    def check(name, positive, reject):
        row = dict(name=name, positive=bool(positive), injected_fault_rejected=bool(reject))
        rows.append(row)
        if not all([positive, reject]):
            raise AssertionError(name)
    t = fixture()
    d = submit(t, np.ones(4) * 1.2)
    good = all_checks(t, d)['validity'] == 'PASS'
    for (k, v) in [('units', 'm'), ('frame', 'forged'), ('outer_vertices', [[0, 0]] * 4), ('inner_vertices', [[0, 0, 0, 1]] * 4), ('faces', [[0, 1, -1], [1, 3, 2]]), ('faces', [[0, 1, 2]]), ('reported_wall_mm', 999)]:
        bad = copy.deepcopy(d)
        bad[k] = v
        check('contract:' + k + ':' + str(v)[:20], good, all_checks(t, bad)['validity'] == 'INVALID')
    for x in [float('nan'), float('inf')]:
        bad = copy.deepcopy(d)
        bad['outer_vertices'][0][2] = x
        check('nonfinite ' + str(x), good, all_checks(t, bad)['validity'] == 'INVALID')
    bad = copy.deepcopy(d)
    bad['outer_vertices'][0][0] += 1
    check('changed footprint', good, all_checks(t, bad)['validity'] == 'INVALID')
    for (z, inner, key) in [(0.9, 0.08, 'wall'), (-1, 0.08, 'nesting'), (1.2, 0.01, 'film_min'), (1.5, 0.2, 'film_max'), (3.0, 0.08, 'antagonist')]:
        b = submit(t, np.full(4, z), np.full(4, inner))
        check(key, good, all_checks(t, b)['checks'][key]['status'] == 'FAIL')
    for (fam, key, change) in [('bridge3', 'connector_area', lambda t: t['requirements'].update(connector_mm2=5)), ('implant_crown', 'channel_angle', lambda t: t.update(implant_axis=[0.6, 0, 0.8])), ('lattice_onlay', 'strut_width', lambda t: t.update(minimum_strut_mm=0.2))]:
        tf = fixture(fam)
        df = submit(tf, np.ones(4) * 1.2)
        positive = all_checks(tf, df)['validity'] == 'PASS'
        change(tf)
        check(key, positive, all_checks(tf, df)['checks'][key]['status'] == 'FAIL')
    tv = fixture()
    positive = feasibility(tv)['status'] == 'FEASIBLE'
    tv['obstacle_b'] = np.full(4, 0.9)
    check('necessary infeasibility', positive, feasibility(tv)['status'] == 'INFEASIBLE' and optimize(tv)['status'] == 'ABSTAIN')
    tv = fixture()
    tv['A'] = csr_matrix((0, 4))
    tv['obstacle_b'] = []
    check('no antagonist cannot PASS', good, all_checks(tv, d)['validity'] == 'UNKNOWN')
    w = section_weights(t['xy'], t['faces'], 0.25)
    check('section closed form', abs(w @ np.full(4, 2.0) - 8) < 1e-09, abs(w @ np.full(4, 3.0) - 8) > 1)
    from preparation import bad_edges
    (v, f) = shell(t['xy'], np.ones(4) * 1.2, np.ones(4) * 0.08, t['faces'])
    check('closed shell', bad_edges(f) == 0, bad_edges(f[:-1]) > 0)
    from legacy_v1.checks import cement, milling, cement_max
    planes = [[1, 0, 0, 2], [-1, 0, 0, 2], [0, 1, 0, 2], [0, -1, 0, 2], [0, 0, 1, 2], [0, 0, -1, 2]]
    pl = [[0, 0, 1, 0]]
    check('v1 two-dimensional preparation', cement.check([[0, 0, 0]], planes, 0.1)['status'] == 'PASS', cement.check([[0, 0]], planes, 0.1)['status'] == 'INVALID')
    check('v1 two-dimensional milling centre', milling.verify_positive(pl, [[0, 0, 0]], 1, 0, [0, 0, -1], [[0, 0, -1]]), not milling.verify_positive(pl, [[0, 0, 0]], 1, 0, [0, 0, -1], [[0, 0]]))
    check('v1 four-dimensional insertion direction', milling.verify_positive(pl, [[0, 0, 0]], 1, 0, [0, 0, -1], [[0, 0, -1]]), not milling.verify_positive(pl, [[0, 0, 0]], 1, 0, [0, 0, 0, 1], [[0, 0, -1]]))
    check('v1 empty cement witness', cement_max.verify_positive(planes, [[0, 0, 0]], [[0, 0, 0]], 0.12), not cement_max.verify_positive(planes, [[0, 0, 100]], [], 0.12))
    check('milling intrusion', milling.verify_positive(pl, [[0, 0, 0]], 0.5, 0, [0, 0, -1], [[0, 0, -0.5]]), not milling.verify_positive(pl, [[0, 0, 0]], 0.5, 0, [0, 0, -1], [[0, 0, 0]]))
    m = np.repeat(np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])[:, None], 4, axis=1)
    a = cluster_summary(m)
    b = cluster_summary(np.repeat(m, 100, axis=1))
    bad = cluster_summary(m.reshape(-1, 1))
    check('cluster replication cannot narrow CI', a['case_cluster_bootstrap_95'] == b['case_cluster_bootstrap_95'], bad['case_cluster_bootstrap_95'][1] - bad['case_cluster_bootstrap_95'][0] < a['case_cluster_bootstrap_95'][1] - a['case_cluster_bootstrap_95'][0])
    from topology import split_fans
    from source_solid import cap_patch
    from source_solid_replay import triangles_canonical
    v = np.array([[0.0, 0, 1], [1, 0, 1], [0, 1, 1], [-1, 0, 1], [0, -1, 1]])
    f = np.array([[0, 1, 2], [0, 3, 4]])
    (vv, ff, info) = split_fans(v, f)
    (closed_v, closed_f, loops) = cap_patch(vv, ff, -1.0)
    caught = False
    try:
        cap_patch(v, f, -1.0)
    except ValueError:
        caught = True
    check('branched source vertex rejected before split', bad_edges(closed_f) == 0 and len(loops) == 2, caught)
    before = triangles_canonical(v[f])
    after = triangles_canonical(vv[ff])
    changed = after.copy()
    changed[0, 0, 0] += 0.001
    check('source triangle identity', np.array_equal(before, after), not np.array_equal(before, changed))
    return rows

def atomic(path, blob):
    p = Path(path)
    temp = p.with_name(p.name + '.probe_tmp')
    temp.write_bytes(blob)
    temp.replace(p)

def linked_copy(src, dst):
    """Own immutable release only. All mutations use replace, never write a link."""
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)
    return dst

def lp_controls(bundle):
    from fast_geometry import optimize as fast_optimize, feasibility as fast_feasibility
    from task_io import attach
    cases = [(p, ts) for (p, ts) in bundle.tasks() if ts[0]['split'] == 'dev'][:6]
    rows = []
    fault = False
    for (path, ts) in cases:
        cached = {}
        for t0 in ts:
            if t0['status'] != 'READY':
                continue
            fn = t0['geometry_file']
            if fn not in cached:
                cached[fn] = bundle.npz('payload/public/' + fn)
            t = attach(t0, cached[fn])
            a = fast_optimize(t)
            b = optimize(t)
            err = None
            if a['status'] == b['status'] == 'DESIGN':
                z = np.asarray(a['outer_vertices'])[:, 2]
                ref = np.asarray(b['outer_vertices'])[:, 2]
                objective = lambda x: float(np.sum(t['weights'] * abs(x - t['prior'])))
                err = abs(objective(z) - objective(ref))
                fault = fault or abs(objective(ref + 1.0) - objective(ref)) > 1e-08
            passed = a['status'] == b['status'] and (err is None or err <= 1e-08) and (fast_feasibility(t)['status'] == feasibility(t)['status'])
            rows.append(dict(task_id=t['task_id'], passed=passed, objective_abs_error_mm3=err, resolution='PER_TOOTH'))
    if not rows or not all((r['passed'] for r in rows)) or (not fault):
        raise AssertionError('actual inherited LP parity/injection failed')
    return dict(count=len(rows), all_pass=True, altered_solution_rejected=fault, rows=rows, claim='Engineering equivalence check; no superiority claim')

def integrity_controls(bundle):
    from scorer import score, score_submission, decode, design
    from release_anchor import BENCHMARK_SHA256, PREDICTIONS_SHA256
    t0 = time.perf_counter()
    rows = []
    root = bundle.root
    scratch = bundle.payload.parent / 'integrity_probe'
    if scratch.exists():
        raise RuntimeError('old probe remains; inspect before rerun')
    scratch.mkdir()
    try:
        for name in ['BENCHMARK_LOCK.json', 'SCORER_PARAMETERS.json', 'PREREG_R1.json', 'PREREG_R2.json', 'DECOMPOSITION.json', 'RUNTIME_LOCK.json', 'SOURCE_REUSE.json', 'FROZEN_PREPARATIONS.json']:
            linked_copy(root / name, scratch / name)
        for folder in ['code', 'vendor', 'raw']:
            (scratch / folder).mkdir(exist_ok=True)
        for name in bundle.files:
            if name.startswith('payload/'):
                continue
            dest = scratch / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            if not dest.exists():
                linked_copy(root / name, dest)
        linked_copy(root / 'code/release_anchor.py', scratch / 'code/release_anchor.py')
        (scratch / 'payload').mkdir()
        folders = [p[len('payload/'):] for p in bundle.lock['closed_directories'] if p.startswith('payload/')] + ['predictions']
        for folder in folders:
            shutil.copytree(bundle.payload / folder, scratch / 'payload' / folder, copy_function=linked_copy)
        baseline = Bundle(scratch)
        witness = None
        for (path, ts) in baseline.tasks():
            if ts[0]['split'] != 'dev':
                continue
            for (i, t) in enumerate(ts):
                if t['status'] != 'READY' or t['family'] in ['bridge3', 'implant_crown', 'lattice_onlay']:
                    continue
                from legacy.generators import parametric
                loaded = load_local(t, scratch / 'payload/public')
                d = parametric(loaded)
                old = all_checks(loaded, d)
                if old['validity'] != 'FAIL':
                    continue
                changed = copy.deepcopy(t)
                changed['requirements']['clearance_mm'] = -1.0
                new = all_checks(load_local(changed, scratch / 'payload/public'), d)
                if new['validity'] == 'PASS':
                    witness = (path, ts, i, old['validity'], new['validity'])
                    break
            if witness:
                break
        if witness is None:
            raise AssertionError('no native C02 regression witness')
        (taskpath, ts, idx, old, new) = witness
        changed = copy.deepcopy(ts)
        changed[idx]['requirements']['clearance_mm'] = -1.0
        first = next((t for t in ts if t['status'] == 'READY'))
        scene = 'payload/public/' + first['geometry_file']
        case = Path(taskpath).stem
        ref = 'payload/private/references/' + case + '.npz'

        def mutate_file(label, rel, blob, entry='verify'):
            p = scratch / rel
            original = p.read_bytes()
            atomic(p, blob)
            caught = False
            reason = ''
            try:
                if entry == 'score':
                    score(root=scratch, output=scratch / 'ignored')
                elif entry == 'submit':
                    score_submission(scratch, scratch / 'payload/predictions', PREDICTIONS_SHA256, scratch / 'ignored')
                elif entry == 'read_after_verify':
                    baseline.json(rel)
                else:
                    Bundle(scratch)
            except IntegrityError as e:
                caught = True
                reason = str(e)
            finally:
                atomic(p, original)
            rows.append(dict(name=label, entry=entry, positive_restored=True, rejected=caught, reason=reason))
            if not caught:
                raise AssertionError('accepted mutation ' + label)
        blob = json.dumps(changed).encode()
        for entry in ['verify', 'score', 'submit', 'read_after_verify']:
            mutate_file('C02 clearance minus1', taskpath, blob, entry)
        for (key, value) in [('wall_mm', 0.0), ('film_min_mm', -0.5), ('film_max_mm', 50.0), ('connector_mm2', 0.0), ('channel_max_deg', 180.0)]:
            bad = copy.deepcopy(ts)
            bad[idx]['requirements'][key] = value
            mutate_file('task parameter ' + key, taskpath, json.dumps(bad).encode())
        bad = copy.deepcopy(ts)
        bad[idx]['units'] = 'm'
        mutate_file('task units', taskpath, json.dumps(bad).encode())
        bad = copy.deepcopy(ts)
        bad[idx]['geometry_file'] = '../../private/references/' + case + '.npz'
        mutate_file('task traversal', taskpath, json.dumps(bad).encode())
        bad = copy.deepcopy(ts)
        bad[idx]['preparation_height_mm'] -= 1.0
        mutate_file('preparation height', taskpath, json.dumps(bad).encode())
        for (name, rel) in [('antagonist geometry', scene), ('reference facit', ref), ('scorer code', 'code/scorer.py'), ('scorer parameter', 'SCORER_PARAMETERS.json'), ('interpreter', 'payload/runtime/bin/python3')]:
            original = (scratch / rel).read_bytes()
            mutate_file(name, rel, original + b'INJECTED')
        lock = read(scratch / 'BENCHMARK_LOCK.json')
        bad = copy.deepcopy(lock)
        bad['files'][taskpath]['sha256'] = h(blob)
        mutate_file('self-rehashed benchmark manifest', 'BENCHMARK_LOCK.json', json.dumps(bad).encode())
        bad = copy.deepcopy(lock)
        bad['schema'] = 'other-benchmark'
        mutate_file('wrong benchmark identity', 'BENCHMARK_LOCK.json', json.dumps(bad).encode())
        for (label, action) in [('extra task', lambda p: atomic(p, b'[]')), ('task symlink', lambda p: p.symlink_to(scratch / taskpath))]:
            p = scratch / 'payload/public/tasks/extra.json'
            action(p)
            try:
                caught = False
                try:
                    Bundle(scratch)
                except IntegrityError as e:
                    caught = True
                    reason = str(e)
                rows.append(dict(name=label, rejected=caught, positive_restored=True, reason=reason))
                assert caught
            finally:
                p.unlink()
        p = scratch / taskpath
        original = p.read_bytes()
        p.unlink()
        try:
            caught = False
            try:
                Bundle(scratch)
            except IntegrityError as e:
                caught = True
                reason = str(e)
            rows.append(dict(name='missing task', rejected=caught, positive_restored=True, reason=reason))
            assert caught
        finally:
            atomic(p, original)
        p = scratch / 'payload/public/tasks'
        hold = scratch / 'held_tasks'
        p.rename(hold)
        p.symlink_to(hold, target_is_directory=True)
        try:
            caught = False
            try:
                Bundle(scratch)
            except IntegrityError as e:
                caught = True
                reason = str(e)
            rows.append(dict(name='parent directory symlink', rejected=caught, positive_restored=True, reason=reason))
            assert caught
        finally:
            p.unlink()
            hold.rename(p)
        cases = [Path(p).stem for (p, _) in baseline.tasks()]
        proot = scratch / 'payload/predictions'
        pm = read(proot / 'FROZEN_PREDICTIONS.json')
        predictions(proot, PREDICTIONS_SHA256, cases, BENCHMARK_SHA256)
        rel = next(iter(pm['files']))
        p = proot / rel
        original = p.read_bytes()
        for action in ['changed', 'missing', 'symlink', 'extra']:
            extra = p.with_name('extra.npz')
            if action == 'changed':
                atomic(p, original + b'changed')
            elif action == 'missing':
                p.unlink()
            elif action == 'symlink':
                p.unlink()
                p.symlink_to(bundle.payload / 'predictions' / rel)
            else:
                atomic(extra, original)
            caught = False
            try:
                predictions(proot, PREDICTIONS_SHA256, cases, BENCHMARK_SHA256)
            except IntegrityError as e:
                caught = True
                reason = str(e)
            finally:
                if p.is_symlink():
                    p.unlink()
                atomic(p, original)
                if extra.exists():
                    extra.unlink()
            rows.append(dict(name='prediction ' + action, rejected=caught, positive_restored=True, reason=reason))
            assert caught
        p = proot / 'FROZEN_PREDICTIONS.json'
        orig = p.read_bytes()
        bad = copy.deepcopy(pm)
        bad['benchmark_sha256'] = '0' * 64
        atomic(p, json.dumps(bad).encode())
        try:
            caught = False
            try:
                predictions(proot, sha(p), cases, BENCHMARK_SHA256)
            except IntegrityError as e:
                caught = True
                reason = str(e)
            rows.append(dict(name='self-selected prediction benchmark', rejected=caught, positive_restored=True, reason=reason))
            assert caught
        finally:
            atomic(p, orig)
        for names in [[], ['parametric', 'parametric'], ['../escape']]:
            bad = copy.deepcopy(pm)
            bad['participants'] = names
            atomic(p, json.dumps(bad).encode())
            caught = False
            try:
                predictions(proot, sha(p), cases, BENCHMARK_SHA256)
            except IntegrityError:
                caught = True
            finally:
                atomic(p, orig)
            rows.append(dict(name='invalid participant domain ' + str(names), rejected=caught, positive_restored=True))
            assert caught
        for text in [b'{"x":1,"x":2}', b'{"x":NaN}', b'{"x":Infinity}', b'{"x":1e999}']:
            caught = False
            try:
                strict_json(text)
            except IntegrityError:
                caught = True
            rows.append(dict(name='strict JSON ' + text.decode(), rejected=caught, positive_restored=strict_json(b'{"x":1}') == {'x': 1}))
            assert caught
        for name in ['../secret', '/absolute', 'a//b', 'a/./b', 'a\\b']:
            caught = False
            try:
                relative(name)
            except IntegrityError:
                caught = True
            rows.append(dict(name='manifest path ' + name, rejected=caught, positive_restored=relative('a/b') == ('a', 'b')))
            assert caught
        buf = io.BytesIO()
        np.save(buf, np.array([1.0]), allow_pickle=False)
        valid_npy = buf.getvalue()
        header = io.BytesIO()
        np.lib.format.write_array_header_1_0(header, dict(descr='<f8', fortran_order=False, shape=(10 ** 12,)))
        attacks = [('duplicate NPZ key', [('x.npy', valid_npy), ('x.npy', valid_npy)]), ('NPZ path traversal', [('../x.npy', valid_npy)]), ('NPY forged allocation', [('x.npy', header.getvalue())]), ('truncated NPY', [('x.npy', valid_npy[:-1])])]
        good = io.BytesIO()
        np.savez_compressed(good, x=np.array([1.0]))
        for (label, members) in attacks:
            buf = io.BytesIO()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                with zipfile.ZipFile(buf, 'w') as z:
                    for (name, blob) in members:
                        z.writestr(name, blob)
            caught = False
            try:
                load_npz(buf.getvalue())
            except IntegrityError:
                caught = True
            rows.append(dict(name=label, rejected=caught, positive_restored=bool(load_npz(good.getvalue())['x'][0] == 1.0)))
            assert caught
        from scorer import evaluate
        for dest in [scratch / 'payload/public/new_results', proot / 'new_results', scratch / 'code']:
            caught = False
            try:
                evaluate(baseline, proot, PREDICTIONS_SHA256, dest)
            except IntegrityError:
                caught = True
            rows.append(dict(name='output protected tree ' + str(dest.relative_to(scratch)), rejected=caught, positive_restored=True))
            assert caught
        p = proot / ('parametric/' + case + '.npz')
        a = load_npz(p.read_bytes())
        decode(a, ts)
        loaded = load_local(ts[idx], scratch / 'payload/public')
        for (name, change) in [('duplicate task ID', lambda a: a['task_ids'].__setitem__(1, a['task_ids'][0])), ('extra array', lambda a: a.update(extra=np.array([1]))), ('missing array', lambda a: a.pop('outer_' + str(idx))), ('fake UNKNOWN', lambda a: a['status'].__setitem__(idx, 'UNKNOWN_SITE'))]:
            bad = {k: v.copy() for (k, v) in a.items()}
            change(bad)
            caught = False
            try:
                decode(bad, ts)
            except IntegrityError:
                caught = True
            rows.append(dict(name=name, rejected=caught, positive_restored=True))
            assert caught
        for badval in [np.zeros((len(loaded['xy']), 2)), np.full(len(loaded['xy']), np.nan), np.full(len(loaded['xy']), np.inf)]:
            bad = {k: v.copy() for (k, v) in a.items()}
            bad['outer_' + str(idx)] = badval
            caught = False
            try:
                design(bad, idx, loaded, 'DESIGN')
            except IntegrityError:
                caught = True
            rows.append(dict(name='malformed coordinates ' + str(badval.shape), rejected=caught, positive_restored=True))
            assert caught
        from scorer import metrics
        positive_metric = metrics(loaded, a['outer_' + str(idx)], baseline.npz(ref)[ts[idx]['family']], 0.1)[0]
        original_prediction = p.read_bytes()
        manifest_path = proot / 'FROZEN_PREDICTIONS.json'
        original_manifest = manifest_path.read_bytes()
        bad = {k: v.copy() for (k, v) in a.items()}
        bad['outer_' + str(idx)] = np.full(len(loaded['xy']), 1e+200)
        buf = io.BytesIO()
        np.savez_compressed(buf, **bad)
        atomic(p, buf.getvalue())
        forged = copy.deepcopy(pm)
        forged['files']['parametric/' + case + '.npz'] = dict(sha256=sha(p), bytes=p.stat().st_size)
        atomic(manifest_path, json.dumps(forged).encode())
        caught = False
        reason = ''
        try:
            score_submission(scratch, proot, sha(manifest_path), scratch / 'numeric_probe_output')
        except IntegrityError as e:
            caught = str(e) == 'unrepresentable metric arithmetic'
            reason = str(e)
        finally:
            atomic(p, original_prediction)
            atomic(manifest_path, original_manifest)
        rows.append(dict(name='finite coordinate overflow through supported submission', rejected=caught, positive_restored=bool(np.isfinite(positive_metric)), reason=reason))
        assert caught and np.isfinite(positive_metric)
        Bundle(scratch)
        predictions(proot, PREDICTIONS_SHA256, cases, BENCHMARK_SHA256)
        return dict(checks=rows, count=len(rows), all_rejected=all((r['rejected'] for r in rows)), all_positive_restored=True, c02_native_witness=dict(before=old, unguarded_after=new, task_split='dev', guarded='REJECTED'), seconds=time.perf_counter() - t0, threat_model='Trusted evaluator and release anchor; untrusted post-freeze data/predictions. Kernel/root compromise or replacing the trusted evaluator is excluded.')
    finally:
        shutil.rmtree(scratch)

def sandbox_control(bundle):
    from sandbox import command
    out = bundle.payload.parent / 'sandbox_probe'
    out.mkdir(exist_ok=True)
    try:
        script = _release_expand("from pathlib import Path\nimport json\nassert list(Path('/inputs/tasks').glob('*.json'))\ndenied=[]\nfor p in ['/private/COHORT.json','/inputs/../private/COHORT.json','/code/scorer.py','/proc/1/root@DENTAL_WORK_ROOT@/X11/TARGET_CONTEXT_SELECTION.json',SOURCE]:\n try: Path(p).read_bytes()\n except (OSError,PermissionError): denied.append(p)\n else: raise RuntimeError('private read succeeded '+p)\ntry: Path('/inputs/injected.txt').write_text('bad')\nexcept OSError: pass\nelse: raise RuntimeError('input mount writable')\nroutes=Path('/proc/net/route').read_text().splitlines()\nassert len(routes)==1\nprint(json.dumps({'public_read':True,'actual_private_reads_denied':len(denied),'input_write_denied':True,'network_route_rows':len(routes)-1}))\n").replace('SOURCE', repr(str(bundle.payload / 'private/COHORT.json')))
        p = subprocess.run(command(out, ['-c', script]), capture_output=True, text=True, timeout=30)
        if p.returncode:
            raise AssertionError('sandbox failed ' + p.stderr)
        return json.loads(p.stdout)
    finally:
        shutil.rmtree(out)
