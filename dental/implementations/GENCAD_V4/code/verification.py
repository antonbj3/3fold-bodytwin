"""Run real rejection probes and process isolation; restore every injected byte."""
from dental_release.paths import expand as _release_expand
from common import *
from integrity import verify, IntegrityError, lock, checked_json_bytes
from generate import command
import subprocess, socket

def run():
    start = time.perf_counter()
    m = lock()
    verify()
    checks = []
    prefixes = ['payload/public/tasks/', 'payload/public/scenes/', 'payload/private/references/', 'payload/predictions/population/', 'payload/whole_inputs/', 'payload/whole_private/', 'payload/whole_predictions/', 'payload/models/', 'payload/literature/', 'payload/participant_code/field/field_generator.py', 'code/quality.py', 'PREREG_R2.json', 'payload/runtime/']
    for prefix in prefixes:
        rel = next((k for k in m['files'] if k.startswith(prefix) and (not k.endswith('bin/python3')) and (not k.endswith('.so'))))
        p = ROOT / rel
        blob = p.read_bytes()
        try:
            p.write_bytes(blob + b'\nINJECTED_WRONG_VALUE\n')
            rejected = False
            try:
                verify()
            except IntegrityError:
                rejected = True
            cp = subprocess.run([str(DATA / 'runtime/bin/python3'), '-s', '-E', '-B', str(ROOT / 'code/score_panel.py')], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            checks.append(dict(name='tamper_' + prefix, relative=rel, manifest_rejected=rejected, score_entry_rejected=cp.returncode != 0, diagnostic=cp.stderr[-250:]))
        finally:
            p.write_bytes(blob)
    extra = DATA / 'public/UNDECLARED_FILE'
    try:
        extra.write_text('wrong extra input')
        try:
            verify()
            rejected = False
        except IntegrityError:
            rejected = True
        checks.append(dict(name='undeclared_input', rejected=rejected))
    finally:
        extra.unlink()
    rel = next((k for k in m['files'] if k.startswith('payload/public/tasks/')))
    p = ROOT / rel
    blob = p.read_bytes()
    verify()
    try:
        p.write_bytes(blob.replace(b'"mm"', b'"cm"', 1))
        try:
            checked_json_bytes(p)
            rejected = False
        except IntegrityError:
            rejected = True
        checks.append(dict(name='read_after_verify_units_mutation', rejected=rejected))
    finally:
        p.write_bytes(blob)
    output = DATA / 'probe_output'
    output.mkdir(exist_ok=True)
    script = _release_expand("import os,json,socket\nfrom pathlib import Path\nout={}\nfor p in ['/inputs/../private','/whole_private','@DENTAL_EXTERNAL_ROOT@/projects','@DENTAL_EXTERNAL_ROOT@/storage/datasets']:\n out[p]=not Path(p).exists()\ntry:\n Path('/inputs/WRITE_PROBE').write_text('wrong');out['input_read_only']=False\nexcept OSError:out['input_read_only']=True\ns=socket.socket();s.settimeout(.2)\ntry:s.connect(('192.0.2.1',9));out['network_blocked']=False\nexcept OSError:out['network_blocked']=True\ns.close()\nPath('/output/ISOLATION.json').write_text(json.dumps(out))\n")
    cmd = command(output)
    cmd = cmd[:-1] + ['-c', script]
    cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    iso = read(output / 'ISOLATION.json') if cp.returncode == 0 else {'process_failed': False}
    checks.append(dict(name='public_generator_isolation', passed=cp.returncode == 0 and all(iso.values()), observations=iso))
    split = read(DATA / 'private/MODEL_SPLIT_PROOF.json')
    checks.append(dict(name='train_dev_test_group_separation', passed=not split['x42_test_overlap'] and (not split['template_test_overlap'])))
    from quality import roof_metrics
    import scipy.sparse as sp
    t = dict(xy=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]), faces=np.array([[0, 1, 2]]), weights=np.ones(3), preparation_z=np.zeros(3), ceiling=np.ones(3))
    for (name, value) in [('nonfinite', np.nan), ('overflow', 1e+308)]:
        try:
            roof_metrics(t, np.full(3, value), np.zeros(3), np.zeros(3))
            rejected = False
        except (ValueError, FloatingPointError):
            rejected = True
        checks.append(dict(name=name + '_metric_rejected', rejected=rejected))
    restored = verify()
    ok = all((all((v for (k, v) in c.items() if k in ['rejected', 'passed', 'manifest_rejected', 'score_entry_rejected'])) for c in checks))
    out = dict(passed=ok, controls=checks, restored_integrity=restored, seconds=time.perf_counter() - start, scope='Execution and rejection tests, not independent scientific review')
    dump(ROOT / 'raw/VERIFICATION.json', out)
    if not ok:
        raise RuntimeError('Control failed')
    print('verification', len(checks), 'PASS')
if __name__ == '__main__':
    run()
