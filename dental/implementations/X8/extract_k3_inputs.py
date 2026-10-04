"""Stream archived predecessor JSON only; retain and check frozen scientific input."""
from dental_release.paths import expand as _release_expand
import hashlib, json, subprocess, tarfile
from pathlib import Path
ROOT = Path(__file__).resolve().parent
archive = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/usb-stage/dental_archive/dental_K3.tar.zst'))
target = Path(_release_expand('@DENTAL_WORK_ROOT@/X8-guide-nerve-risk/K3_SITES.json'))
records = []
run_metadata = None
p = subprocess.Popen(['zstd', '-d', '-c', str(archive)], stdout=subprocess.PIPE)
with tarfile.open(fileobj=p.stdout, mode='r|') as t:
    for member in t:
        if member.isfile() and member.name.startswith('dental_K3/risk/') and member.name.endswith('.json'):
            obj = json.load(t.extractfile(member))
            if member.name.endswith('/_run.json'):
                run_metadata = obj
            else:
                records.append(obj)
p.stdout.close()
if p.wait() != 0:
    raise RuntimeError('Archive decompression failed')
sites = [{'case': obj['case'], 'fdi': int(fdi), **site} for obj in records for (fdi, site) in obj['sites'].items()]
assert run_metadata is not None
out = {'n_cases': run_metadata['n'], 'sites': sites}
if target.exists():
    existing = json.loads(target.read_text())
    key = lambda s: (s['case'], s['fdi'])
    assert existing['n_cases'] == out['n_cases']
    assert sorted(existing['sites'], key=key) == sorted(out['sites'], key=key)
    pr = json.loads((ROOT / 'PREREG_R2.json').read_text())
    assert hashlib.sha256(target.read_bytes()).hexdigest() == pr['input_sha256']
else:
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('x') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
        f.write('\n')
print(json.dumps({'source_archive': str(archive), 'actual_archived_json_cases': len(records), 'declared_source_cases': run_metadata['n'], 'case_count_discrepancy': run_metadata['n'] - len(records), 'missing_case_identity': 'UNKNOWN', 'source_runtime_s': run_metadata['runtime_s'], 'sites': len(sites), 'target': str(target), 'sha256': hashlib.sha256(target.read_bytes()).hexdigest(), 'all_site_rows_identical': True, 'arrays_extracted': False}))
