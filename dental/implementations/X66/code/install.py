import subprocess, json
from common import *

def main():
    preview = ROOT / 'information_links.preview.jsonl'
    target = DENTAL / 'notes/expansion/information_links_x66_20261003.jsonl'
    payload = preview.read_bytes()
    native = {}
    for p in [DENTAL / 'GRAPH.json', DENTAL / 'MERGED_GRAPH.json', DENTAL / 'PROFILE.json']:
        if p.is_file():
            native[str(p)] = sha(p)
    existing = {str(p): sha(p) for p in (DENTAL / 'notes/expansion').glob('*.jsonl')}
    if target.exists():
        if target.read_bytes() != payload:
            raise SystemExit('Refuse to overwrite different shared expansion file')
    else:
        with target.open('xb') as f:
            f.write(payload)
    proc = subprocess.run(['./graph', 'working', 'build'], cwd=DENTAL, capture_output=True, text=True)
    (ROOT / 'GRAPH_BUILD.stdout').write_text(proc.stdout)
    (ROOT / 'GRAPH_BUILD.stderr').write_text(proc.stderr)
    if proc.returncode:
        raise SystemExit('working build failed; logs preserved')
    expected = [json.loads(l)['id'] for l in payload.decode().splitlines()]
    found = []
    missing = []
    details = []
    for tid in expected:
        chk = subprocess.run(['./graph', 'working', 'packet', '--id', tid], cwd=DENTAL, capture_output=True, text=True)
        (ROOT / 'source_review' / ('INSTALLED_' + tid + '.stdout')).write_text(chk.stdout)
        (ROOT / 'source_review' / ('INSTALLED_' + tid + '.stderr')).write_text(chk.stderr)
        details.append({'id': tid, 'exitcode': chk.returncode, 'packet_bytes': len(chk.stdout)})
        if chk.returncode == 0 and tid in chk.stdout:
            found.append(tid)
        else:
            missing.append(tid)
    after = {str(p): sha(p) for p in (DENTAL / 'notes/expansion').glob('*.jsonl')}
    changed = [p for (p, h) in existing.items() if after.get(p) != h]
    native_change = [p for (p, h) in native.items() if sha(p) != h]
    put('EXPANSION_HASHES_AFTER.json', after)
    result = {'installed_path': str(target), 'sha256': sha(target), 'bytes': target.stat().st_size, 'expected_ids': expected, 'imported_ids': found, 'missing_ids': missing, 'build_exitcode': proc.returncode, 'build_stdout': proc.stdout, 'other_expansion_files_changed_during_install': changed, 'native_graph_files_changed': native_change, 'authored_shared_files': [str(target)], 'all_imported': not missing, 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'packet_checks': details}
    put('GRAPH_INSTALL_CHECK.json', result)
    assert not missing and (not changed) and (not native_change)
    print(json.dumps(result, ensure_ascii=False, indent=2))
if __name__ == '__main__':
    main()
