"""Check delivery bindings without rerunning numerical experiments."""
from dental_release.paths import expand as _release_expand
import ast
import datetime
import hashlib
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_CAD_10X'))

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda : f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def read(path):
    return json.loads((ROOT / path).read_text())

def main():
    checks = {}
    for p in (ROOT / 'code').glob('*.py'):
        ast.parse(p.read_text())
    checks['python_source_syntax'] = True
    excluded = []
    for p in ROOT.rglob('*.json'):
        if p == ROOT / 'raw/GRAPH_PACKET.json' and p.stat().st_size == 0:
            excluded.append({'path': str(p.relative_to(ROOT)), 'reason': 'Preserved empty stdout of rejected graph packet command; see COMMANDS.md (stderr was not saved as a separate file)'})
            continue
        json.loads(p.read_text())
    checks['json_parse_except_preserved_failed_command_stdout'] = True
    frozen = []
    for p in ROOT.glob('*.sha256'):
        candidate = p.with_suffix('.json')
        if candidate.exists():
            assert sha(candidate) == p.read_text().strip(), str(candidate)
            frozen.append(str(candidate.relative_to(ROOT)))
    checks['frozen_hashes'] = len(frozen)
    deps = read('DEPENDENCIES.json')['files']
    for r in deps:
        assert sha(r['path']) == r['sha256'], r['path']
    checks['read_only_dependency_hashes'] = len(deps)
    rec = read('raw/DEMO_RECEIPT_FULL.json')
    assert rec['results_sha256'] == sha(ROOT / 'results.json')
    checks['full_replay_result_binding'] = True
    feedback = read('GRAPH_FEEDBACK.json')
    assert feedback['sha256'] == sha(ROOT / 'REVIEW_SNAPSHOT.json')
    assert Path(read('raw/GRAPH_FEEDBACK_RECEIPT.json')['receipt']).exists()
    checks['graph_snapshot_and_receipt'] = True
    controls = read('raw/CONTROLS.json')
    assert controls['all_pass'] and all((r['pass_'] for r in controls['checks']))
    checks['controls_pass'] = controls['count']
    for (name, count) in [('REGENERATION_ALL', 58), ('MARGIN_REPLAY_ALL', 130)]:
        r = read('raw/' + name + '.json')
        assert r['all_identical'] and r['requested'] == count
    checks['identical_replays'] = {'crowns': 58, 'margin_inputs': 130}
    result = read('results.json')
    assert result['external_referent']['refutes_us']
    assert result['summary']['full_no_manual_physical_rate'] == 'UNKNOWN_NOT_MEASURED'
    checks['negative_external_result_and_physical_unknown'] = True
    out = dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), checks=checks, excluded_failed_command_outputs=excluded)
    (ROOT / 'raw/DELIVERY_AUDIT.json').write_text(json.dumps(out, indent=2) + '\n')

    def items(paths, root=None):
        return [dict(path=str(p.relative_to(root) if root else p), bytes=p.stat().st_size, sha256=sha(p)) for p in paths]
    files = items(sorted((p for p in ROOT.rglob('*') if p.is_file() and p.name != 'ARTIFACT_MANIFEST.json' and ('.tmp' not in p.name))), ROOT)
    datafiles = items(sorted((p for p in DATA.rglob('*') if p.is_file())))
    n = sum((p['bytes'] for p in datafiles))
    assert n < 3000000000
    (ROOT / 'ARTIFACT_MANIFEST.json').write_text(json.dumps(dict(files=files, data_files=datafiles, data_bytes=n, root=str(DATA)), indent=2) + '\n')
    print(json.dumps(dict(checks=checks, lane_files=len(files), data_files=len(datafiles), data_bytes=n)))
if __name__ == '__main__':
    main()
