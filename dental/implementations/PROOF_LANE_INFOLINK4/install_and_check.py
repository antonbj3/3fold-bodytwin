"""Only shared authored write: the new expansion file explicitly requested by user."""
from pathlib import Path
import json, hashlib, subprocess, datetime
H = Path(__file__).resolve().parent
R = H.parent.parent
p = R / 'notes/expansion/information_links_astra4_20261003.jsonl'
b = (H / 'information_links.preview.jsonl').read_bytes()

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
val = json.loads((H / 'VALIDATION.json').read_text())
assert val['all_passed'] and val['all_faults_rejected']
rows = [json.loads(x) for x in b.decode().splitlines()]
assert 10 <= len(rows) <= 15
ids = [x['id'] for x in rows]
assert len(ids) == len(set(ids))
existing = []
for f in p.parent.glob('*.jsonl'):
    if f == p:
        continue
    for line in f.read_text().splitlines():
        try:
            existing.append(json.loads(line))
        except json.JSONDecodeError:
            pass
oldids = {x['id'] for x in existing if 'id' in x}
assert not set(ids) & oldids
chains = {r['chain_id'] for r in []}
import csv
with (R / 'results/PROOF_LANE_INFOLINK/CHAIN_COVERAGE.csv').open() as f:
    for r in csv.DictReader(f):
        chains.add(r.get('chain_id', r.get('id')))
for x in rows:
    assert x['claim_type'] == 'information_link' and x['review_state'] == 'PENDING_INDEPENDENT_REVIEW' and (x['status'] == 'OPEN') and (x['depends_on'] == [])
    rel = x['relations'][0]
    if rel['target_kind'] == 'draft_node':
        assert x['consumer'] in oldids
    else:
        assert x['consumer'] in chains, (x['consumer'], chains)
    if 'source_node_id' in x:
        assert x['source_node_id'] in oldids
    assert x['resolution_level'] in ('POPULATION', 'PER_POINT')
    assert x['time_scale'] in ('SIMULTANEOUS', 'HANDOVER')
    assert x['source_contract']
before = {f.name: sha(f) for f in p.parent.glob('*.jsonl') if f != p}
view = R / 'results/GRAPH_WORKING_VIEW_20260923/WORKING_VIEW.json'
oldview = json.loads(view.read_text())

def dictionaries(obj):
    if isinstance(obj, dict):
        yield obj
        for v in obj.values():
            yield from dictionaries(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from dictionaries(v)
old_seen = {v.get('id') for v in dictionaries(oldview) if isinstance(v.get('id'), str)}
if p.exists():
    assert p.read_bytes() == b, 'Refuse to overwrite different authored content'
else:
    with p.open('xb') as f:
        f.write(b)
r = subprocess.run(['./graph', 'working', 'build'], cwd=R, capture_output=True, text=True)
(H / 'GRAPH_BUILD.stdout').write_text(r.stdout)
(H / 'GRAPH_BUILD.stderr').write_text(r.stderr)
assert r.returncode == 0, (r.returncode, r.stderr)
v = json.loads(view.read_text())
seen = {a.get('id') for a in dictionaries(v) if isinstance(a.get('id'), str)}
assert set(ids) <= seen
allrows = list(dictionaries(v))
matched = []
for row in rows:
    hits = [a for a in allrows if a.get('id') == row['id']]
    assert hits
    assert any((a.get('claim') == row['claim'] for a in hits))
    matched.append(row['id'])
after = {f.name: sha(f) for f in p.parent.glob('*.jsonl') if f != p}
assert before == after, 'Concurrent change in other expansions requires attribution check'
start = json.loads((H / 'EXPANSION_HASHES_BEFORE.json').read_text())
check = {'created_file': str(p), 'created_file_sha256': sha(p), 'link_count': len(rows), 'all_imported': True, 'imported_ids': matched, 'ids_new_relative_to_preinstall_view': sorted(set(ids) - old_seen), 'other_expansion_files_unchanged_during_install': before == after, 'other_expansion_files_changed_since_round_start': [k for k in start if start[k] != after.get(k)], 'other_expansion_files_added_since_round_start': sorted(set(after) - set(start)), 'build_receipt': json.loads(r.stdout), 'checked_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'write_policy': 'Authored one new expansion file; working build is expressly authorized, no native refresh or source-graph mutation.'}
(H / 'GRAPH_INSTALL_CHECK.json').write_text(json.dumps(check, indent=2, ensure_ascii=False) + '\n')
(H / 'EXPANSION_HASHES_AFTER.json').write_text(json.dumps(after, indent=2) + '\n')
print(json.dumps(check, indent=2))
