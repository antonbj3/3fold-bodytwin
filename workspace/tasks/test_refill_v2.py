import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('refill_v2', ROOT / 'tasks/refill_v2.py')
refill = importlib.util.module_from_spec(spec)
spec.loader.exec_module(refill)


def test_graph_bindings_and_typed_pairs():
    mapping = json.loads((ROOT / 'results/CX-QGRAPH/Q_GRAPH_MAP.json').read_text())
    pairs = json.loads((ROOT / 'results/CX-QGRAPH/PAIRS.json').read_text())['pairs']
    assert len(mapping['questions']) == 168
    assert all(bool(q['matches']) != bool(q.get('new_node_proposed')) for q in mapping['questions'])
    assert pairs
    native = json.loads((ROOT / 'MERGED_GRAPH.json').read_text())
    dep = {(e['from'], e['to']) for e in native['links'] if e['kind'] == 'depends_on'}
    couples = {(a['from'], t) for a in json.loads((ROOT / 'COUPLING_ANNOTATIONS.json').read_text()) for t in a.get('resolved_targets', [])}
    for pair in pairs:
        assert len(pair['questions']) == 2
        for e in pair['evidence']:
            path = e['path']
            if e['kind'].startswith('depends_on'):
                assert all((a, b) in dep or (b, a) in dep for a, b in zip(path, path[1:]))
            else:
                assert e['kind'] == 'coupling_annotation'
                assert (path[0], path[1]) in couples or (path[1], path[0]) in couples


def test_criterion_requires_explicit_pass(tmp_path):
    d = tmp_path / 'packet'
    d.mkdir()
    (d / 'RESULTS.md').write_text('criterion met')
    (d / 'results.json').write_text('{"criteria": {"all_pass": false}}')
    assert not refill.criterion_met(d)
    (d / 'results.json').write_text('{"criteria": {"all_pass": true}}')
    assert refill.criterion_met(d)


def test_preview_is_read_only_and_audit_has_ten():
    queue = ROOT / 'tasks/lanes/bt_queue.txt'
    before = queue.read_bytes()
    selected = refill.plan(80, now=1_800_000_000)
    assert queue.read_bytes() == before
    audit = [x for x in selected if x[0].startswith('BT-AUD-')]
    assert len(audit) == 1
    names = [name for _, name in audit[0][2] if name.endswith('_RESULTS.md')]
    assert len(names) == 10
    assert all(x[0].startswith(('BT-DAT-', 'BT-DATX-', 'BT-DATG-', 'BT-DATQ-', 'BT-INT-', 'BT-VAL-', 'BT-VAL2-', 'BT-CPL-', 'BT-DENT2-', 'BT-AUD-')) for x in selected)
    assert refill.robust_models() == ['Q014', 'Q018', 'Q029', 'Q036']


def test_new_families_have_capacity_and_keep_public_data_separate():
    candidates = refill.graph_data_candidates()
    assert len(candidates) >= 300
    assert len({refill.graph_data_id(node) for node, _, _ in candidates}) == len(candidates)
    selected = refill.plan(80, now=1_800_000_000)
    assert len(selected) >= 80
    names = [item[0] for item in selected]
    for prefix, expected in (('BT-DATG-', candidates),
                             ('BT-DATQ-', refill.QUALITY_TARGETS),
                             ('BT-VAL2-', refill.integrated_models()),
                             ('BT-DENT2-', refill.DENTAL_FOLLOWUPS)):
        assert expected
        assert any(j.startswith(prefix) for j in names) or any(
            d.name.startswith(prefix) for d in refill.R.iterdir())
    assert sum(j.startswith(('BT-DAT-', 'BT-DATX-', 'BT-DATG-', 'BT-DATQ-')) for j in names[:48]) >= 45
    for j, brief, files in selected:
        if j.startswith('BT-DATG-'):
            assert not files
            assert 'the collaborator/GC/restricted model data-data' in brief
