from pathlib import Path
import json, hashlib
H = Path(__file__).resolve().parent
R = H.parent.parent
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
required = ['CURRENT_WORK_STATE.json', 'PREREG_INFOLINK4_R1.json', 'PREREG_INFOLINK4_R2.json', 'PREREG_INTERACTION_R3.json', 'FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R3.json', 'RESULTS.md', 'results.json', 'COMMANDS.md', 'HANDOFF.md', 'README_DEMO.md', 'GRAPH_FEEDBACK.json', 'FACIT.csv', 'information_links_figure.png', 'run_all.sh', 'SOURCE_MANIFEST.json', 'SEMANTIC_REVIEW_R1.json', 'DECOMPOSITION_PER_LINK.json', 'GRAPH_INSTALL_CHECK.json']
assert all(((H / x).is_file() for x in required))
j = json.loads((H / 'results.json').read_text())
g = json.loads((H / 'GRAPH_INSTALL_CHECK.json').read_text())
f = json.loads((H / 'GRAPH_FEEDBACK.json').read_text())
v = json.loads((H / 'VALIDATION.json').read_text())
s = json.loads((H / 'CURRENT_WORK_STATE.json').read_text())
assert g['created_file_sha256'] == sha(Path(g['created_file'])) == sha(H / 'information_links.preview.jsonl')
assert f['sha256'] == sha(H / 'results.json')
assert s['phase'] == 'COMPLETE_PENDING_INDEPENDENT_REVIEW' and s['graph_installed']
assert j['established_links'] == g['link_count'] == 11
assert j['attrition']['assessed_candidates'] == j['established_links'] + j['attrition']['rejected'] == 21
assert j['fault_injections_rejected_total'] == 132 and v['fault_injections'] == 130
assert all((x['rejected'] for x in json.loads((H / 'FAULT_INJECTION.json').read_text())))
for p in H.glob('*.sha256'):
    assert sha(H / p.name[:-7]) == p.read_text().strip(), p.name
before = json.loads((H / 'EXPANSION_HASHES_BEFORE.json').read_text())
assert all((sha(R / 'notes/expansion' / n) == h for (n, h) in before.items()))
assert max((p.stat().st_size for p in H.rglob('*') if p.is_file())) < 50000000
assert sum((p.stat().st_size for p in H.rglob('*') if p.is_file())) < 3000000000
check = {'status': 'PASS', 'new_link_count': 11, 'required_artifacts': len(required), 'source_value_and_semantic_review': 'PENDING_INDEPENDENT_REVIEW', 'feedback_matches_result_hash': True, 'expansion_matches_checked_preview': True, 'other_expansions_unchanged': True, 'frozen_files_unchanged': True, 'all_numeric_controls_and_counterexamples_passed': True, 'package_bytes_at_check': sum((p.stat().st_size for p in H.rglob('*') if p.is_file())), 'large_artifacts_over50MB': 0}
(H / 'DELIVERY_CHECK.json').write_text(json.dumps(check, indent=2) + '\n')
print(json.dumps(check, indent=2))
