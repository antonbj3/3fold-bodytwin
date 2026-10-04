from pathlib import Path
import json, hashlib, datetime, shutil
H = Path(__file__).resolve().parent
for n in ['results.json', 'information_links.preview.jsonl', 'VALIDATION.json', 'FAULT_INJECTION.json', 'SUFFICIENCY_TESTS.json']:
    shutil.copyfile(H / n, H / (n + '.R1'))
review = {'reviewer': 'producer self-review, NOT independent review', 'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'numeric_extraction': 'PASS', 'semantic_gate': 'FAIL 2 time-scale assignments', 'failures': [{'link': k, 'field': 'time_scale', 'before': 'HANDOVER', 'after': 'SIMULTANEOUS', 'reason': 'Source measures wear and concurrent support/counterface under the same endpoint protocol. No separate fast-to-slow process coupling is established by this observation link.'} for k in ['L06', 'L07']], 'thresholds_changed': False, 'original_contract_preserved': 'LINK_CONTRACTS.json', 'revised_contract': 'LINK_CONTRACTS_V2.json'}
(H / 'SEMANTIC_REVIEW_R1.json').write_text(json.dumps(review, indent=2) + '\n')
c = json.loads((H / 'LINK_CONTRACTS.json').read_text())
for x in c:
    if x['id'] in ['L06', 'L07']:
        x['time_scale'] = 'SIMULTANEOUS'
p = H / 'LINK_CONTRACTS_V2.json'
p.write_text(json.dumps(c, indent=2, ensure_ascii=False) + '\n')
(H / (p.name + '.sha256')).write_text(hashlib.sha256(p.read_bytes()).hexdigest() + '\n')
print('Original semantic failure saved; V2 contracts frozen before rerun.')
