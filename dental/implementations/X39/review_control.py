"""Reject any drift in the independently reviewed manuscript and numeric results."""
from pathlib import Path
import hashlib, json, re
P = Path(__file__).resolve().parent

def sha(blob):
    return hashlib.sha256(blob).hexdigest()

def projection(q):
    return {k: v for (k, v) in q.items() if k not in ['full_cost', 'review_numeric_evidence']}

def compare_text(blob, expected):
    return sha(blob) == expected

def compare_result(actual, expected):
    return projection(actual) == expected

def main():
    lock = json.loads((P / 'REVIEW_NUMERIC_LOCK.json').read_text())
    r = json.loads((P / 'results.json').read_text())
    failures = []
    faults = 0
    for (name, expected) in lock['reviewed_document_sha256'].items():
        blob = (P / name).read_bytes()
        if not compare_text(blob, expected):
            failures.append(name + ' differs from reviewed text')
        for m in re.finditer(b'\\d+(?:[.,]\\d+)*', blob):
            bad = blob[:m.start()] + b'999999' + blob[m.end():]
            if compare_text(bad, expected):
                raise AssertionError('numeric mutation accepted')
            faults += 1
    if not compare_result(r, lock['result_projection']):
        failures.append('results.json differs from reviewed numerical projection')
    bad = json.loads(json.dumps(r))
    bad['opposition'][0]['hits'] = 159
    assert not compare_result(bad, lock['result_projection'])
    if r['publication_ready'] is not False:
        failures.append('publication readiness promoted')
    out = {'status': 'FAIL' if failures else 'PASS', 'failures': failures, 'numeric_token_injections_rejected': faults, 'actual_result_injection_rejected': True, 'scope': 'Whole reviewed output identity plus independently recounted fixed results; no clinical validation', 'review_numeric_evidence': 'REVIEW_NUMERIC_EVIDENCE.json', 'retrospective_review_lock': True}
    (P / 'REVIEW_CONTROL_RESULT.json').write_text(json.dumps(out, indent=2) + '\n')
    if failures:
        raise AssertionError('; '.join(failures))
    print('Review control PASS:', faults, 'numeric text mutations and one result mutation rejected.')
if __name__ == '__main__':
    main()
