"""Recover original receipt bytes only when their captured SHA256 agrees exactly."""
import hashlib, json, itertools
from pathlib import Path
from prepare_replay import HERE, write

def recover(obj, expected):
    for (indent, ascii_only, sorted_keys, newline) in itertools.product([2, None, 4], [False, True], [False, True], ['\n', '']):
        data = (json.dumps(obj, indent=indent, ensure_ascii=ascii_only, sort_keys=sorted_keys) + newline).encode()
        if hashlib.sha256(data).hexdigest() == expected:
            return data
    return None

def main():
    sources = list((HERE / 'raw/prior_runs').glob('*/receipts/*.json')) + list((HERE / 'raw/receipts').glob('*.json'))
    records = []
    missed = []
    for p in sources:
        r = json.loads(p.read_text())
        if not r.get('validator_sha256'):
            continue
        raw = recover(r['validator_receipt'], r['validator_sha256'])
        if raw is None:
            missed.append({'receipt': str(p), 'sha256': r['validator_sha256']})
            continue
        dest = p.parent.parent / 'original_validator_bytes' / p.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            assert dest.read_bytes() == raw
        else:
            dest.write_bytes(raw)
        records.append({'wrapper_receipt': str(p), 'original_bytes_locator': str(dest), 'original_sha256': r['validator_sha256'], 'bytes': len(raw), 'scope': 'Byte-exact reconstruction from captured original JSON; accepted only by the independently captured original byte hash'})
    result = {'captured_original_byte_hashes': len(records) + len(missed), 'exact_original_bytes_recovered': len(records), 'unrecovered': missed, 'records': records, 'pass': not missed, 'fallback': 'Never synthesize a receipt with a different hash; retain semantic receipt and report the missing exact-byte provenance'}
    write(HERE / 'raw/ORIGINAL_VALIDATOR_BYTE_AUDIT.json', result)
    print(json.dumps({k: result[k] for k in ['captured_original_byte_hashes', 'exact_original_bytes_recovered', 'unrecovered', 'pass']}))
    if missed:
        raise SystemExit(1)
if __name__ == '__main__':
    main()
