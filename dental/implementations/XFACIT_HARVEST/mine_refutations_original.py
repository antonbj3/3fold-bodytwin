"""Pull every swarm report that says it refutes one of our values, and rank them.

Measured 2026-10-04 morning: of 340 completed swarm reports in 12 hours, 338 carry a `refutes_us`
field and 129 say TRUE. I had opened six. Counting reports is not reading them, and a self-declared
refutation that nobody opens is worse than no report, because it looks like coverage.

Ranking is by what can be checked fastest: a report that names OUR value and THEIR value with a
locator is checkable in minutes; one that only asserts a refutation is not.
"""
from __future__ import annotations
from dental_release.paths import expand as _release_expand
import json
import pathlib
import re
STORE = pathlib.Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/storage/research/bunny48_20260926/bodytwin'))
OUT = pathlib.Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/coupled-model/results/ASSEMBLY_REFUTATIONS'))
NUM = re.compile('-?\\d+\\.?\\d*(?:[eE][+-]?\\d+)?')

def main() -> None:
    rows = []
    for f in STORE.glob('*/results.json'):
        try:
            d = json.loads(f.read_text())
        except Exception:
            continue
        ext = d.get('external_referent') or {}
        if not isinstance(ext, dict) or ext.get('refutes_us') is not True:
            continue
        ours = str(ext.get('our_value', ''))
        theirs = str(ext.get('compared_quantity', '')) + ' ' + str(ext.get('locator', ''))
        rows.append({'job': f.parent.name, 'target': d.get('target_id'), 'claim_type': d.get('claim_type'), 'decision': str(d.get('decision', ''))[:300], 'our_value': ours[:400], 'compared_quantity': str(ext.get('compared_quantity', ''))[:300], 'locator': str(ext.get('locator', ''))[:300], 'has_doi_or_pmid': bool(re.search('10\\.\\d{4}|PMID|pmid', theirs)), 'our_numbers': NUM.findall(ours)[:6], 'kind': str(ext.get('kind', ''))})
    rows.sort(key=lambda r: (bool(r['our_numbers']) and r['has_doi_or_pmid'], bool(r['our_numbers'])), reverse=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'REFUTATIONS.json').write_text(json.dumps({'total': len(rows), 'rows': rows, 'review_state': 'PENDING_INDEPENDENT_REVIEW'}, indent=2))
    checkable = [r for r in rows if r['our_numbers'] and r['has_doi_or_pmid']]
    print(f'refutes_us TRUE: {len(rows)}')
    print(f'  with numerical values AND external DOI/PMID (directly verifiable): {len(checkable)}')
    print(f"  with own numbers but without external locator: {len([r for r in rows if r['our_numbers'] and (not r['has_doi_or_pmid'])])}")
    print(f"  without own numbers: {len([r for r in rows if not r['our_numbers']])}")
    print()
    for r in checkable[:12]:
        print(f"{r['job'][-14:]}  {str(r['target'])[:22]:<22} {r['our_numbers']}")
        print(f"""   against: {r['compared_quantity'][:120]}""")
    print('\nwrote', OUT / 'REFUTATIONS.json')
if __name__ == '__main__':
    main()
