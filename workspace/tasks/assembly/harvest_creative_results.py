"""Read all completed creative job reports in the WORKSPACE and pull out consumable numbers.

WHY THIS EXISTS, stated because the miss was mine. `tasks/build_night/scan_swarm.py` reads the
controller's store under external_mount, which holds only the self-generated BT-FW48 follow-up
chain. The creative jobs built from our own net -- BT-NET, BT-CONN, BT-OBST, BT-CAP, BT-IDEA,
BT-ANOM, BT-2ND -- run in this workspace's results/ directory instead. Measured 2026-10-04 morning:
383 such directories exist, 375 carry a RESULTS.md, and all 375 carry a DOI or a PMID. None had been
read, because every scan pointed at the other directory. Counting one directory and calling it the
swarm is the error; this file reads the other one.

WHAT IT PULLS. Per report: the external locators, every number with a unit, and whether the report
names one of our net variables or edge ids. A number is only useful if it can be traced, so a line
without a locator in the same report is kept separately rather than mixed in.

NO GATE IS APPLIED HERE beyond admission. Admission needs the unit-and-range gate and a
decision that consumes the number; this is the reading step that was missing.
"""
from __future__ import annotations

import json
import pathlib
import re

W = pathlib.Path('')
FAMILIES = ('BT-NET-', 'BT-CONN-', 'BT-OBST-', 'BT-CAP-', 'BT-IDEA-', 'BT-ANOM-', 'BT-2ND-')
DOI = re.compile(r'10\.\d{4,9}/[^\s)\]",;]+')
PMID = re.compile(r'PMID[:\s]*(\d{6,9})', re.I)
# a number followed by a unit, which is what makes it consumable rather than decorative
UNIT = r'(?:D|mm|cm|µm|um|nm|pm|mmHg|MPa|kPa|Pa|N|mN|s|ms|µs|us|ns|min|h|%|µM|uM|mM|nM|mol|L|mL|K|°C|Hz|kHz|MHz|GHz|J|W|g|kg|mg|µg)'
NUMUNIT = re.compile(r'(?<![\w.])(-?\d+(?:[.,]\d+)?(?:[eE][+-]?\d+)?)\s*(' + UNIT + r')(?![A-Za-z])')


def net_names() -> set[str]:
    d = json.loads((W / 'CONSTRAINT_NETS.json').read_text())['bodytwin']['tissue_constraint_net']
    out = {v['id'] for v in d['variables']} | {e['id'] for e in d['edges']}
    return out


def main() -> None:
    names = net_names()
    reports, rows = [], []
    for d in sorted(W.glob('results/BT-*')):
        if not d.name.startswith(FAMILIES):
            continue
        f = d / 'RESULTS.md'
        if not f.is_file() or f.stat().st_size == 0:
            continue
        text = f.read_text(errors='replace')
        dois = sorted(set(DOI.findall(text)))[:6]
        pmids = sorted(set(PMID.findall(text)))[:6]
        hits = sorted({n for n in names if n in text})
        nums = NUMUNIT.findall(text)
        reports.append({
            'job': d.name, 'family': d.name.split('-')[1], 'bytes': len(text),
            'dois': dois, 'pmids': pmids, 'net_names_mentioned': hits[:8],
            'numbers_with_units': len(nums),
            'traceable': bool(dois or pmids),
        })
        if dois or pmids:
            seen = set()
            for val, unit in nums:
                key = (val, unit)
                if key in seen:
                    continue
                seen.add(key)
                rows.append({'job': d.name, 'value': val, 'unit': unit,
                             'doi': dois[0] if dois else None,
                             'pmid': pmids[0] if pmids else None,
                             'net_names': hits[:4]})
    out = W / 'results/ASSEMBLY_CREATIVE_HARVEST_V2'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'REPORTS.json').write_text(json.dumps(reports, indent=1, ensure_ascii=False))
    (out / 'NUMBERS.json').write_text(json.dumps(rows, indent=1, ensure_ascii=False))
    fam = {}
    for r in reports:
        fam[r['family']] = fam.get(r['family'], 0) + 1
    print(f'reports read: {len(reports)}  traceable: {sum(r["traceable"] for r in reports)}')
    print('per family:', fam)
    print(f'numbers with units, total: {sum(r["numbers_with_units"] for r in reports)}'
          f'  unique in traceable reports: {len(rows)}')
    linked = [r for r in reports if r['net_names_mentioned']]
    print(f'reports naming a net variable or edge: {len(linked)}')
    for r in sorted(reports, key=lambda x: -x['numbers_with_units'])[:8]:
        print(f'  {r["job"][:46]:<48} {r["numbers_with_units"]:>4} numbers  '
              f'{"DOI" if r["dois"] else ""}{"/PMID" if r["pmids"] else ""}  {r["net_names_mentioned"][:2]}')
    print('wrote', out)


if __name__ == '__main__':
    main()
