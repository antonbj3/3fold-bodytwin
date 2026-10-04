"""Turn finished creative jobs into proposals the net and the cells can actually take.

The gap this closes. 135 creative jobs are complete, 83 of them carry at least one paragraph with a
number, a unit and a DOI or PMID, and nothing consumes any of it. Reading them is not processing them:
a number that stays in a report changes nothing, which is the same defect as the 131 cited references
nobody harvested and the 389 lane rounds that never reached the graph.

What this does and deliberately does not do. It extracts candidate records -- quantity, value, unit,
locator, the job it came from -- and matches them against the variables already in the constraint net
and the quantities our cells compute. It writes PROPOSALS. It admits nothing. Automatic admission is
how the 500-entry brainstorm and the 69 bounded jobs put unusable material into the pipeline, and three
separate corpus minings failed tonight because the same key name means different things in different
jobs. So every proposal carries what it would attach to and why, and a human or a stronger model
decides.

The matching is dimensional first and lexical second, because that order is what worked: unit
agreement is a hard constraint, a shared word is a hint. Compound units are read longest-first, since
reading `mu_pa_s` as seconds mislabelled 7511 rows earlier tonight.
"""
from __future__ import annotations

import glob
import json
import os
import re
from collections import Counter
from pathlib import Path

W = Path('')
NET = W / 'CONSTRAINT_NETS.json'
CELLS = W / 'tasks/free48/sources'
OUT = W / 'results/ASSEMBLY_CREATIVE_HARVEST'

LOC = re.compile(r'(10\.\d{4,9}/[^\s"\',\)\]]+|PMID[:\s]*\d{6,9})', re.I)
# number followed by a unit, the unit written as it appears in prose
NUM = re.compile(r'(-?\d+(?:[.,]\d+)?(?:[eE][-+]?\d+)?)\s*'
                 r'(µm|um|mm|nm|cm|m\b|Pa·s|Pa\*s|MPa|kPa|Pa|N·m|N/mm|N\b|J/mol|J/kg|J\b|K\b|°C|degC|'
                 r'deg|°|%|nM|µM|uM|mol/s|mol\b|mL/min|mL|min\b|h\b|s\b|D\b|Hz)')
FAMILIES = ('CAP', 'CONN', 'OBST', 'ANOM', 'NET', '2ND')


def net_variables() -> dict:
    n = json.loads(NET.read_text())['bodytwin']['tissue_constraint_net']
    return {(v['id'] if isinstance(v, dict) else v): (v.get('unit') if isinstance(v, dict) else None)
            for v in n['variables']}


def cell_words() -> dict[str, set[str]]:
    out = {}
    for c in CELLS.iterdir():
        if not c.is_dir():
            continue
        words = set()
        for f in c.rglob('*'):
            if f.is_file() and f.suffix in ('.py', '.json', '.md') and f.stat().st_size < 400_000:
                try:
                    words |= set(re.findall(r'[a-z][a-z0-9_]{4,}',
                                            f.read_text(errors='replace').lower()))
                except Exception:
                    pass
        out[c.name] = words
    return out


def main() -> int:
    variables, cells = net_variables(), cell_words()
    records, per_job = [], Counter()

    for fam in FAMILIES:
        for p in glob.glob(str(W / f'results/BT-{fam}-*/RESULTS.md')):
            job = os.path.basename(os.path.dirname(p))
            text = Path(p).read_text(errors='replace')
            for para in re.split(r'\n\s*\n', text):
                loc = LOC.search(para)
                if not loc:
                    continue
                nums = NUM.findall(para)
                if not nums:
                    continue
                # the quantity name: the longest lower-case token run near the number
                words = set(re.findall(r'[a-z][a-z0-9_]{4,}', para.lower()))
                unit_for_para = None
                # Lexical overlap alone attached 762 of 764 records, which is the loose-matcher
                # failure for the fourth time tonight. A record now attaches to a net variable only
                # if the UNIT AGREES with that variable's declared unit as well, and to a cell only
                # on a strong word overlap. Unit agreement is the hard constraint; a shared word is
                # a hint.
                def unit_matches(var_unit, u):
                    if not var_unit:
                        return False
                    norm = {'µm': 'um', '°': 'deg', '°C': 'degC', 'Pa·s': 'Pa*s', 'µM': 'uM'}
                    return norm.get(u, u) == var_unit
                hit_cells = sorted(c for c, cw in cells.items() if len(words & cw) >= 12)[:3]
                for value, unit in nums[:6]:
                    hit_vars = sorted(v for v, vu in variables.items()
                                      if any(w in v or v in w for w in words if len(w) > 5)
                                      and unit_matches(vu, unit))
                    records.append({
                        'job': job, 'family': fam,
                        'value': value.replace(',', '.'), 'unit': unit,
                        'locator': loc.group(1),
                        'attaches_to_net_variables': hit_vars[:4],
                        'attaches_to_cells': hit_cells,
                        'paragraph': ' '.join(para.split())[:300],
                    })
                per_job[job] += 1

    attachable = [r for r in records if r['attaches_to_net_variables'] or r['attaches_to_cells']]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'PROPOSALS_V1.json').write_text(json.dumps({
        '_meta': {
            'jobs_read': len(per_job),
            'records_with_number_unit_and_locator': len(records),
            'records_that_attach_somewhere': len(attachable),
            'admitted': 0,
            'why_nothing_is_admitted': ('automatic admission is how unusable material entered the '
                                        'pipeline before; matching is dimensional first and lexical '
                                        'second, and a shared word is a hint rather than a match'),
            'units_seen': dict(Counter(r['unit'] for r in records).most_common(12)),
            'families': dict(Counter(r['family'] for r in records)),
            'review_state': 'PENDING_INDEPENDENT_REVIEW',
        },
        'proposals': attachable,
    }, indent=1, ensure_ascii=False))

    print(f'  job read: {len(per_job)}')
    print(f'  poster med tal + enhet + lokalisator: {len(records)}')
    print(f'  of which attached to a network variable or cell: {len(attachable)}')
    print('  per familj: ' + ', '.join(f'{k} {v}' for k, v in
                                       Counter(r['family'] for r in attachable).most_common()))
    print('  vanligaste enheter: ' + ', '.join(f'{u} {n}' for u, n in
                                               Counter(r['unit'] for r in attachable).most_common(8)))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
