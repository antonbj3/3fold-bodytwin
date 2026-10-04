"""Briefs for the axis the mandate actually names: what becomes askable when two cells are CONNECTED.

Why this is the supply line and not more of the same. The creative queue holds 178 jobs and drains in
about two hours at the measured rate, and the three generators behind it are each exhausted by
construction: one brief per undetermined edge, one per lane obstacle, one per cell. Those are all
one-per-thing. Connections are one-per-PAIR, so the supply is combinatorial rather than finite, and it
is the axis the project is judged on — things get connected and the twin decides.

How a pair is chosen, and why not all pairs. Two cells are worth asking about when they already share a
dimension, because a shared unit is the cheapest evidence that a quantity could pass between them. That
is a weak filter and it is meant to be: a strong filter would encode my guess about which couplings
matter, and my guesses are what the operator has been correcting all evening. 45 cells give 990 pairs;
the shared-unit filter cuts that to the ones where a port is even conceivable, and the pairs are then
ranked by how many dimensions they share, on the assumption that more shared dimensions means more
surface to connect across. That ranking is an assumption and is stated as one.

What the brief carries: both cells' own descriptions, the dimensions they share, and the quantities each
holds in those dimensions. Nothing about why the pair might matter, because that is the question.
"""
from __future__ import annotations

import itertools
import json
import re
from collections import defaultdict
from pathlib import Path

W = Path('.')
CELLS = W / 'tasks/free48/sources'
QUEUE = W / 'tasks/lanes/bt_queue.txt'
BLOCK = re.compile(r'excluded_category|excluded_category|excluded_category|sinusoid|excluded_category|excluded_category|excluded_category', re.I)
QUANT = re.compile(r'["\']([a-z0-9_]{4,}?_(um|mm|nm|pa|mpa|kpa|j|n|k|deg|degc|s|h|mol|percent|'
                   r'ml_min|pa_s|mol_s|kg_m3|j_mol|m_s|per_s))["\']', re.I)
MAX_PAIRS = 120


def main() -> int:
    info = {}
    for cell in sorted(p for p in CELLS.iterdir() if p.is_dir()):
        text = ''
        for f in sorted(cell.rglob('*')):
            if f.is_file() and f.suffix in ('.py', '.md', '.json') and f.stat().st_size < 400_000:
                try:
                    text += f.read_text(errors='replace')
                except Exception:
                    pass
        if not text or BLOCK.search(text):
            continue
        by_unit = defaultdict(set)
        for m in QUANT.finditer(text):
            by_unit[m.group(2).lower()].add(m.group(1))
        if not by_unit:
            continue
        purpose = ''
        for f in sorted(cell.glob('*.py')):
            m = re.search(r'"""(.{20,300}?)(?:\n\n|""")', f.read_text(errors='replace'), re.S)
            if m:
                purpose = ' '.join(m.group(1).split())[:280]
                break
        info[cell.name] = {'units': by_unit, 'purpose': purpose}

    pairs = []
    for a, b in itertools.combinations(sorted(info), 2):
        shared = set(info[a]['units']) & set(info[b]['units'])
        if len(shared) < 2:
            continue
        pairs.append((len(shared), a, b, sorted(shared)))
    pairs.sort(reverse=True)

    made = []
    for nshared, a, b, shared in pairs[:MAX_PAIRS]:
        jid = f'BT-CONN-{a}--{b}'[:80]
        d = W / 'results' / jid
        if (d / 'RESULTS.md').exists():
            continue
        d.mkdir(parents=True, exist_ok=True)
        lines = [f'# {a} and {b}', '',
                 "Two cells in the same twin. This is how they describe themselves:", '']
        for c in (a, b):
            if info[c]['purpose']:
                lines += [f'**{c}**', '', f'> {info[c]["purpose"]}', '']
        lines += [f'They carry greatness in {nshared} gemensamma dimensioner:', '']
        for u in shared[:6]:
            qa = sorted(info[a]['units'][u])[:4]
            qb = sorted(info[b]['units'][u])[:4]
            lines += [f'- **{u}** — {a}: ' + ', '.join(f'`{x}`' for x in qa)
                      + f' | {b}: ' + ', '.join(f'`{x}`' for x in qb)]
        lines += ['',
                  "Vad blir **Adjustable** if the two are connected together, which none of them can answer separately? A connection is not to write an adapter — It is that a quantity from one becomes input to the other and that something can thus be determined.", '',
                  "Write three different ways of looking at what the connection would make possible, where at least one is not about tissue or the human body at all. Drive apart first, then choose.", '',
                  "For what you choose: which quantity passes, in which direction, and what is outside us to score the result against?", '',
                  "No special answer form. DOI eller PMIDIf those two can't connect, tell me what's missing between them.", '',
                  'Status PENDING_INDEPENDENT_REVIEW.', '']
        # Same arithmetic requirement as the other generators: it costs a minute, needs no source
        # and no network, and it fell eight edges in this net on 2026-10-04.
        lines = lines + ['', "## Count on the material before searching (obligatoriskt)", '',
            "You have bash and python3. List each number with unit in the material above. Form each product and quota that gives another unit that also stands there: print × area ger kraft, styvhet × length gives force, concentration × volume provides quantity, flow × tid ger volym, hastighet × time gives length, effect × time gives energy. Compare with the number the material itself indicates and print the ratio.", '',
            "Direction determines: a quota UNDER 1 when a peak value is multiplied by its surface against a total quantity is physically impossible, and then at least one of three legs is wrong — Tell me which and why. 1 is a multiplier; name it if you can. Is all 1 within rounding, write a line about it and move on.", '',
            "In addition, for a cellPAR: if the number of the two cells gives the same quantity in different units, convert to the same unit and compare. If they differ, it is either a convention or a finding, and you should say which one."]
        (d / 'BRIEF.md').write_text('\n'.join(lines))
        (d / 'ALLOW_WEB').write_text('1\n')
        json.dump({'id': jid, 'kind': 'creative_connection_between_cells', 'cells': [a, b],
                   'shared_dimensions': shared, 'shared_dimension_count': nshared,
                   'ranking_assumption': ('pairs ranked by shared dimension count, on the assumption '
                                          'that more shared dimensions means more surface to connect '
                                          'across; that is an assumption, not a measurement'),
                   'brief_style': 'capability_framing_no_contradiction_no_interpretation',
                   'required_output_form': None, 'category': 'integration',
                   'claim_type': 'capability',
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  (d / 'JOB.json').open('w'), ensure_ascii=False, indent=1)
        made.append(jid)

    if made:
        slots = ('A', 'B', 'C', 'D')
        lines = [f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)]
        existing = QUEUE.read_text().splitlines() if QUEUE.exists() else []
        have = {l.split()[-1] for l in existing if l.split()}
        lines = [l for l in lines if l.split()[-1] not in have]
        QUEUE.write_text('\n'.join(lines + existing) + '\n')
    print(f'{len(info)} celler, {len(pairs)} par med minst 2 gemensamma dimensioner, {len(made)} briefer queued (tak {MAX_PAIRS})')
    for m in made[:8]:
        print('  ' + m)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
