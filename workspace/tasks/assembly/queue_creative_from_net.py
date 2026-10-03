"""One brief per undetermined edge, in the style that worked, with the net supplying the material.

What worked, measured. The previous batch of 69 bounded jobs returned zero quantities any cell
consumes and 15 of 15 verified-absence answers. Six briefs rewritten three times on the operator's
corrections returned 8 to 24 DOIs each with PMIDs, and the first one I consumed produced a floor under
the whole power decision: ISO 11979-2 permits the lens label to deviate by more than the 0.5 D grid step
we choose on, in 89 of 89 eyes.

The three corrections, which are the whole recipe:

1. No interpretation. Earlier I wrote "that is not noise, it is two devices agreeing with themselves
   and disagreeing with each other" and handed over my own reading, which is hard to argue with once
   it is on the page. The numbers go in bare.
2. No named analogy domains. I had written "metrology, surveying, astronomy, manufacturing control",
   which is exactly the part with the value, so naming them collects only what I already thought of.
3. Diverge before searching. Three framings first, at least one not about the subject at all, and no
   required output form. A single pointed question is a funnel.

So each undetermined edge becomes a brief carrying its own quantified gap and nothing else. The net is
assisting, which is what it is for: it already holds the pairs of numbers that cannot both be
comfortable, each with the file the number came from.

WHERE THIS FORM DOES NOT APPLY, from the dental lane's own audit and sharper than my generalisation.
The bare form is for a job that HUNTS something outside us. A job that COMPUTES a number we then
consume keeps the strict rules, because dental measured 12 of 30 such jobs failing on rigour -- a wrong
tail, and a universal quantifier taken from a grid. Those are failure modes of computation, not of
search, and removing the scaffolding there removes the only thing catching them. I was about to carry
the free form across everything on the strength of one good batch.

One detail that carries the whole result and is easy to lose when adapting this: `constraint` must
CONTAIN the number, because the brief is literally that text quoted. An edge reading "uncertain, needs
measurement" produces a worthless brief; one reading "identical surface history to 0.0 K leaves 0.3806 K
at target depth" produces a good one. Numeric gaps on the edges are the precondition, not an improvement.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

W = Path('.')
NET = W / 'CONSTRAINT_NETS.json'
QUEUE = W / 'tasks/lanes/bt_queue.txt'
BLOCK = re.compile(r'excluded_category|excluded_category|excluded_category|sinusoid|excluded_category|excluded_category|excluded_category', re.I)

ASK = "\n\nBefore you look for anything: write down **three different ways of looking at this**,\nwith at least one that has nothing to do with tissue or eyes at all. Diverge first.\n\nThen choose what you think goes furthest and pursue it. What do you find?\n\nNo particular answer format. Give DOI or PMID for what you build on. If you conclude\nthat the material above points somewhere other than the question, say so.\n\nStatus PENDING_INDEPENDENT_REVIEW.\n"


def main() -> int:
    net = json.loads(NET.read_text())['bodytwin']['tissue_constraint_net']
    queued = {l.split()[-1] for l in QUEUE.read_text().splitlines() if l.split()} \
        if QUEUE.exists() else set()
    made, skipped = [], 0
    for e in net['edges']:
        if e['status'] not in ('OPEN', 'UNKNOWN'):
            continue
        if BLOCK.search(e['constraint']) or BLOCK.search(' '.join(e['between'])):
            skipped += 1
            continue
        tag = re.sub(r'^T-E\d+-', '', e['id']).upper().replace('_', '-')[:46].strip('-')
        jid = f'BT-NET-{tag}'
        d = W / 'results' / jid
        if (d / 'RESULTS.md').exists():
            continue
        if jid in queued:
            # Idempotent on job CREATION was not enough: brief_cycle.sh runs unattended, and a job
            # whose answer has not come back yet has no RESULTS.md, so every cycle queued all twelve
            # again. Measured at 13 duplicate rows after two cycles.
            continue
        d.mkdir(parents=True, exist_ok=True)

        # The neighbourhood, as plain facts. No status vocabulary, no framing, no interpretation.
        near = [x for x in net['edges'] if x is not e
                and set(x['between']) & set(e['between'])
                and not BLOCK.search(x['constraint'])]
        body = [f"# {' and '.join(e['between'])}", '',
                "This is what our own measurements say:", '',
                f"> {e['constraint']}", '']
        if near:
            body += ["The same quantities appear here:", '']
            body += [f"> {x['constraint'][:230]}" for x in near[:3]]
            body += ['']
        (d / 'BRIEF.md').write_text('\n'.join(body).rstrip() + ASK)
        (d / 'ALLOW_WEB').write_text('1\n')
        json.dump({'id': jid, 'kind': 'creative_from_net', 'edge': e['id'],
                   'edge_status': e['status'], 'between': e['between'],
                   'brief_style': 'raw_material_no_interpretation_no_named_analogies_diverge_first',
                   'required_output_form': None, 'category': 'exploratory',
                   'claim_type': 'information_link',
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  (d / 'JOB.json').open('w'), ensure_ascii=False, indent=1)
        made.append(jid)

    # Front of the queue: the previous batch sat at position 13839 and would have waited days.
    if made:
        slots = ('A', 'B', 'C', 'D')
        lines = [f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)]
        existing = QUEUE.read_text().splitlines() if QUEUE.exists() else []
        QUEUE.write_text('\n'.join(lines + existing) + '\n')
    print(f'{len(made)} briefer out of the net; {skipped} uteslutna av filtret')
    for m in made:
        print('  ' + m)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
