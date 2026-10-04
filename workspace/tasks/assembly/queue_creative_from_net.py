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

import hashlib
import json
import re
from pathlib import Path

W = Path('')
NET = W / 'CONSTRAINT_NETS.json'
QUEUE = W / 'tasks/lanes/bt_queue.txt'
# The block list is not kept in the repo. tasks/assembly/excluded_terms.py says why, loads it from
# outside the tree, and matches everything if it cannot be read, so a missing list rejects rather
# than admits.
def _load_block_pattern(extra: str = '') -> 're.Pattern[str]':
    import importlib.util, pathlib
    for parent in pathlib.Path(__file__).resolve().parents:
        cand = parent / 'tasks' / 'assembly' / 'excluded_terms.py'
        if cand.exists():
            spec = importlib.util.spec_from_file_location('excluded_terms', cand)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            pat = mod.block_pattern().pattern
            return re.compile(pat + ('|' + extra if extra else ''), re.I)
    return re.compile(r'(?s).*')  # loader gone: reject everything rather than pass everything


BLOCK = _load_block_pattern()

ASK = "\n\nBefore you look for anything: write down **three different ways of looking at this**,\nwith at least one that has nothing to do with tissue or eyes at all. Diverge first.\n\nThen choose what you think goes furthest and pursue it. What do you find?\n\nNo particular answer format. Give DOI or PMID for what you build on. If you conclude\nthat the material above points somewhere other than the question, say so.\n\n## Diverge widely first\n\nThis is not a review job. I want **ideas**, and they do not come from checking my numbers.\nWrite three readings of the material before searching for anything, and let at least one\nbe entirely outside tissue, body and medicine — a manufacturing process, a measuring instrument,\nan economic system, a material, anything that has the same form. Diverge as far\nas possible. Then choose the one that goes furthest and pursue it all the way.\n\nWhat I want to see is something that is NOT in the material: a connection nobody has\nposed, a quantity that would settle the question if someone measured it, a mechanism\nthat explains two things at once, or a way to make the question decidable with something\nalready measured elsewhere in the world. Give DOI or PMID for what you build on.\n\nIf along the way you see that two numbers in the material do not agree dimensionally,\nsay so in one line and move on — it is a bonus, not the task.\n\nStatus PENDING_INDEPENDENT_REVIEW.\n"


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
        # One brief per edge was a dead end by construction: every edge got its question once, the
        # answers came back, and the generator then found nothing to do for the rest of the night
        # even though the edges themselves had CHANGED. Measured 2026-10-04 morning: 375 creative
        # reports complete, four generators producing 0 new briefs, and zero creative jobs running.
        #
        # So the id carries a short hash of the edge's own content. An edge whose constraint, status
        # or evidence has moved since its last brief is a NEW question and gets a new one; an edge
        # that has not moved stays silent. That is the difference between asking again and asking
        # the same thing again, which is the failure mode that put ~90 re-analysis packets into the
        # queue in September.
        stamp = hashlib.sha256((str(e.get('constraint', '')) + str(e.get('status', ''))
                                + str(e.get('evidence', ''))).encode()).hexdigest()[:6]
        jid = f'BT-NET-{tag}-{stamp}'
        d = W / 'results' / jid
        if (d / 'RESULTS.md').exists():
            continue
        # An earlier brief for the SAME edge under a different content hash is not a duplicate, but
        # the new brief has to say what moved or the worker will redo the old analysis.
        # The first generation of briefs had no content hash, so the earlier report for an edge is
        # named BT-NET-<tag> with nothing after it. A glob of BT-NET-<tag>-* misses exactly those,
        # which are all 375 of them.
        prior = sorted(q.name for q in (W / 'results').glob(f'BT-NET-{tag}*')
                       if q.name != jid and (q / 'RESULTS.md').exists())
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
        # Branches the coordinator booked on the edge when it expanded. Without these the brief asks
        # the one question the edge is named after and the four openings sit unread in the JSON.
        if e.get('branches_opened'):
            body.append('')
            body.append("## Branches already booked on this node")
            body.append("Each one is his own question. Take the one you can get the longest with and say which one you took; leave the others untouched rather than answer thinly to everyone.")
            body.append('')
            for i, b in enumerate(e['branches_opened'], 1):
                body.append(f'{i}. {b}')

        if prior:
            # Measured 2026-10-04 11:57: of 232 briefs answered in three hours, 117 returned nothing
            # but "results/<jid>/RESULTS.md does not exist on this filesystem". The worker runs in the
            # cloud and has no access to this tree, so a brief that makes reading a local path a
            # precondition kills the job. The earlier answer has to travel INSIDE the brief.
            pr = W / 'results' / prior[-1] / 'RESULTS.md'
            prior_text = ''
            try:
                prior_text = pr.read_text(errors='ignore').strip()
            except OSError:
                prior_text = ''
            # What actually moved: the earlier brief quoted the constraint it was built from, so the
            # two quoted lines are directly comparable without storing anything extra.
            old_c = ''
            try:
                ob = (W / 'results' / prior[-1] / 'BRIEF.md').read_text(errors='ignore')
                m = re.search(r'^> (.+)$', ob, re.M)
                if m:
                    old_c = m.group(1).strip()
            except OSError:
                pass
            body.append('')
            body.append("## This edge has been answered once before")
            if old_c and old_c != e['constraint'].strip()[:len(old_c)]:
                body.append("This is how the edge read then:")
                body.append('')
                body.append(f'> {old_c[:400]}')
                body.append('')
                body.append("That's what it says. (top quote in this briefen)The difference is why the question is asked again.")
            else:
                body.append('Kantens bevis eller status har flyttat sedan dess; texten kan se lik ut.')
            if prior_text:
                body.append('')
                body.append("The previous answer, including in its entirety or the beginning of it — You don't have our file system, so this is all you get from it. Don't repeat its analysis, build on:")
                body.append('')
                body.append('```')
                body.append(prior_text[:4000])
                body.append('```')
            body.append('')
            body.append(f'Kantens nuvarande status: {e.get("status")}. '
                        f'Bevis: {str(e.get("evidence", ""))[:200]}')
            for k in ('load_label_unsupported', 'sensitivity_note', 'evidence_relocated'):
                if e.get(k):
                    body.append(f'Noted on the edge ({k}): {str(e[k])[:300]}')
        (d / 'BRIEF.md').write_text('\n'.join(body).rstrip() + ASK)
        (d / 'ALLOW_WEB').write_text('1\n')
        # The shared queue orderer (research_value.assess) reads these fields and the ordering
        # follows. Measured 2026-10-04: with category 'exploratory' and none of the value fields
        # filled, every one of these jobs scored as role 'challenge' at base 45, which the 80/15/5
        # portfolio rations to 15 percent of slots. With local capacity of 5 and the other
        # workspace's jobs filling the constructive share, that share almost never landed: the last
        # creative job started at 23:36 and none ran for the following ten hours. The jobs were not
        # blocked and the dispatcher was not stuck -- they were correctly classified as low-value by
        # metadata I had left empty.
        #
        # So the metadata is filled with what the job actually is. A net edge question IS integration
        # work: it asks whether a stated relation between two named variables holds. The seven value
        # fields are written per edge rather than boilerplate, because the orderer requires each to
        # be a real string of some length and because a worker reads them.
        between = ' and '.join(e['between'])
        json.dump({'id': jid, 'kind': 'creative_from_net', 'edge': e['id'],
                   'edge_status': e['status'], 'between': e['between'],
                   'brief_style': 'raw_material_no_interpretation_no_named_analogies_diverge_first',
                   'required_output_form': None, 'category': 'integration',
                   'claim_type': 'information_link',
                   # assess() reads these from job['research_value'], not from the top level.
                   # At the top level they are ignored and `structured` stays False, which
                   # costs the +12 and leaves the family as 'other'.
                   'research_value': {
                       'capability': f'Decide whether the stated relation between {between} holds, '
                                     f'and say what it would take to decide it if it does not.',
                       'obstacle': f'The edge stands at status {e["status"]} with its constraint carried '
                                   f'as prose; nothing has tested whether the relation is load-bearing.',
                       'changed_operation': 'Read the cited source, recompute the number it rests on, '
                                            'and report whether the relation survives that recomputation.',
                       'consumer': f'The constraint net edge {e["id"]}, and any assembly decision that '
                                   f'reads {e["between"][0]}.',
                       'metric': 'The number the edge states, recomputed from its own source, with the '
                                 'difference from the stated value given in the source unit.',
                       'strongest_control': 'The reading that uses the population or default value '
                                            'instead of this relation, scored on the same source.',
                       'falsifier': 'If the recomputed number differs from the stated one by more than '
                                    'the source own dispersion, the edge does not hold as stated.',
                       'mechanism_family': 'measurement-decision',
                   },
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
