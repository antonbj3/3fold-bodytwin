#!/usr/bin/env python3
"""Open brainstorming on the free tier, with the consumability filter inside the task.

The operator's point, which I conceded: the resource is free, so a low hit rate is still worth having as
long as the output is filtered before it costs anything downstream. My earlier objection was about base
rate rather than value, and the evidence actually supports him — the 500-entry brainstorm produced 40
candidates that passed a strict consumability test, and the four false premises among the first ten
admitted were caught cheaply by the jobs themselves. The filter worked; what was wrong was treating
unfiltered output as if it carried weight.

So the only change here is that the filter moves INTO the task. A proposal is not accepted unless it
names a quantity with a unit that one of our 45 cells already computes, or an external measurement with
a locator. Everything else is returned as a rejected candidate with the reason, which is cheap and
honest, and the batch reports its own hit rate so the base rate becomes measured rather than argued.

Each job gets a different angle so the batch does not converge on one idea restated. The angles come
from where the night's measurements left open ground, not from a generic list.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('.')
QUEUE = W / 'tasks/lanes/bt_queue.txt'
RESULTS = W / 'results'

ANGLES = [
    ('DECISIONS', 'Which other decisions could the twin make, rather than predictions it could report? '
                  'The chain already chooses an implant power and a toric cylinder and is scored against '
                  'what happened to real patients. Name decisions of the same shape in other tissue: an '
                  'inverse question, a target, and an external referent that already exists.'),
    ('ESTIMATED', 'Where does practice use a POPULATION ESTIMATE for something that can be measured in '
                  'the individual? One such swap cut an error from 0.626 to 0.282 D. Name others, each '
                  'with the estimated quantity, its unit, and who measures it individually.'),
    ('CEILING', 'Where is our target precision finer than the measurement can resolve? One instrument '
                'repeats worse than the clinical threshold in 43 of 300 eyes. Name other quantities '
                'where the instrument spread is published and may exceed what we aim for.'),
    ('LOCALITY', 'Which of our summaries throw away a location, a time or a co-factor? A summary '
                 'suffices only when what the readout asks for sits at the same place, instant and '
                 'weighting. Name candidates and the coordinate each discards.'),
    ('UNUSED', 'What do we already hold and never use? Measured examples: 131 cited references nobody '
               'harvested, a declared constraint-net slot left empty, 500 candidate tracks unadmitted. '
               'Name more of this shape — identified but not consumed — with the path.'),
    ('EXTERNAL', 'What becomes askable now that an independent finite-element solver is installed? Six '
                 'interface items are required to make it answer our question and two we cannot supply. '
                 'Name questions it could settle that we cannot settle analytically.'),
    ('NEGATIVE', 'Which of our claims would be most damaging if wrong, and what single measurement would '
                 'refute it? One of our values was 3.33x a held-out measurement and its constant turned '
                 'out to have been chosen to produce a target. Name our other vulnerable values.'),
    ('ROBOT', 'What does a tool-tissue contact model need that tissue mechanics can supply? Zero of 202 '
             'robot ports carry a computed number and 109 have no source. Name quantities with units '
             'that would turn a port into a computed one.'),
]

FILTER = (
    '\n## The filter, and it is part of the task\n'
    'A proposal counts only if it carries ONE of these:\n'
    '  (a) a quantity with a unit that one of our 45 cells under `tasks/free48/sources/` already '
    'computes — name the cell and the variable; or\n'
    '  (b) an external measurement with a locator (DOI, PMID or dataset URL) and its unit.\n\n'
    'Anything else is a REJECTED candidate. Reject it yourself, in a list, with the reason in one line. '
    'A rejected candidate costs nothing; an unanchored proposal that travels downstream costs a reader '
    'their time and eventually a wrong number.\n\n'
    '**Report your own hit rate**: proposals considered, accepted under (a), accepted under (b), '
    'rejected. That makes the base rate for this kind of work measured rather than argued, which is the '
    'reason this batch exists.\n\n'
    '## Deliver\n'
    '| proposal | anchor type | quantity | unit | cell or locator | what it would change |\n'
    '|---|---|---|---|---|---|\n\n'
    'Five to fifteen accepted proposals is a good outcome; one well-anchored proposal beats ten ideas. '
    'Do not pad.\n\nStatus `PENDING_INDEPENDENT_REVIEW`. No claim of biological validation. Write '
    '`RESULTS.md` starting with the job id, plus `results.json`.\n'
)


def main() -> int:
    made = []
    for tag, question in ANGLES:
        jid = f'BT-IDEA-{tag}'
        d = RESULTS / jid
        if (d / 'RESULTS.md').exists():
            continue
        d.mkdir(parents=True, exist_ok=True)
        (d / 'BRIEF.md').write_text(
            f"# Improvement proposals, one angle: {tag}\n\n"
            f"{question}\n\n"
            f"You are looking for ADDITIONS to a tissue twin that already makes two surgical decisions "
            f"scored against 89 real patients. Read `RESULTS_2026-10-03.md` in the workspace root first "
            f"for what is established, so you do not propose what exists.\n"
            + FILTER)
        (d / 'ALLOW_WEB').write_text('1\n')
        json.dump({'id': jid, 'kind': 'improvement_brainstorm', 'angle': tag,
                   'category': 'exploratory', 'claim_type': 'information_link',
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  (d / 'JOB.json').open('w'), ensure_ascii=False, indent=1)
        made.append(jid)

    slots = ('A', 'B', 'C', 'D')
    if made:
        with QUEUE.open('a') as f:
            f.write('\n'.join(f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)) + '\n')
    print(f'queued {len(made)} brainstorm angles: {", ".join(m.split("-")[-1] for m in made)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
