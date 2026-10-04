#!/usr/bin/env python3
"""Re-queue the briefs that died on a precondition the worker could not satisfy.

Measured 2026-10-04 12:05 over the 231 briefs answered in the preceding three hours: **2** came back
with nothing but "results/<jid>/RESULTS.md does not exist on this filesystem", not the 104 a first
pass claimed. The first pass searched the whole head of the report, and the phrase appears inside
long reports that note the missing path and then do the work anyway (median report 17.6 kB, the two
dead ones 4.8 and 6.1 kB). The test below therefore requires BOTH a short report and the refusal in
its opening, and that pair is what separates a dead brief from a working one.

This script takes each such job, pastes the earlier report INTO the brief in place of the pointer,
and queues it again under a new id. It does not touch jobs that did real work, and it refuses to
re-queue a job whose earlier report it cannot read either -- in that case the pointer gets dropped
and the question is asked without the history, which is still an answerable question.
"""
import pathlib
import re
import sys

W = pathlib.Path(__file__).resolve().parents[2]
QUEUE = W / 'tasks/lanes/bt_queue.txt'
# The phrases as a report use it vaguely. They stood in Swedish and English mixed
# Because the reports did; now the reports are English, and the monster is flying with.
# Monsters who don't match the next tiger in the stables to fail, and then koas nothing about.
REFUSAL = re.compile(r'(does not exist|is \*\*not\*\* present|not present|not reachable|'
                     r'cannot be satisfied|could not be read|was not found)')
POINTER = re.compile(r'`results/(BT-[A-Za-z0-9_.-]+)/RESULTS\.md`')
# The headings written by it have the file sjalv, and it must then be able to find again.
# stod tidigare pa svenska i bada andar; nar briefarna oversattes maste monstret folja med,
# otherwise it doesn't match anything and the jobs are rearranged without the previous report being attached.
SECTION = re.compile(r'\n## (?:This is not the first question about this edge|'
                     r'This edge has been answered once before)\n.*?(?=\n## |\Z)', re.S)


def main() -> int:
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    made = []
    for d in sorted((W / 'results').glob('BT-*')):
        res, brief = d / 'RESULTS.md', d / 'BRIEF.md'
        if not res.exists() or not brief.exists():
            continue
        full = res.read_text(errors='ignore')
        head = full[:1500]
        if len(full) >= 8000 or not REFUSAL.search(head):
            continue
        b = brief.read_text(errors='ignore')
        m = POINTER.search(b) or POINTER.search(head)
        if not m:
            continue
        prior = W / 'results' / m.group(1) / 'RESULTS.md'
        text = prior.read_text(errors='ignore').strip() if prior.exists() else ''
        # The pointer section is replaced whether or not the prior text was found: leaving it in is
        # what killed the job the first time.
        if text:
            repl = ('\n## This edge has been answered once before\n'
                    'The earlier answer follows pasted below. You do not have our filesystem, so this '
                    'is all you get of it. Do not redo its analysis, build on it:\n\n```\n'
                    + text[:4000] + '\n```\n')
        else:
            repl = ''
        # repl holds the earlier report verbatim, which contains backslashes; a plain
        # replacement string would be read as escapes and raise.
        nb, n = SECTION.subn(lambda _m: repl, b)
        if not n:
            # Older briefs carried the pointer as a loose paragraph rather than a section.
            nb = '\n'.join(l for l in b.splitlines() if 'RESULTS.md`' not in l) + repl
        jid = d.name.split('--')[0] + '--r2'
        nd = W / 'results' / jid
        if (nd / 'RESULTS.md').exists():
            continue
        nd.mkdir(parents=True, exist_ok=True)
        (nd / 'BRIEF.md').write_text(nb.rstrip() + '\n')
        (nd / 'ALLOW_WEB').write_text('1\n')
        for f in ('JOB.json', 'EGRESS_ROUTE.json', 'MODEL_ROUTE.json'):
            if (d / f).exists() and not (nd / f).exists():
                (nd / f).write_text((d / f).read_text(errors='ignore').replace(d.name, jid))
        made.append(jid)
        if len(made) >= limit:
            break
    if made:
        slots = ('A', 'B', 'C', 'D')
        lines = [f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)]
        existing = QUEUE.read_text().splitlines() if QUEUE.exists() else []
        QUEUE.write_text('\n'.join(lines + existing) + '\n')
    print(f'{len(made)} briefs requeued with earlier answers pasted in')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
