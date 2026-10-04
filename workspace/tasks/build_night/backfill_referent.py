#!/usr/bin/env python3
"""Append the external-facit requirement to the BRIEFS of queued BodyTwin jobs that predate it.

Why this exists, measured rather than assumed. The requirement reaches a job through
free_controller._normalize, which runs in accepted() and therefore only touches jobs accepted after
the wrapper was deployed at 2026-10-02 05:11. On a random sample of 400 queued BT-FW jobs, 16 of 16
accepted after that time carry the requirement and 33 of 384 accepted before it do. Those jobs will
be run, will cost compute, and cannot produce a result comparable to anything outside this project.
Deleting them would throw away real identifiability work; asking them for an external facit costs
nothing, because an absent field already scores zero in the scorer and an honest `synthetic_only` is
an accepted answer.

It writes BRIEF.md and NOT JOB.json. That is the whole point of the file and it was nearly got
wrong: the dental controller composes BRIEF.md from the body once, at creation (controller.py:223),
and local_worker.sh then runs the agent with the instruction "Read BRIEF.md and execute the bounded
task". Nothing re-reads JOB.json at dispatch, so editing the body would have changed a field no
agent ever opens.

Deliberately narrow, because the queue is shared with another lane and a running agent must never
have its brief rewritten underneath it:
  * only ids that appear in the queue file and start with BT-FW, never BT-DW, which is not ours,
  * never a job with a result (results.json, or a non-empty RESULTS.md),
  * never a job under a live claim; a claim whose pid is gone is stale and does not protect it,
  * never a brief that already asks for the field.
A job that started and exited WITHOUT a result is a retry candidate and is included on purpose: 346
of them exited 75 on a provider rate limit having produced nothing, and the retry re-reads the
brief, so this is exactly when the requirement can still be added.

The appended text is free_controller._SN_REFERENT_BLOCK, read from that file rather than copied, so
backfilled jobs are asked for exactly the same field as live ones.

Usage:
  python3 tasks/build_night/backfill_referent.py              # dry run, prints what would change
  python3 tasks/build_night/backfill_referent.py --apply      # writes, after tarring the originals
"""
from __future__ import annotations

import os
import re
import sys
import tarfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
QUEUE = ROOT / 'tasks' / 'lanes' / 'bt_queue.txt'
RESULTS = ROOT / 'results'
CONTROLLER = ROOT / 'tasks' / 'free_controller.py'


def referent_block() -> str:
    """Read the literal out of the controller instead of importing it.

    Importing free_controller pulls in the shared dental controller, which has side effects on
    module load; this file must be safe to run while that controller is live.
    """
    src = CONTROLLER.read_text()
    m = re.search(r"_SN_REFERENT_BLOCK = \((.*?)\)\ndef _normalize", src, re.S)
    if not m:
        sys.exit('could not find _SN_REFERENT_BLOCK in free_controller.py; refusing to guess the text')
    block = eval('(' + m.group(1) + ')')  # a string-literal concatenation, nothing else
    if not isinstance(block, str) or 'external_referent' not in block:
        sys.exit('extracted block is not the expected requirement text')
    return block


def queued_ids() -> list[str]:
    ids = []
    for line in QUEUE.read_text().splitlines():
        parts = line.split()
        if len(parts) > 2 and not line.startswith('#') and parts[2].startswith('BT-FW'):
            ids.append(parts[2])
    return list(dict.fromkeys(ids))


def claim_is_live(d: Path) -> bool:
    """A claim protects a job only while its process exists; the dispatcher prunes stale ones too."""
    if (d / '.ovh_claim').exists():
        return True
    p = d / '.local_claim'
    if not p.exists():
        return False
    try:
        os.kill(int(p.read_text().strip()), 0)
        return True
    except (OSError, ValueError):
        return False


def has_result(d: Path) -> bool:
    try:
        if (d / 'results.json').exists():
            return True
        f = d / 'RESULTS.md'
        return f.exists() and f.stat().st_size > 0
    except OSError:
        return True  # unreadable state counts as finished rather than risk a rewrite


def main() -> int:
    apply = '--apply' in sys.argv
    block = referent_block()
    todo: list[Path] = []
    counts = {'no_dir': 0, 'no_brief': 0, 'has_result': 0, 'claimed_live': 0,
              'brief_already_asks': 0, 'unreadable': 0, 'to_change': 0}

    for jid in queued_ids():
        d = RESULTS / jid
        if not d.is_dir():
            counts['no_dir'] += 1
            continue
        brief = d / 'BRIEF.md'
        if not brief.exists():
            counts['no_brief'] += 1
            continue
        if has_result(d):
            counts['has_result'] += 1
            continue
        if claim_is_live(d):
            counts['claimed_live'] += 1
            continue
        try:
            text = brief.read_text()
        except (OSError, UnicodeDecodeError):
            counts['unreadable'] += 1
            continue
        if 'external_referent' in text:
            counts['brief_already_asks'] += 1
            continue
        counts['to_change'] += 1
        todo.append(brief)

    for k, v in counts.items():
        print(f'{k:20s} {v}')
    if not apply:
        print(f'\nDRY RUN. {len(todo)} BRIEF.md would change. Re-run with --apply.')
        for p in todo[:5]:
            print('  e.g.', p.parent.name)
        return 0

    stamp = time.strftime('%Y%m%d_%H%M%S')
    backup = ROOT / 'tasks' / 'build_night' / f'backfill_referent_backup_{stamp}.tar.gz'
    with tarfile.open(backup, 'w:gz') as tar:
        for p in todo:
            tar.add(p, arcname=f'{p.parent.name}/BRIEF.md')
    print(f'originals tarred to {backup} ({backup.stat().st_size} bytes)')

    written = skipped = 0
    for p in todo:
        d = p.parent
        # Re-check immediately before writing: the dispatcher may have claimed or finished the job
        # while this ran, and a brief must never be rewritten under a running agent.
        if has_result(d) or claim_is_live(d):
            skipped += 1
            continue
        tmp = p.with_suffix('.md.tmp')
        tmp.write_text(p.read_text() + block + '\n')
        tmp.replace(p)
        written += 1
    print(f'written {written}, skipped because they started or finished mid-run {skipped}')
    manifest = ROOT / 'tasks' / 'build_night' / f'backfill_referent_manifest_{stamp}.txt'
    manifest.write_text('\n'.join(p.parent.name for p in todo) + '\n')
    print(f'manifest {manifest}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
