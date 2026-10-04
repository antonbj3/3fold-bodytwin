#!/usr/bin/env python3
"""Move jobs that have been dispatched repeatedly without ever producing a result out of the queue.

Measured 2026-10-04 14:55: 44 job ids sitting in bt_queue.txt had been started three or more times
each with no RESULTS.md ever written, consuming 169 dispatch slots between them. One of them,
BT-DW48-AUTO-94170636bc8fa7, had been started nine times, carried an AGENT_EXIT from 2026-10-03 and
sat at the front of the queue through every reordering -- which is why the front looked frozen on the
same id for hours.

The queue prune in brief_cycle.sh only removes lines that HAVE a RESULTS.md, so a job that fails
instead of finishing is immortal: it is re-read, re-dispatched and re-failed on every pass, and it
holds the position it occupies. That is a throughput leak, not a correctness problem, and it falls on
the shared worker pool rather than on any one lane.

Nothing is deleted. Quarantined lines are appended to bt_queue.quarantine with their start count and
the date, so another session can put any of them back after fixing whatever makes them fail. The
threshold is deliberately high: three failed starts is past the point where a transient explains it.
"""
from __future__ import annotations

import collections
import pathlib
import re
import sys
from datetime import date

W = pathlib.Path(__file__).resolve().parents[2]
QUEUE = W / 'tasks/lanes/bt_queue.txt'
LOG = W / 'tasks/lanes/bt_queue.log'
QUAR = W / 'tasks/lanes/bt_queue.quarantine'
MIN_FAILED_STARTS = 3
# A row can also jam the pool without ever starting. Measured 2026-10-05 15:26: the dispatcher log
# held 15 491 lines of "local admission deferred ... status=75", the top offender deferred 207 times,
# and the job at the head of the live queue had been deferred in a tight retry loop for 13 minutes
# while only 2 of 5 worker slots were in use. A deferred row never accumulates a start, so the
# start-count rule above cannot see it. One of the worst offenders was a row this coordinator itself
# re-queued, so this is not someone else's mess.
# Deferrals arrive about every 18 s, so 20 of them is roughly six minutes of one row holding the
# head of the queue -- long enough that a transient resource dip is ruled out, short enough that the
# pool is not idle for a quarter of an hour. The first threshold of 40 missed the live offender by
# six deferrals while it blocked the pool for fourteen minutes.
MIN_DEFERRALS = 20


def main() -> int:
    if not (QUEUE.exists() and LOG.exists()):
        return 0
    text = LOG.read_text(errors='ignore')
    starts = collections.Counter(m.group(1) for m in re.finditer(r'\] start (\S+)', text))
    deferred = collections.Counter(m.group(1) for m in
                                   re.finditer(r'admission deferred (\S+)', text))
    keep, moved = [], []
    for line in QUEUE.read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        jid = parts[-1]
        n = starts.get(jid, 0)
        dn = deferred.get(jid, 0)
        done = (W / 'results' / jid / 'RESULTS.md').exists()
        if not done and (n >= MIN_FAILED_STARTS or dn >= MIN_DEFERRALS):
            moved.append((max(n, dn), line + f'  # starter={n} deferrals={dn}'))
        else:
            keep.append(line)
    if not moved:
        return 0
    QUEUE.write_text('\n'.join(keep) + '\n')
    with QUAR.open('a') as f:
        for n, line in sorted(moved, reverse=True):
            f.write(f'{date.today().isoformat()} vikt={n} {line}\n')
    print(f'quarantine: {len(moved)} jobs moved, queue {len(keep) + len(moved)} -> {len(keep)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
