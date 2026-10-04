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


def main() -> int:
    if not (QUEUE.exists() and LOG.exists()):
        return 0
    starts = collections.Counter(m.group(1) for m in
                                 re.finditer(r'\] start (\S+)', LOG.read_text(errors='ignore')))
    keep, moved = [], []
    for line in QUEUE.read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        jid = parts[-1]
        n = starts.get(jid, 0)
        if n >= MIN_FAILED_STARTS and not (W / 'results' / jid / 'RESULTS.md').exists():
            moved.append((n, line))
        else:
            keep.append(line)
    if not moved:
        return 0
    QUEUE.write_text('\n'.join(keep) + '\n')
    with QUAR.open('a') as f:
        for n, line in sorted(moved, reverse=True):
            f.write(f'{date.today().isoformat()} starter={n} {line}\n')
    print(f'quarantine: {len(moved)} jobb flyttade, {sum((n for n, _ in moved))} wasted starts, queue {len(keep) + len(moved)} -> {len(keep)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
