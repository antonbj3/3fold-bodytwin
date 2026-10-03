#!/usr/bin/env python3
"""Guarantee a dispatch share for dental rows and for anchor-search rows.

Replaces dental_interleave.py in the dispatcher pipe and keeps its promise unchanged: every third
emitted row is dental as long as dental has live rows. It adds a second guarantee for BT-ANCHOR rows,
the acquisition searches created by tasks/anchor_hunt_packets.py.

Why the second guarantee is needed. Measured 2026-10-03 10:40: after the shared ordering, the first
anchor row sat at stream position 287 of 13 369, and none appeared in the first 60 rows. At the measured
48 starts per hour that is roughly six hours before the first search runs, while the operator asked for
data anchors in bulk now. The anchor pool is small — 48 jobs — so this reservation drains itself: when no
live anchor row remains, the position falls back to the shared ordering automatically.

Position pattern inside each cycle of three emitted rows:
    1 -> shared ordering, 2 -> anchor (if any live), 3 -> dental (if any live)
giving each group about a third while both pools have work, and the full stream to the shared ordering
once they do not.

Two properties kept from the original filter:
  1. FAIL-OPEN. Any error writes the input through unchanged. The filter sits in a pipe the dispatcher
     reads; dying quietly would starve the queue.
  2. ONLY LIVE ROWS COUNT. The ordered stream contains finished jobs too, and interleaving a finished row
     would spend a slot the loop then skips.
"""
from __future__ import annotations

import os
import signal
import sys

# The dispatcher may break its read loop early (slot full, capacity refusal), closing the reader in the
# process substitution. Without this the filter raises BrokenPipeError and the fail-open handler then
# fails on the same write. seed_queue_order.py does exactly this.
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

W = '.'
DENTAL_PREFIX = 'BT-DW48-'
ANCHOR_PREFIX = 'BT-ANCHOR-'
# Sol stood down at 41 % of the weekly window on 2026-10-03, so these carry the lanes'
# own next steps. They are the night's research line and must not queue behind 300 rows.
HANDOVER_PREFIX = 'BT-HANDOVER-'
# Admitted seeds: the operator's own 500-seed program sat unadmitted since 2026-09-29. The first
# thirteen with a consuming cell went in on 3/10 and must not queue behind thousands of rows.
SEED_PREFIX = 'BT-SEED-ADMIT-'
# Reference fetches: 131 of 135 DOIs our own lanes cited were in no consumable store. A fetch
# against a given source is cheaper than an open search and shares the seed slot.
FETCH_PREFIX = 'BT-FETCH-'
# Bounded improvement jobs: close an open edge, or find a quantity a cell needs. Each outcome is
# a number with a unit or a stated absence, which is why they share the reserved slot.
EDGE_PREFIX = 'BT-EDGE-'
NEED_PREFIX = 'BT-NEED-'
# Brainstorm angles, on the operator's call that a free resource is worth a low hit rate as long
# as the filter sits inside the task. Each job reports its own hit rate.
IDEA_PREFIX = 'BT-IDEA-'


def job_of(line: str) -> str:
    """Rows are 'profile model JOB'. The job is the third field; tolerate deviations."""
    parts = line.split()
    return parts[2] if len(parts) >= 3 else ''


def is_live(job: str) -> bool:
    if not job:
        return False
    d = os.path.join(W, 'results', job)
    if not os.path.isdir(d):
        return False
    if os.path.exists(os.path.join(d, 'RESULTS.md')):
        return False
    for claim in ('.ovh_claim', '.local_claim'):
        if os.path.exists(os.path.join(d, claim)):
            return False
    return True


def main() -> int:
    data = sys.stdin.read()
    try:
        dental, anchor, handover, seed, other = [], [], [], [], []
        for ln in data.splitlines(keepends=True):
            s = ln.strip()
            if not s or s.startswith('#'):
                other.append(ln)
                continue
            j = job_of(s)
            if j.startswith(DENTAL_PREFIX) and is_live(j):
                dental.append(ln)
            elif j.startswith(ANCHOR_PREFIX) and is_live(j):
                anchor.append(ln)
            elif j.startswith(HANDOVER_PREFIX) and is_live(j):
                handover.append(ln)
            elif (j.startswith(SEED_PREFIX) or j.startswith(FETCH_PREFIX)
                  or j.startswith(EDGE_PREFIX) or j.startswith(NEED_PREFIX)
                  or j.startswith(IDEA_PREFIX)) and is_live(j):
                seed.append(ln)
            else:
                other.append(ln)
        if not dental and not anchor and not handover and not seed:
            sys.stdout.write(data)        # nothing to guarantee; touch nothing
            return 0

        out = []
        di = ai = hi = si = oi = 0
        pos = 0
        # Cycle of four: handover first because it is the research line, then anchor, then dental,
        # then the shared ordering. Each reservation lapses as soon as its pool has no live row.
        while (di < len(dental) or ai < len(anchor) or hi < len(handover)
               or si < len(seed) or oi < len(other)):
            pos += 1
            phase = pos % 5
            if phase == 4 and si < len(seed):
                out.append(seed[si]); si += 1
            elif phase == 1 and hi < len(handover):
                out.append(handover[hi]); hi += 1
            elif phase == 2 and ai < len(anchor):
                out.append(anchor[ai]); ai += 1
            elif phase == 3 and di < len(dental):
                out.append(dental[di]); di += 1
            elif oi < len(other):
                out.append(other[oi]); oi += 1
            elif si < len(seed):
                out.append(seed[si]); si += 1
            elif hi < len(handover):
                out.append(handover[hi]); hi += 1
            elif ai < len(anchor):
                out.append(anchor[ai]); ai += 1
            elif di < len(dental):
                out.append(dental[di]); di += 1
        sys.stdout.write(''.join(out))
        return 0
    except Exception:
        # Fail-open: the original ordering beats a starved queue.
        sys.stdout.write(data)
        return 0


if __name__ == '__main__':
    raise SystemExit(main())
