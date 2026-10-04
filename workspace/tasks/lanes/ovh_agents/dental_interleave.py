#!/usr/bin/env python3
"""Guarantee one dental row every third start, without touching the shared ordering function.

Why this file exists. Measured 2026-10-03 over 441 OVH starts in 9,17 h: dental got 30 starts, thus
6,8 %, against bodytwin's 351 = 79,6 %. Availability is not the cause — the shared queue has thousands of
dental rows — and my dispatcher already has a CEILING on dental (12 concurrent) but no FLOOR. What
determines the distribution is order() in the shared research_value.py, and I do not touch that file: dental
owns the impact term that lowered their rows. This filter sits AFTER ordering and only moves
positions, never content.

Dental explicitly requested this (3/10) after Anton asked for a reasonable distribution.

Two things make the filter safe:
  1. FAIL-OPEN. On every error, input is written out unchanged. The filter sits in a pipe the dispatcher
     reads; if it died silently, the queue would starve, so no exception may reach the top.
  2. ONLY LIVE ROWS COUNT. The ordered stream contains all queue entries, including completed ones. To
     interleave a completed dental row would cost a slot position that the loop then skips, so
     rows with RESULTS.md or an existing claim do not count as dental availability.

Order within each project is preserved, so the shared research valuation still applies — the
only change is the interleaving between projects.
"""
from __future__ import annotations

import os
import signal
import sys

# The dispatcher can BREAK the while loop early (full slot, capacity rejection), closing the reader in
# the process substitution. Without this the filter gets BrokenPipeError and spews a traceback into the dispatcher's
# stderr — and the fail-open handler in turn crashes on the same write. The sibling script
# seed_queue_order.py does exactly this; the same pattern here.
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

W = ''
EVERY = 3          # every third emitted row should be dental as long as dental has live rows
DENTAL_PREFIX = 'BT-DW48-'


def job_of(line: str) -> str:
    """The rows are 'profile model JOB'. The job is the third field; tolerate deviations."""
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
        lines = data.splitlines(keepends=True)
        dental, other = [], []
        for ln in lines:
            s = ln.strip()
            if not s or s.startswith('#'):
                other.append(ln)
                continue
            j = job_of(s)
            (dental if j.startswith(DENTAL_PREFIX) and is_live(j) else other).append(ln)
        if not dental:
            sys.stdout.write(data)        # nothing to guarantee; touch nothing
            return 0
        out = []
        di = oi = 0
        pos = 0
        while di < len(dental) or oi < len(other):
            pos += 1
            take_dental = (pos % EVERY == 0) and di < len(dental)
            if take_dental:
                out.append(dental[di]); di += 1
            elif oi < len(other):
                out.append(other[oi]); oi += 1
            elif di < len(dental):
                out.append(dental[di]); di += 1
        sys.stdout.write(''.join(out))
        return 0
    except Exception:
        # Fail-open: prefer the original order to a starved queue.
        sys.stdout.write(data)
        return 0


if __name__ == '__main__':
    raise SystemExit(main())
