#!/usr/bin/env python3
"""Put the creative families at the front of the swarm queue. Idempotent, safe to run often.

Why it runs every five minutes instead of once an hour. On 2026-10-04 the ordering was set by hand at
12:48 and the front of the queue was BT-DW48-AUTO again by 13:05, with three freshly generated briefs
sitting at position 1177. The dispatcher (tasks/lanes/bt_queue2.sh) reads the file top-down and does
not rewrite it, so the rotation comes from the refill generators prepending auto-packets when the
queue runs low. Those generators are shared with other sessions and are not touched here; the
ordering is simply reasserted often enough that an auto-packet never sits ahead of a creative brief
for more than one driver tick.

Silent unless it changes something, so it can run from the driver without filling the log.
"""
from __future__ import annotations

import pathlib
import re
import sys

CREATIVE = re.compile(r'BT-(ANOM|OBST|NET|2ND|CAP|IDEA|XSEED|CONN)')
Q = pathlib.Path(__file__).resolve().parents[2] / 'tasks/lanes/bt_queue.txt'


def main() -> int:
    if not Q.exists():
        return 0
    lines = [l for l in Q.read_text().splitlines() if l.strip()]
    cre = [l for l in lines if CREATIVE.search(l)]
    if not cre or lines[:len(cre)] == cre:
        return 0
    rest = [l for l in lines if not CREATIVE.search(l)]
    Q.write_text('\n'.join(cre + rest) + '\n')
    was = rest[0].split()[-1] if rest else '-'
    print(f'creative first: {len(cre)} av {len(lines)} (gender baseline was {was})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
