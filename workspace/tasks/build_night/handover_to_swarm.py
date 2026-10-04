#!/usr/bin/env python3
"""Move each waiting lane's next step from Sol to the swarm as a runnable job.

Why. The operator's rule is that above 40 % of the lane_runner weekly window, Sol stands down and the work
continues on the free models. The window reached 41.0 % at 11:00 on 2026-10-03, so tasks/build_night/STOP
is in place and no new Sol round starts. The lanes themselves do not stop: each one has a NEXT_ROUND.md
written by the round that just finished, and that is exactly a brief.

The timing is good. Two capacity fixes landed the same hour: the admission guard was counting 8.6 GiB of
reclaimable page cache as occupied memory, and reclaiming it took concurrency from 12 to 17; and
ling-3.1-flash joined swarm and swarm_worker, which matters because provider rate limits are per model and
12.5 % of fetches were returning them. So the swarm has more room than it had all night.

What a handover job carries: the lane's own next step, the control it must beat, the falsifier, and the
numbers the lane already holds. What it must not carry: anything internal. The brief is written from the
lane's files, and a job is skipped if its text matches the private filter.
"""
from __future__ import annotations

import json
import os
import re
import sys

W = ''
SLOTS = os.path.join(W, 'tasks/build_night/slots.txt')
STATE = os.path.join(W, 'tasks/build_night/state')
RESULTS = os.path.join(W, 'results')
QUEUE = os.path.join(W, 'tasks/lanes/bt_queue.txt')

def _load_block_pattern(extra: str = '') -> 're.Pattern[str]':
    """Admission pattern for the harvest. Configuration is not part of the tree."""
    import importlib.util, os, pathlib
    cand = pathlib.Path(os.environ.get('BODYTWIN_ADMISSION',
                                       pathlib.Path.home() / '.bodytwin' / 'tools' / 'admission.py'))
    if cand.exists():
        spec = importlib.util.spec_from_file_location('admission', cand)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return re.compile(mod.block_pattern().pattern + ('|' + extra if extra else ''), re.I)
    return re.compile(r'(?s).*')


BLOCK = _load_block_pattern(r'device dosimetry')
MAX_BRIEF = 9000


def waiting_lanes() -> list[str]:
    """A lane that ended and has no .go is one the coordinator has not refilled."""
    out = []
    for lane in (l.strip() for l in open(SLOTS) if l.strip()):
        if os.path.exists(os.path.join(STATE, lane + '.ended')) and \
           not os.path.exists(os.path.join(STATE, lane + '.go')):
            out.append(lane)
    return out


def read(path: str, limit: int = 6000) -> str:
    try:
        with open(path, encoding='utf-8', errors='replace') as f:
            return f.read(limit)
    except OSError:
        return ''


def last_round(lane: str) -> dict:
    d = os.path.join(RESULTS, lane, 'night_rounds')
    if not os.path.isdir(d):
        return {}
    rounds = []
    for fn in os.listdir(d):
        m = re.fullmatch(r'r(\d+)\.json', fn)
        if m:
            rounds.append((int(m.group(1)), fn))
    if not rounds:
        return {}
    rounds.sort()
    try:
        return json.load(open(os.path.join(d, rounds[-1][1])))
    except Exception:
        return {}


def brief(lane: str, nxt: str, rnd: dict) -> str:
    obstacle = str(rnd.get('obstacle') or '').strip()
    control = str(rnd.get('control') or '').strip()
    claim = str(rnd.get('claim_type') or '').strip()
    gate = str(rnd.get('gate') or '').strip()
    parts = [
        f'# {lane} — continue this research line',
        '',
        'This line was being run by a stronger model which has stood down for quota reasons. You are '
        'picking it up mid-stream, so read the state below before deciding anything. Do the next step, '
        'not a survey of the field.',
        '',
        '## The next step, as the previous round defined it',
        nxt.strip()[:4000] or '(no NEXT_ROUND.md found; derive the next step from the state below)',
        '',
        '## State you inherit',
    ]
    if claim:
        parts.append(f'- **Claim type:** `{claim}`. An *algorithm* claim needs an equally informed '
                     f'control; an *information_link* needs practice WITHOUT the new information; a '
                     f'*capability* needs an answer plus an external facit. Use the right one.')
    if gate:
        parts.append(f'- **Last gate:** `{gate}`. Do not report a gate as passed when its scope '
                     f'excludes the physical question.')
    if control:
        parts.append(f'- **Control to beat:** {control[:600]}')
    if obstacle:
        parts.append(f'- **Obstacle the last round hit:** {obstacle[:600]}')
    parts += [
        '',
        '## Rules that decide whether your result counts',
        '1. **One decisive number, with its unit**, and say what it is measured against.',
        '2. **Name your control explicitly.** A result without a stated control is not a result here.',
        '3. **State a falsifier before you compute**, and report it even when it fires. A fired '
        'falsifier is a finding, not a failure.',
        '4. **If a summary statistic decides something, test whether it is sufficient**: construct two '
        'states with the identical summary to machine precision and measure the downstream difference. '
        'This test failed in six independent subsystems on 2–3 October, so it is not a formality.',
        '5. **A missing measurement leaves as an orderable item** — quantity, unit, what it decides — '
        'not as a sentence. Write it to `ACQUISITION_TARGETS.json` in your job directory.',
        '6. **Declare the validity range** of every number you take from a source: species, tissue, '
        'temperature, method, n. A value used outside its range is the most common way a record is '
        'unusable.',
        '7. Status `PENDING_INDEPENDENT_REVIEW`. No claim of biological validation, no "breakthrough" '
        'without an equally informed control.',
        '',
        'Write `RESULTS.md` starting with the job id line, plus `results.json` with your numbers.',
    ]
    return '\n'.join(parts)[:MAX_BRIEF]


def main() -> int:
    lanes = waiting_lanes()
    if not lanes:
        print('no waiting lanes')
        return 0
    made, blocked = [], []
    for lane in lanes:
        base = os.path.join(RESULTS, lane)
        nxt = read(os.path.join(base, 'NEXT_ROUND.md'))
        rnd = last_round(lane)
        text = brief(lane, nxt, rnd)
        if BLOCK.search(text):
            blocked.append(lane)
            continue
        jid = 'BT-HANDOVER-' + re.sub(r'^LANE_', '', lane)[:40]
        d = os.path.join(RESULTS, jid)
        if os.path.exists(os.path.join(d, 'RESULTS.md')):
            continue
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, 'BRIEF.md'), 'w').write(text)
        json.dump({'id': jid, 'kind': 'sol_handover', 'category': 'mechanism',
                   'lane': lane, 'claim_type': rnd.get('claim_type'),
                   'breakthrough_priority': True,
                   'review_state': 'PENDING_INDEPENDENT_REVIEW'},
                  open(os.path.join(d, 'JOB.json'), 'w'), ensure_ascii=False, indent=1)
        made.append(jid)

    slots = ('A', 'B', 'C', 'D')
    if made:
        with open(QUEUE, 'a') as f:
            f.write('\n'.join(f'{slots[i % 4]} swarm {j}' for i, j in enumerate(made)) + '\n')
    print(f'handover jobs: {len(made)}   blocked by private filter: {len(blocked)}')
    for j in made:
        print('  ' + j)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
