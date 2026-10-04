#!/usr/bin/env python3
"""Emission rate of verdict-named fields across finished jobs — the probe a reference count cannot be.

Why this exists. The graph lane built an instrument that walks every produced dict key and asks whether
any OTHER file references it, which finds fields that are produced and never consumed. That probe has a
blind spot I ran into on our side: a verdict field that is almost never EMITTED looks consumed, because
the few records that do carry it are read normally. `negative_result` is exactly that case here — it is
read wherever it appears and is still absent from 99.6 % of results. Reference count cannot see it;
emission rate can.

So this measures the other axis. For every field whose NAME claims a verdict, it reports how many
finished jobs emit it at all, and with which types. A field we ask for in briefs and receive in 0.4 %
of results is a quality control that is not happening, whatever the brief says.

The promise filter is deliberately the same idea as the vacuous-check scanner: a name that announces a
decision. Without it the output is dominated by ordinary report fields, which are written for a person
to read and whose absence means nothing.

Two things this does NOT claim. A low emission rate is not a defect by itself: most fields are optional
and most jobs have no occasion to set them. It becomes a defect only where a brief REQUIRES the field,
and that pairing is the last column, computed by looking for the field name in the job's own BRIEF.md.
And a high emission rate says nothing about whether the value is correct.

Usage:
  python3 tasks/build_night/verdict_emission.py            # whole results tree
  python3 tasks/build_night/verdict_emission.py --min 3     # only fields emitted at least 3 times
"""
from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'results'

# A name that announces a decision. Kept close to the vacuous-scanner's promise list on purpose.
PROMISE = ('negative_result', 'refutes_us', 'falsif', 'refut', 'contradict', 'admissib', 'violat',
           'refused', 'unverified', 'degenerate', 'aliased', 'reversed', 'certif', 'abstain',
           'gate', 'guard', 'valid', 'holds', 'passed', 'failed', 'verdict', 'decision',
           'excluded', 'survives', 'identifiab')


def is_verdict(name: str) -> bool:
    low = name.lower()
    return any(word in low for word in PROMISE)


def walk_keys(obj, depth: int = 0, limit: int = 2):
    """Top two levels only. Deeper nesting is per-job structure, not a shared contract."""
    if depth > limit or not isinstance(obj, dict):
        return
    for k, v in obj.items():
        if isinstance(k, str):
            yield k, v
        if isinstance(v, dict):
            yield from walk_keys(v, depth + 1, limit)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--min', type=int, default=1, help='only report fields emitted at least this often')
    args = ap.parse_args()

    total = 0
    emit: Counter = Counter()
    types: defaultdict[str, Counter] = defaultdict(Counter)
    required: Counter = Counter()
    asked_sample: list = []

    for p in RESULTS.glob('*/results.json'):
        try:
            d = json.loads(p.read_text())
        except (OSError, ValueError, UnicodeDecodeError):
            continue
        if not isinstance(d, dict):
            continue
        total += 1
        seen = set()
        for k, v in walk_keys(d):
            if not is_verdict(k) or k in seen:
                continue
            seen.add(k)
            emit[k] += 1
            types[k][type(v).__name__] += 1
        # Deliberately NOT done here. 2/10 17:45: the first version asked, per job, whether each of
        # the accumulated field names appeared in that job's BRIEF.md, with a fresh regex each time.
        # That is O(files x distinct_fields) with the inner cost growing as the field set grows, and
        # it had produced no output after 35 minutes over ~12k jobs. The asked-but-absent question is
        # worth answering but it is a second pass over a sample, not an inner loop over everything.
        asked_sample.append(p.parent)

    # Bounded second pass: ask the asked-but-absent question on a sample, with the field names
    # compiled once. A rate from 400 jobs is enough to say whether a requirement is being ignored.
    if asked_sample:
        import random
        random.seed(0)
        # 2/10 17:50: an alternation over EVERY distinct field name was the real hot spot - the field
        # set runs to thousands, and a regex with thousands of branches over each brief is far slower
        # than parsing all 12k results files (measured: parse 7 s, key walk 10 s, this: minutes).
        # The report only prints the common fields anyway, so the question is asked about those.
        names = [k for k, _ in emit.most_common(40)]
        pattern = None
        for d in random.sample(asked_sample, min(400, len(asked_sample))):
            brief = d / 'BRIEF.md'
            if not brief.exists():
                continue
            try:
                text = brief.read_text(errors='replace')
            except OSError:
                continue
            try:
                emitted = set(json.loads((d / 'results.json').read_text()))
            except (OSError, ValueError, UnicodeDecodeError):
                emitted = set()
            for k in names:
                # Word boundaries, not substring. 2/10 17:52: plain `k in text` reported "gate" as
                # asked in 309 of 400 briefs, because it matches inside aggregate, mitigate and
                # investigate. The column is only worth printing if a short common name cannot
                # borrow counts from unrelated words.
                if k not in emitted and re.search(r'(?<![A-Za-z_])' + re.escape(k) + r'(?![A-Za-z_])', text):
                    required[k] += 1
        print(f'asked-but-absent measured on a sample of {min(400,len(asked_sample))} jobs\n')

    if not total:
        print('no results.json found')
        return 1
    print(f'finished jobs with a parseable results.json: {total}\n')
    print(f'{"field":34s} {"emitted":>8s} {"rate":>7s}  {"asked-but-absent":>16s}  types')
    for k, n in emit.most_common():
        if n < args.min:
            continue
        t = ', '.join(f'{ty}:{c}' for ty, c in types[k].most_common(4))
        print(f'{k[:34]:34s} {n:8d} {100*n/total:6.2f}%  {required.get(k,0):16d}  {t}')
    mixed = {k: dict(c) for k, c in types.items() if len([x for x in c if x != 'NoneType']) > 1}
    if mixed:
        print('\nfields with inconsistent type (a consumer cannot tell "false" from "prose explaining why"):')
        for k, c in sorted(mixed.items(), key=lambda kv: -emit[kv[0]]):
            print(f'  {k:34s} {c}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
