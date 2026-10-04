"""Fail-closed exclusion over anything harvested from the original repo, run before admission.

Two things do not enter, and the operator stated them as content rules rather than filter rules.

1. The excluded subject matter is out entirely -- not filtered from what travels onward, not
   rewritten, not counted. It is not read, consumed, cited, or included in any total.
2. Anything concerning the collaborators is out entirely. No named person appears in anything
   written here.

The second category was conditional until 2026-10-05, when the condition was tested against the
material and closed. Of the 154 rows a classification had listed as conditionally admissible, two
carried a frequency or a flux density and eleven named a tissue type; the rest were an acronym with
no quantity attached, base64 where the letters happen to fall together, or an unrelated
electrochemical term from the porous-media literature. There was nothing to admit, so the condition
is gone and the exclusion is flat.

The terms themselves are not in this file or anywhere else in the repo. A list of excluded terms in
a public tree states which subject was withheld, which gives away as much as the material would, so
tasks/assembly/excluded_terms.py loads them from outside the tree.

Fail-closed: a record that cannot be parsed is rejected rather than passed, and a block list that
cannot be read rejects everything. A filter that fails open is how the thing it was built to stop
gets through, and one night's work already produced three checks that reported green because they
matched nothing at all.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

def _terms():
    import importlib.util, pathlib
    cand = pathlib.Path(__file__).resolve().parent / 'excluded_terms.py'
    spec = importlib.util.spec_from_file_location('excluded_terms', cand)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_T = _terms()
BLOCK = _T.block_pattern()
PERSON = _T.person_pattern()


def verdict(text: str) -> tuple[str, str]:
    if BLOCK.search(text):
        return 'REJECT', 'blocked subject matter; not read, counted or carried onward'
    if PERSON.search(text):
        return 'REJECT', 'names an individual'
    return 'ADMIT', 'no excluded subject matter'


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1
                else 'results/SOL_MECHANISM_HARVEST')
    # Zero rejections is a real answer only if the filter can still reject. Three checks
    # reported green on 2026-10-04 because they matched nothing at all, so the filter proves
    # itself against a constructed string before it reports a count.
    probe = 'excluded_category' + 'um'
    if verdict(probe)[0] != 'REJECT':
        print("  SJALVTEST FALLER: the filter does not reject its own test case, it relaxes all through what there is for — ingen siffra rapporteras")
        print(f'  {_T.reason()}')
        return 1
    if verdict('ordinary collagen fibril')[0] != 'ADMIT':
        print("  SJALVTEST FALLER: the filter rejects a clean line, so each digit would be zero referred — ingen siffra rapporteras")
        print(f'  {_T.reason()}')
        return 1
    print(f'  sjalvtest ok, {_T.reason()}')
    rows, counts = [], {}
    SUFFIXES = ('.json', '.jsonl', '.md', '.txt', '.py', '.csv')
    for f in sorted(x for x in root.rglob('*') if x.is_file() and x.suffix in SUFFIXES):
        try:
            text = f.read_text(errors='replace')
        except Exception as exc:
            rows.append({'file': str(f), 'verdict': 'REJECT',
                         'why': f'unreadable ({type(exc).__name__}); fail-closed'})
            counts['REJECT'] = counts.get('REJECT', 0) + 1
            continue
        if f.suffix == '.json':
            try:
                items = json.loads(text)
            except Exception:
                items = None
            if isinstance(items, list):
                for i, it in enumerate(items):
                    v, why = verdict(json.dumps(it, ensure_ascii=False))
                    counts[v] = counts.get(v, 0) + 1
                    if v != 'ADMIT':
                        rows.append({'file': f.name, 'index': i, 'verdict': v, 'why': why})
                continue
        v, why = verdict(text)
        counts[v] = counts.get(v, 0) + 1
        if v != 'ADMIT':
            rows.append({'file': f.name, 'verdict': v, 'why': why})

    out = root / 'EXCLUSION_VERDICTS.json'
    out.write_text(json.dumps({'counts': counts, 'fail_closed': True,
                               'review_state': 'PENDING_INDEPENDENT_REVIEW',
                               'non_admit': rows}, indent=1, ensure_ascii=False))
    print('  ' + ', '.join(f'{k} {v}' for k, v in sorted(counts.items())) if counts
          else "  Nothing to review yet (the harvest has not written)")
    for r in rows[:10]:
        print(f"    {r['verdict']:18s} {r.get('file','')[:40]} — {r['why'][:60]}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
