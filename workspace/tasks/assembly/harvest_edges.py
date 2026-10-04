"""Harvest edges from the 389 rounds that measured something and never reached the net.

The finding that forced this. The net carried 22 edges from 6 lanes, while 49 of 55 lanes had completed
389 rounds between them and appeared in it nowhere. So briefs generated from the net came out
eye-heavy, and I read that as the body being eye-heavy. It is not: it is a census of what I personally
worked on, mistaken for a map. The operator's objection -- calling the other domains unrelated while
claiming to make the twin strong -- is what exposed it.

What makes an edge, and why this can be automated at all. A round report already carries the two things
an edge needs: a quantified relation, and a file the number can be read back from. The staleness
checker's own criterion is the gate -- if the value cannot be re-read from `file :: key = value`, it is
not provenance and the candidate is dropped rather than written. That keeps the automation from
inflating the net with prose.

What this deliberately does NOT do. It does not decide that a relation is TIGHT. A harvested edge
enters as UNKNOWN with its number, because asserting a relation is settled is a judgement that needs a
control, and no automated pass has one. The number is what makes the edge useful to a brief; the status
is what a human or a lane has to earn.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

W = Path('')
NET = W / 'CONSTRAINT_NETS.json'
# The block list is not kept in the repo. tasks/assembly/excluded_terms.py says why, loads it from
# outside the tree, and matches everything if it cannot be read, so a missing list rejects rather
# than admits.
def _load_block_pattern(extra: str = '') -> 're.Pattern[str]':
    import importlib.util, pathlib
    for parent in pathlib.Path(__file__).resolve().parents:
        cand = parent / 'tasks' / 'assembly' / 'excluded_terms.py'
        if cand.exists():
            spec = importlib.util.spec_from_file_location('excluded_terms', cand)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            pat = mod.block_pattern().pattern
            return re.compile(pat + ('|' + extra if extra else ''), re.I)
    return re.compile(r'(?s).*')  # loader gone: reject everything rather than pass everything


BLOCK = _load_block_pattern()
# A number that means something: a gap, a difference, an error, a threshold, a count of violations.
INTERESTING = re.compile(r'gap|diff|error|residual|violat|threshold|margin|ratio|mismatch|undetermined'
                         r'|discrep|shortfall|excess|spread|bias', re.I)
SKIP = re.compile(r'count$|_n$|^n_|seed|version|round|elapsed|_s$|epoch|index|size|tokens', re.I)


def numeric_leaves(obj, prefix=''):
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            out += numeric_leaves(v, f'{prefix}.{k}' if prefix else k)
    elif isinstance(obj, (int, float)) and not isinstance(obj, bool):
        out.append((prefix, float(obj)))
    return out


def main() -> int:
    full = json.loads(NET.read_text())
    net = full['bodytwin']['tissue_constraint_net']
    have_ev = ' '.join(str(e.get('evidence', '')) for e in net['edges'])
    existing_ids = {e['id'] for e in net['edges']}

    candidates, dropped = [], {'no_number': 0, 'blocked': 0, 'already': 0}
    for lane_dir in sorted((W / 'results').glob('LANE_*')):
        lane = lane_dir.name
        if lane in have_ev:
            dropped['already'] += 1
            continue
        rounds = sorted(lane_dir.glob('night_rounds/r*.json'),
                        key=lambda p: int(re.sub(r'\D', '', p.stem) or 0))
        if not rounds:
            continue
        r = rounds[-1]
        try:
            d = json.loads(r.read_text())
        except Exception:
            continue
        text = json.dumps(d, ensure_ascii=False)
        if BLOCK.search(text):
            dropped['blocked'] += 1
            continue
        picks = [(k, v) for k, v in numeric_leaves(d)
                 if INTERESTING.search(k) and not SKIP.search(k.split('.')[-1])
                 and v not in (0.0, 1.0)]
        if not picks:
            dropped['no_number'] += 1
            continue
        key, val = max(picks, key=lambda kv: abs(kv[1]))
        obstacle = str(d.get('obstacle') or d.get('gate_scope') or '')[:220]
        candidates.append({
            'lane': lane, 'round': r.stem, 'key': key, 'value': val,
            'evidence': f'{r.relative_to(W)} :: {key} = {val}',
            'obstacle': obstacle,
            'operation': str(d.get('operation') or '')[:200],
        })

    added = []
    for i, c in enumerate(candidates):
        var = re.sub(r'^LANE_', '', c['lane']).lower()
        eid = f'H-E{i + 1}-{var}'
        if eid in existing_ids:
            continue
        net['variables'].append(var) if var not in net['variables'] else None
        net['edges'].append({
            'id': eid,
            'between': [var],
            'constraint': (f"{c['key']} = {c['value']} in the lane's own latest round. "
                           f"Obstacle as the lane states it: {c['obstacle']}"),
            'status': 'UNKNOWN',
            'evidence': c['evidence'],
            'provenance': 'auto_harvested_2026_10_03',
            'why_unknown': ('a harvested number is not a settled relation; status has to be earned by a '
                            'control, and this pass has none'),
        })
        added.append(eid)

    NET.write_text(json.dumps(full, indent=1, ensure_ascii=False))
    print(f'lanes without edge carrying a number: {len(candidates)}')
    print(f'  kanter tillagda: {len(added)}')
    print(f"  net now: {len(net['variables'])} variabler, {len(net['edges'])} kanter")
    print(f'  uteslutna: {dropped}')
    for c in candidates[:10]:
        print(f"    {c['lane']:32s} {c['key'][:44]:44s} = {c['value']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
