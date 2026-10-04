#!/bin/bash
# Copy the live constraint net into the tracked file. Run on EVERY tick, not hourly.
#
# 2026-10-04, found by the graph lane: the live CONSTRAINT_NETS.json is a symlink into
# generations/, and .gitignore line 85 (/*.json) excludes it, so every edit lives outside git until
# copied. The hourly sync in brief_cycle.sh kept the edge COUNT current, which is what made this
# hard to see: the tracked copy had all 158 edges and all three new ASM edges, and looked complete.
# What it did not have was the eight load-label corrections made between two hourly runs. A commit
# saying "twelve lanes" then reads as if the net work is accounted for.
#
# Per-tick instead of hourly, because a tick is where the edits happen.
set -u
cd .
python3 - <<'PY'
import json, pathlib
live = json.load(open('CONSTRAINT_NETS.json'))['bodytwin']['tissue_constraint_net']
p = pathlib.Path('data/CONSTRAINT_NET_TISSUE.json')
d = json.loads(p.read_text())
old = d['bodytwin']['tissue_constraint_net']
# Report what is ABOUT to change, so a silent no-op is distinguishable from a silent miss.
def fingerprint(net):
    return (len(net['edges']), len(net['variables']),
            sum(1 for e in net['edges'] if e.get('load_label_unsupported')),
            sum(1 for e in net['edges'] if e.get('sensitivity_note')),
            sum(1 for e in net['edges'] if e.get('evidence_relocated')))
before, after = fingerprint(old), fingerprint(live)
d['bodytwin']['tissue_constraint_net'] = live
p.write_text(json.dumps(d, indent=2, ensure_ascii=False))
names = ('edges', 'variables', 'load labels', 'sensitivity notes', 'relocated evidence')
if before == after:
    print('net synced: unchanged (' + ', '.join(f'{n} {v}' for n, v in zip(names, after)) + ')')
else:
    print('net synced: ' + ', '.join(
        f'{n} {b}->{a}' for n, b, a in zip(names, before, after) if b != a))
PY
