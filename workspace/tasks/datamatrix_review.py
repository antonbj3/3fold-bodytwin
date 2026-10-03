#!/usr/bin/env python3
"""Regenerate the coordinator's sample and unqueued priority proposal."""
import json
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results'; LANE=RESULTS/'CX-DATAMATRIX'
manifest=json.loads((LANE/'MANIFEST.json').read_text())['packets']
groups=defaultdict(list)
for item in manifest:groups[item['protocol']].append(item)
lines=["# CX-DATAMATRIX — sample for the coordinator",'',"Two complete BRIEF+FILTER per protocol. The packets contain raw data and a separate load test.",'']
queue=[]
for proto,items in sorted(groups.items()):
    lines+=['## '+proto,''];chosen=[];seen=set()
    for item in items:
        if item['dataset'] not in seen:chosen.append(item);seen.add(item['dataset'])
        if len(chosen)==2:break
    for item in chosen:
        packet=RESULTS/item['packet']
        lines+=['### '+item['packet'],'','#### BRIEF.md','',(packet/'BRIEF.md').read_text().rstrip(),'','#### FILTER.md','',(packet/'FILTER.md').read_text().rstrip(),'']
    for i,item in enumerate(sorted(items,key=lambda x:(x['dataset'],x['unit']))):
        queue.append('ABC'[i%3]+' swarm '+item['packet'])
(LANE/'SAMPLE.md').write_text('\n'.join(lines)+'\n')
(LANE/'QUEUE_PROPOSAL.txt').write_text('\n'.join(queue)+'\n')
print(json.dumps({'sample_packets':16,'queue_proposals':len(queue)}))
