#!/usr/bin/env python3
"""Re-load deterministic data-matrix samples and enforce the packet filter."""
import argparse, hashlib, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import datamatrix as dm

root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'results/CX-DATAMATRIX/MANIFEST.json').read_text())['packets']
parser=argparse.ArgumentParser();parser.add_argument('--per-protocol',type=int,default=5);args=parser.parse_args()
groups=defaultdict(list)
for x in manifest:groups[x['protocol']].append(x)
report=[]
for proto,items in sorted(groups.items()):
    # Round-robin datasets to expose format changes.
    byset=defaultdict(list)
    for item in items:byset[item['dataset']].append(item)
    selected=[]
    while len(selected)<args.per_protocol and any(byset.values()):
        for dataset in sorted(byset):
            if byset[dataset] and len(selected)<args.per_protocol:selected.append(byset[dataset].pop(0))
    for item in selected:
        path=root/'results'/item['packet'];spec=json.loads((path/'inputs/PROTOCOL.json').read_text());raw=(path/'inputs/raw.csv').read_bytes();errors=[]
        if hashlib.sha256(raw).hexdigest()!=spec['raw_sha256']:errors.append('SHA')
        try:
            fields,check=dm.measure(proto,*dm.raw_table(raw))
            if list(fields)!=spec['fields'] or not np.isfinite(check):errors.append('fields/check')
            if check>.05:errors.append('split-sample > 0.05 SD')
        except Exception as exc:errors.append(str(exc));check=None
        for f in ['BRIEF.md','FILTER.md','DATA_SUFFICIENCY.md','inputs/reader.py','inputs/NIGHT_PREAMBLE.md']:
            if not (path/f).exists():errors.append('missing '+f)
        report.append({'packet':item['packet'],'protocol':proto,'dataset':item['dataset'],'status':'PASS' if not errors else 'FAIL','check':check,'errors':errors})
output=root/'results/CX-DATAMATRIX/CHECK_PACKETS.json';output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'tested':len(report),'passed':sum(x['status']=='PASS' for x in report),'failed':sum(x['status']=='FAIL' for x in report)}))
sys.exit(any(x['status']=='FAIL' for x in report))
