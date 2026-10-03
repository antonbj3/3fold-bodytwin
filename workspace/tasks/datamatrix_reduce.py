#!/usr/bin/env python3
"""Reduce finished BT-DM packet tables into one CSV per protocol."""
import csv, hashlib, json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PACKETS=ROOT/'results'
OUT=PACKETS/'DATAMATRIX_TABLES'
MANIFEST=PACKETS/'CX-DATAMATRIX/MANIFEST.json'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    manifest=json.loads(MANIFEST.read_text())['packets']; grouped={};stats=Counter()
    for item in manifest:grouped.setdefault(item['protocol'],[]).append(item)
    for proto,items in sorted(grouped.items()):
        fields=None;rows=[]
        for item in items:
            packet=PACKETS/item['packet'];spec=json.loads((packet/'inputs/PROTOCOL.json').read_text())
            fields=spec['fields'];table=packet/'TABLE.csv'
            if not table.exists():stats['unfinished']+=1;continue
            if hashlib.sha256((packet/'inputs/raw.csv').read_bytes()).hexdigest()!=spec['raw_sha256']:
                stats['hash_fail']+=1;continue
            with table.open(newline='') as f:loaded=list(csv.DictReader(f))
            if len(loaded)!=1 or any(k not in loaded[0] for k in fields):stats['schema_fail']+=1;continue
            try:values={k:float(loaded[0][k]) for k in fields}
            except (ValueError,TypeError):stats['schema_fail']+=1;continue
            rows.append({'dataset':spec['dataset'],'unit':spec['unit'],'packet':item['packet'],'raw_source':spec['raw_source'],'raw_sha256':spec['raw_sha256'],'check_value':spec['check_value'],**values});stats['included']+=1
        path=OUT/(proto+'.csv')
        with path.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=['dataset','unit','packet','raw_source','raw_sha256','check_value']+fields)
            writer.writeheader();writer.writerows(sorted(rows,key=lambda x:(x['dataset'],x['unit'])))
    print(json.dumps(stats))
if __name__=='__main__':main()
