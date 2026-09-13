"""Compare explicit gate-contract verdicts with frozen adapter gate vectors."""
import argparse
import hashlib
import json
from pathlib import Path
from dataclasses import asdict
from bodytwin.framework.numerical_gate_contract_v1 import GateContract,selftest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('records',nargs='+',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();results=[];hashes={}
    for leg in (0,1):
        rows=[]
        for path in args.records:
            raw=path.read_bytes();hashes[path.name]=hashlib.sha256(raw).hexdigest()
            report=json.loads(raw)
            for pair in report['rows']:
                for record in pair:
                    values=record['gate_values'];names=tuple(f'gate_{i:03}' for i in range(len(values)))
                    verdict=GateContract(names,names,'accepted').assess({'gates':dict(zip(names,values)),'accepted':record['required_pass']})
                    rows.append(dict(module=record['module'],verdict=asdict(verdict),expected=record['required_pass']))
        results.append(rows)
    gates=dict(selftest=all(selftest().values()),full_repeat=results[0]==results[1],
               verdict_parity=all(r['verdict']['accepted']==r['expected'] for r in results[0]),
               negative_retained=any(not r['expected'] and r['verdict']['required_passed']==12 and r['verdict']['required_total']==19 for r in results[0]))
    report=dict(rows=results,gates=gates,input_hashes=hashes,scope='Parity with archived typed-adapter gate vectors; ordinal keys do not revalidate original source key inventories or biological claims.')
    args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(gates=gates,rows_per_repeat=len(results[0]))))
    return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
