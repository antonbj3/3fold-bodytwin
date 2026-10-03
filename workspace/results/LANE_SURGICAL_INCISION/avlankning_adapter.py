"""SURGICAL_INCISION R4 entry point for the optional reviewed model gate.

Usage from BodyTwin root after review/application:
PYTHONDONTWRITEBYTECODE=1 python results/LANE_SURGICAL_INCISION/avlankning_adapter.py requests.json --mode scenario
The native R4 empirical material/depth/injury ports remain unknown.
"""
from pathlib import Path
import argparse,hashlib,json,sys

def main():
    sibling=Path(__file__).resolve().parents[1]/'LANE_SURGICAL_SYNTHESIS'
    sys.path.insert(0,str(sibling))
    from surgical_chain.path_preference import incision_path_gate
    p=argparse.ArgumentParser()
    p.add_argument('request_json',type=Path)
    p.add_argument('--mode',choices=('evidence','scenario'),default='evidence')
    args=p.parse_args()
    result=incision_path_gate(json.loads(args.request_json.read_text()),mode=args.mode)
    result['consumer']='LANE_SURGICAL_INCISION R4'
    ports=Path(__file__).with_name('PORTS.json')
    result['r4_ports_sha256']=hashlib.sha256(ports.read_bytes()).hexdigest() if ports.exists() else None
    result['r4_transfer']='No needle effective J or scenario Gamma converted into native interface/layer work'
    print(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False))

if __name__=='__main__':main()
