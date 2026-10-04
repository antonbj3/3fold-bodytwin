"""Print the reviewed observation contract; no outcomes or target prediction supplied."""
import argparse
import json
from .capabilities import fixture, verify_freeze

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('protocol',choices=['R4P','R4C','R4F','CROWN_DIAG'])
    a=p.parse_args();verify_freeze()
    if a.protocol=='CROWN_DIAG':
        result=fixture('crown_diagnosis.json')
    else:
        result=next(x for x in fixture('lab_protocols.json')['protocols'] if x['id']==a.protocol)
    print(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False))

if __name__=='__main__':main()
