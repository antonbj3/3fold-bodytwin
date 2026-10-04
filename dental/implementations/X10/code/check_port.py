"""Evaluate an explicitly bounded signed-contact JSON port; no assumed sigma."""
import argparse, json
from signed_gap import certify_port
if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('port_json')
    a = p.parse_args()
    with open(a.port_json) as f:
        port = json.load(f)
    result = certify_port(port)
    result['evidence_kind'] = port.get('evidence_kind', 'UNKNOWN')
    result['scientific_admission'] = 'PENDING_INDEPENDENT_REVIEW'
    print(json.dumps(result, indent=2, allow_nan=False))
