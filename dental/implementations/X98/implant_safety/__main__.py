import argparse
import json
from pathlib import Path
from .module import SafetyModule, ContractError

def main():
    p = argparse.ArgumentParser(description='Research digital canal clearance with explicit physical unknowns')
    p.add_argument('--input', type=Path, required=True, help='JSON query; see examples/query.json')
    p.add_argument('--bundle-root', type=Path)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    try:
        result = SafetyModule(a.bundle_root).query(**json.loads(a.input.read_text()))
    except (ContractError, KeyError, ValueError) as e:
        p.exit(2, str(e) + '\n')
    text = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + '\n'
    if a.output:
        a.output.parent.mkdir(exist_ok=True, parents=True)
        a.output.write_text(text)
    else:
        print(text, end='')
if __name__ == '__main__':
    main()
