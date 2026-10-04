import json
from common import ROOT, sha, check_lock

def main():
    for p in ROOT.glob('PREREG*.json'):
        check_lock(p.name)
    for p in ROOT.glob('FROZEN_PREDICTIONS*.json'):
        check_lock(p.name)
    manifest = json.loads((ROOT / 'INPUT_LOCK.json').read_text())
    for row in manifest['inputs']:
        if sha(ROOT / row['local']) != row['sha256']:
            raise ValueError('Input changed: ' + row['local'])
    print('Frozen inputs and preregistrations PASS')
if __name__ == '__main__':
    main()
