import json
from common import ROOT, sha

def check():
    for x in json.loads((ROOT / 'FROZEN_INPUT_MANIFEST.json').read_text())['inputs']:
        if sha(x['path']) != x['sha256']:
            raise ValueError('frozen source identity changed; preserve prior round and preregister a new input: ' + x['path'])
if __name__ == '__main__':
    check()
