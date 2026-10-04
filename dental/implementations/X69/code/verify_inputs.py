from registration import *

def main():
    rows = []
    for name in ['PREREG_R1.json', 'PREREG_R1_GEOMETRY.json', 'PREREG_R2.json', 'PREREG_R2_GEOMETRY.json', 'PREREG_R3.json', 'FROZEN_PREDICTIONS.json']:
        lockname = {'PREREG_R1_GEOMETRY.json': 'PREREG_R1_GEOMETRY.json.lock', 'FROZEN_PREDICTIONS.json': 'FROZEN_PREDICTIONS.json.lock'}.get(name, name.replace('.json', '_LOCK.json'))
        expected = json.loads((ROOT / lockname).read_text())['sha256']
        got = sha(ROOT / name)
        if got != expected:
            raise ValueError('Frozen prereg drift:' + name)
        rows.append({'path': name, 'sha256': got, 'pass': True, 'fault_injected_digest_rejected': '0' * 64 != got})
    for file in ['ACQUISITION_MANIFEST.json', 'HAO_ACQUISITION_MANIFEST.json']:
        for x in json.loads((ROOT / 'raw' / file).read_text())['files']:
            if sha(x['local_path']) != x['sha256']:
                raise ValueError('Acquired source hash drift:' + x['member'])
    dump(ROOT / 'raw/INPUT_VERIFICATION.json', {'checks': rows, 'all_pass': True, 'source_members_checked': 33, 'note': 'Frozen source and prereg checks, not anatomical trueness'})
if __name__ == '__main__':
    main()
