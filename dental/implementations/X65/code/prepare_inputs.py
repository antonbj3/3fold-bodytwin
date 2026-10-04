"""Select a few primary sources and a stress profile. Never unpack the archive."""
from dental_release.paths import expand as _release_expand
import csv, hashlib, json, shutil, subprocess, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SRC = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))
ARCHIVE = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/usb-stage/dental_archive/scratch_K1b.tar.zst'))
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    t = time.perf_counter()
    out = ROOT / 'inputs'
    (out / 'primary').mkdir(parents=True, exist_ok=True)
    copies = {'anchors_legacy.csv': SRC / 'results/K1b_iso14801_assembly/anchors.csv', 'k1b_eval_legacy.json': SRC / 'results/K1b_iso14801_assembly/k1b_eval.json', 'run_k1b_eval_legacy.py': SRC / 'results/K1b_iso14801_assembly/run_k1b_eval.py', 'fe_summary_legacy.json': SRC / 'results/K1_lpbf_implant_fatigue/fe_summary.json', 'literature_candidates.jsonl': SRC / 'tasks/swarm48/sources/LIT_ISO14801/iso14801_papers.jsonl'}
    manifest = []
    for (name, p) in copies.items():
        dst = out / name
        shutil.copyfile(p, dst)
        manifest.append({'path': str(dst.relative_to(ROOT)), 'source': str(p), 'sha256': sha(dst), 'bytes': dst.stat().st_size})
    ids = sorted({r['pmcid'] for r in csv.DictReader((out / 'anchors_legacy.csv').open())})
    ids += ['PMC10560161', 'PMC10820107', 'PMC10976688']
    members = set(subprocess.check_output(['tar', '--zstd', '-tf', str(ARCHIVE)], text=True).splitlines())
    for pmc in ids:
        dst = out / 'primary' / f'{pmc}.xml'
        local = CORPUS / f'{pmc}.xml'
        member = f'scratch_K1b/fulltext/{pmc}.xml'
        if local.exists():
            shutil.copyfile(local, dst)
            source = str(local)
        elif member in members:
            dst.write_bytes(subprocess.check_output(['tar', '--zstd', '-xOf', str(ARCHIVE), member]))
            source = f'{ARCHIVE}::{member}'
        else:
            manifest.append({'pmcid': pmc, 'status': 'MISSING'})
            continue
        manifest.append({'path': str(dst.relative_to(ROOT)), 'source': source, 'sha256': sha(dst), 'bytes': dst.stat().st_size})
    for (name, member) in [('thr_harm.jsonl', 'scratch_K1b/kt/thr_harm.jsonl'), ('g_med_post.json', 'scratch_K1b/asm/g_med_post.json')]:
        dst = out / name
        dst.write_bytes(subprocess.check_output(['tar', '--zstd', '-xOf', str(ARCHIVE), member]))
        manifest.append({'path': str(dst.relative_to(ROOT)), 'source': f'{ARCHIVE}::{member}', 'sha256': sha(dst), 'bytes': dst.stat().st_size})
    result = {'files': manifest, 'wall_s': time.perf_counter() - t, 'total_bytes': sum((r.get('bytes', 0) for r in manifest)), 'cap_bytes': 3000000000}
    (out / 'SOURCE_MANIFEST.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for (k, v) in result.items() if k != 'files'}))
if __name__ == '__main__':
    main()
