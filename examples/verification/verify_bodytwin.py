"""Explicit verification recipes; unmapped fresh declarations refuse full verification."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory(root):
    rows={}
    for line in (root/'docs/RUNNING.md').read_text().splitlines():
        cells=[c.strip().strip('`') for c in line.split('|')]
        if len(cells)>3 and cells[2]=='VERIFIED-FRESH':
            if cells[1] in rows:raise ValueError('Duplicate fresh declaration')
            rows[cells[1]]=hashlib.sha256(line.encode()).hexdigest()
    if not rows:raise ValueError('No fresh declarations found')
    return rows


def local_file(root,name):
    path=(root/name).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():raise ValueError('Recipe path is not a repository file')
    return path


def validate(root,recipes,rows):
    if not recipes:raise ValueError('No explicit recipes supplied')
    covered=set();identities=set()
    for recipe in recipes:
        if recipe['id'] in identities:raise ValueError('Duplicate recipe identity')
        identities.add(recipe['id'])
        budget=recipe.get('timeout_seconds',60)
        if type(budget) is not int or not 1<=budget<=900:raise ValueError('Invalid explicit recipe resource limit')
        if not recipe['members']:raise ValueError('Empty recipe membership')
        for name,expected in recipe['members'].items():
            if name in covered or rows.get(name)!=expected:raise ValueError('Missing, duplicated or changed status binding')
            covered.add(name)
            if name not in recipe['sources']:raise ValueError('Member source is not pinned')
        for group in ('sources','inputs'):
            for name,expected in recipe[group].items():
                if digest(local_file(root,name))!=expected:raise ValueError('Source or input binding mismatch')
        if digest(local_file(root,recipe['expected_report']))!=recipe['expected_report_sha256']:
            raise ValueError('Expected numerical report binding mismatch')
    return covered


def controls(root,recipes,rows):
    refused=0
    for kind in ('sources','inputs','members','expected','duplicate','unknown','empty','timeout'):
        changed=copy.deepcopy(recipes)
        recipe=changed[0]
        if kind in ('sources','inputs','members'):
            recipe[kind][next(iter(recipe[kind]))]='0'*64
        elif kind=='expected':recipe['expected_report_sha256']='0'*64
        elif kind=='duplicate':changed.append(copy.deepcopy(recipe))
        elif kind=='unknown':recipe['members']['missing.py']='0'*64
        elif kind=='timeout':recipe['timeout_seconds']=True
        else:changed=[]
        try:validate(root,changed,rows)
        except ValueError:refused+=1
    return refused==8


def run(root,recipes,artifacts=None):
    environment=os.environ.copy()
    environment.update(PYTHONPATH=str(root/'src'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    results=[]
    with tempfile.TemporaryDirectory(prefix='bodytwin-verify-') as directory:
        for recipe in recipes:
            hashes=[];passed=[]
            for leg in (0,1):
                output=Path(directory)/f"{recipe['id']}-{leg}.json"
                process=subprocess.run([sys.executable,*recipe['argv'],'--output',str(output)],cwd=root,
                                       env=environment,capture_output=True,text=True,timeout=recipe.get('timeout_seconds',60))
                if process.returncode!=0 or not output.is_file():raise ValueError('Numerical recipe failed or emitted no report')
                report=json.loads(output.read_text())
                if artifacts:
                    for artifact in (output,output.with_suffix('.npz')):
                        if artifact.exists():(artifacts/artifact.name).write_bytes(artifact.read_bytes())
                passed.append(bool(report.get('gates')) and all(v is True for v in report['gates'].values()))
                hashes.append(digest(output))
            results.append(dict(id=recipe['id'],hashes=hashes,passed=all(passed) and hashes[0]==hashes[1]==recipe['expected_report_sha256']))
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selected',action='store_true')
    parser.add_argument('--family',choices=('all','detectors','cells'),default='all')
    parser.add_argument('--inventory',action='store_true')
    parser.add_argument('--recipe',help='Execute only this explicitly declared recipe identity')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--artifacts-dir',type=Path)
    args=parser.parse_args()
    if args.output:
        args.output.write_text(json.dumps(dict(gates=dict(completed=False),scope='Verification has not completed.'))+'\n')
    if args.artifacts_dir:
        args.artifacts_dir.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[2]
    rows=inventory(root)
    recipes=json.loads((root/'docs/verification_recipes.json').read_text())['recipes']
    covered=validate(root,recipes,rows)
    missing=sorted(set(rows)-covered)
    report=dict(fresh_declarations=len(rows),covered=len(covered),unmapped=missing,
                scope='Explicit CPU recipes only; full verification refuses unmapped fresh declarations.')
    if args.inventory:
        report['gates']=dict(bindings=True)
    elif missing and not args.selected:
        report['gates']=dict(complete_coverage=False)
        report['numerical_children']=0
    else:
        selected=[r for r in recipes if (args.family=='all' or r['family']==args.family) and (args.recipe is None or r['id']==args.recipe)]
        if not selected:raise ValueError('Selected recipe family is empty')
        report['recipes']=run(root,selected,args.artifacts_dir)
        report['executed_family']=args.family
        report['gates']=dict(recipe_outputs=all(r['passed'] for r in report['recipes']),
                             refusal_controls=controls(root,recipes,rows))
    if args.output:args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({**report,'unmapped_count':len(missing),'unmapped':missing if args.inventory else 'See explicit inventory/output report.'}))
    return 0 if all(report['gates'].values()) else 2


if __name__=='__main__':raise SystemExit(main())
