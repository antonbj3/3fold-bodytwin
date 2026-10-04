#!/usr/bin/env python3
"""Fail a release containing private references, data payloads or unlisted files.

The trusted RELEASE_FILES.json is a review anchor, not a defense against an
attacker who replaces both the release and its anchor.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess

# Encoded patterns keep the scanner's own source free of restricted vocabulary.
PATTERNS = [bytes.fromhex(s).decode() for s in (
    '686f6c6c6f77', '6d696e647477696e', '626f64797477696e2d706572736f6e6c696774',
    '50454d46', '70656e696c', '6f73617964', '436f2d417574686f7265642d4279', '436c61756465',
    '43686174475054', '41492d67656e657261746564')]
BAD = re.compile('|'.join(map(re.escape, PATTERNS)), re.I)
PRIVATE_PATH = re.compile(r'(?<![\w:/@])(?:/home/[A-Za-z0-9._-]+|/mnt/(?:games-240|shared_data|INTENSO)|/media/[A-Za-z0-9._-]+)')
EMAIL = re.compile(r'[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}')
PERSON_FIELDS = re.compile(r'\b(?:patient_name|patient_id|date_of_birth|social_security_number|personal_identity_number)\b', re.I)
ALLOWED = {'.py','.sh','.md','.json','.toml','.txt','.svg','.c','.cpp','.h','.hpp'}
ALLOWED_NAMES = {'.gitignore', 'LICENSE', 'NOTICE'}
# The designated public address is allowed in Git identity fields only.
PUBLIC_IDENTITY_EMAIL = bytes.fromhex('616e746f6e626a334075736572732e6e6f7265706c792e6769746875622e636f6d').decode()
EXCLUDED = {'.git','__pycache__','.venv'}
MAX_BYTES = 10_000_000

def content_findings(name, raw, commit=False):
    issues=[]
    if len(raw)>MAX_BYTES:return ['file_over_10_MB']
    try: value=raw.decode('utf-8')
    except UnicodeError:return issues+['binary_payload']
    if '\x00' in value:issues.append('binary_payload')
    if BAD.search(name) or BAD.search(value):issues.append('restricted_vocabulary')
    if PRIVATE_PATH.search(value):issues.append('private_absolute_path')
    if '@' in value and EMAIL.search(value):issues.append('embedded_email')
    if PERSON_FIELDS.search(value) and Path(name).suffix not in {'.py','.c','.cpp','.h','.hpp'}:
        issues.append('personal_metadata')
    return issues

def python_literal_findings(raw):
    issues=[]
    try:
        tree=ast.parse(raw.decode())
        for node in ast.walk(tree):
            if isinstance(node,ast.Dict):
                for key,value in zip(node.keys,node.values):
                    if (isinstance(key,ast.Constant) and isinstance(key.value,str) and PERSON_FIELDS.fullmatch(key.value)
                        and isinstance(value,ast.Constant) and value.value not in (None,'','synthetic','synthetic-wrong-identity','example','test','UNKNOWN')):
                        issues.append('personal_metadata_literal')
            if isinstance(node,(ast.List,ast.Tuple)) and sum(isinstance(x,ast.Constant) and isinstance(x.value,(int,float)) for x in ast.walk(node))>128:
                issues.append('large_embedded_numeric_payload');break
    except (SyntaxError,UnicodeError):issues.append('invalid_python')
    return issues

def scan(root, verify_manifest=True, inspect_git=True):
    root=Path(root);findings=[];files={};bytes_total=0
    for p in sorted(root.rglob('*')):
        rel=p.relative_to(root)
        if any(part in EXCLUDED for part in rel.parts):continue
        name=rel.as_posix()
        if p.is_symlink():
            findings.append({'file':name,'reason':'symlink'});continue
        if not p.is_file():continue
        size=p.stat().st_size;bytes_total+=size
        if size>MAX_BYTES:
            files[name]='OVERSIZED'
            findings.append({'file':name,'reason':'file_over_10_MB'})
            continue
        raw=p.read_bytes();files[name]=hashlib.sha256(raw).hexdigest()
        if p.suffix not in ALLOWED and p.name not in ALLOWED_NAMES:
            findings.append({'file':name,'reason':'unapproved_payload_extension'})
        for reason in content_findings(name,raw):findings.append({'file':name,'reason':reason})
        if p.suffix=='.py':
            for reason in python_literal_findings(raw):findings.append({'file':name,'reason':reason})
    if verify_manifest:
        anchor=root/'RELEASE_FILES.json'
        if not anchor.exists():findings.append({'file':'RELEASE_FILES.json','reason':'missing_release_anchor'})
        else:
            expected=json.loads(anchor.read_text())['files']
            actual={k:v for k,v in files.items() if k!='RELEASE_FILES.json'}
            for name in sorted(set(expected)|set(actual)):
                if expected.get(name)!=actual.get(name):
                    findings.append({'file':name,'reason':'unlisted_or_changed_release_file'})
    if inspect_git and (root/'.git').exists():
        q=subprocess.run(['git','-C',str(root),'rev-list','--all'],capture_output=True,text=True,check=True)
        for commit in q.stdout.splitlines():
            identity=subprocess.run(['git','-C',str(root),'show','-s','--format=%an%n%ae%n%cn%n%ce',commit],capture_output=True,text=True,check=True).stdout.splitlines()
            for index,value in enumerate(identity):
                if index in (1,3) and value == PUBLIC_IDENTITY_EMAIL:
                    continue
                for reason in content_findings('commit identity',value.encode(),True):
                    findings.append({'file':'commit '+commit,'reason':reason,'identity_field':('author_name','author_email','committer_name','committer_email')[index]})
            message=subprocess.run(['git','-C',str(root),'show','-s','--format=%B',commit],capture_output=True,check=True).stdout
            for reason in content_findings('commit '+commit,message,True):findings.append({'file':'commit '+commit,'reason':reason})
            listing=subprocess.run(['git','-C',str(root),'ls-tree','-rz',commit],capture_output=True,check=True).stdout
            for entry in listing.split(b'\x00'):
                if not entry:continue
                meta,rawname=entry.split(b'\t',1);obj=meta.decode().split()[2];name=rawname.decode()
                if meta.startswith(b'120000 '):
                    findings.append({'file':name,'commit':commit,'reason':'symlink'})
                size=int(subprocess.run(['git','-C',str(root),'cat-file','-s',obj],capture_output=True,check=True).stdout)
                if size>MAX_BYTES:
                    findings.append({'file':name,'commit':commit,'reason':'file_over_10_MB'})
                    continue
                raw=subprocess.run(['git','-C',str(root),'cat-file','blob',obj],capture_output=True,check=True).stdout
                for reason in content_findings(name,raw):findings.append({'file':name,'commit':commit,'reason':reason})
                if Path(name).suffix=='.py':
                    for reason in python_literal_findings(raw):findings.append({'file':name,'commit':commit,'reason':reason})
                if Path(name).suffix not in ALLOWED and Path(name).name not in ALLOWED_NAMES:
                    findings.append({'file':name,'commit':commit,'reason':'unapproved_payload_extension'})
    return {'status':'FAIL' if findings else 'PASS','file_count':len(files),'bytes':bytes_total,'findings':findings,
            'scope':'working release inventory, text contents, Python literals and all reachable commit identities/messages/blobs; designated public no-reply identity allowed only in Git email fields'}

def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);a.add_argument('--no-manifest',action='store_true');a.add_argument('--output',type=Path)
    args=a.parse_args();r=scan(args.root,not args.no_manifest)
    payload=json.dumps(r,indent=2)+'\n'
    if args.output:args.output.write_text(payload)
    print(payload,end='');return 0 if r['status']=='PASS' else 1

if __name__=='__main__':raise SystemExit(main())
