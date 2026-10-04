"""Import every retained Python source file with isolated legacy module names.

Source-only drivers may declare missing external inputs. Missing Python modules
are failures, including when a predecessor was omitted from the catalogue.
Imports do not run a data-dependent pipeline or establish physical validity.
"""
import argparse
import ast
import contextlib
import importlib.util
import importlib
from importlib.machinery import PathFinder, BuiltinImporter, FrozenImporter
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]

class ImportSideEffect(RuntimeError):
    pass

def files(entry):
    return sorted(p for p in (ROOT/'implementations'/entry).rglob('*.py') if '__pycache__' not in p.parts)

def import_name(path):
    package=[]
    parent=path.parent
    while (parent/'__init__.py').is_file():
        package.insert(0,parent.name);parent=parent.parent
    name='.'.join(package+[path.stem]) if package else 'release_import_probe_'+path.stem
    if path.name=='__init__.py':name='.'.join(package)
    return name,parent

def declared_missing(module, entry, contracts):
    top=module.split('.')[0]
    external=contracts['external_code'].get(entry, [])
    if any(module==name or module.startswith(name+'.') for name in external):return True
    return top in contracts['optional_packages'] and importlib.util.find_spec(top) is None

def module_exists(name):
    parts=name.split('.')
    search=sys.path
    full=''
    for part in parts:
        full=part if not full else full+'.'+part
        # Resolve a child from the parent's filesystem locations without
        # creating a namespace spec that expects an executed parent module.
        spec=(BuiltinImporter.find_spec(full) or FrozenImporter.find_spec(full) or PathFinder.find_spec(part,search))
        if spec is None:return False
        search=spec.submodule_search_locations
        if search is None and full!=name:return False
    return True

def preflight(path, entry, contracts):
    """Check top-level imports even after a missing-input initializer.

    Function-local runtime imports are exercised only by the full external
    pipeline. Do not execute package initializers while resolving source files.
    """
    tree=ast.parse(path.read_text())
    names=[]
    def visit(node):
        if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef)):return
        if isinstance(node,ast.If) and ast.unparse(node.test) in {"__name__ == '__main__'", "'__main__' == __name__"}:
            for child in node.orelse:visit(child)
            return
        if isinstance(node,ast.Import):names.extend(alias.name for alias in node.names)
        if isinstance(node,ast.ImportFrom) and node.module and node.level==0:names.append(node.module)
        for child in ast.iter_child_nodes(node):visit(child)
    visit(tree)
    return sorted({name for name in names if not module_exists(name) and not declared_missing(name,entry,contracts)})

def worker(entry):
    sys.path.insert(0,str(ROOT))
    from dental_release.paths import MissingInput
    baseline=list(sys.path)
    report=[]
    contracts=json.loads((ROOT/'provenance/IMPORT_CONTRACTS.json').read_text())
    extra=[str(ROOT/'implementations'/p) for p in contracts['paths'].get(entry, [])]
    # Plotting library initialization is allowed once in a temporary config
    # directory; guards below concern producer driver side effects.
    try:
        import matplotlib
        import matplotlib.pyplot
    except ModuleNotFoundError:
        pass
    # Driver imports must never launch tools, write files, or start acquisition.
    def guard(event,args):
        if event in {'subprocess.Popen','os.system','os.mkdir','os.rename','os.remove','os.rmdir','os.symlink','socket.connect','socket.bind'}:
            raise ImportSideEffect('Source driver performs '+event+' during import')
        if event=='open':
            mode,flags=args[1:3]
            if isinstance(mode,str) and any(c in mode for c in 'wax+') or isinstance(flags,int) and flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND):
                raise ImportSideEffect('Source driver writes during import')
    sys.addaudithook(guard)
    def timeout(signum,frame):raise ImportSideEffect('Source driver computes during import beyond 3 seconds')
    signal.signal(signal.SIGALRM,timeout)
    for p in files(entry):
        # Clear only source modules; installed dependencies stay warm, preserving
        # their normal singleton state. No other implementation is used as a fallback.
        for key,module in list(sys.modules.items()):
            location=getattr(module,'__file__','') or ''
            if str(ROOT/'implementations') in location:sys.modules.pop(key,None)
        name,parent=import_name(p)
        sys.path[:]=[str(p.parent),str(parent),str(ROOT/'implementations'/entry),str(ROOT/'implementations'/entry/'code')]+extra+baseline
        row=dict(file=str(p.relative_to(ROOT)),entry=entry)
        try:
            missing=preflight(p,entry,contracts)
            if missing:raise ModuleNotFoundError('Undeclared top-level import(s): '+', '.join(missing),name=missing[0])
            with contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()):
                signal.alarm(3)
                spec=importlib.util.spec_from_file_location(name,p,submodule_search_locations=[str(p.parent)] if p.name=='__init__.py' else None)
                if not name.startswith('release_import_probe_'):
                    # Normal package loading initializes the parent before the
                    # child. Pre-inserting a half-initialized child creates an
                    # artificial cycle in packages that export their public API.
                    importlib.import_module(name)
                else:
                    module=importlib.util.module_from_spec(spec)
                    sys.modules[name]=module
                    spec.loader.exec_module(module)
            row['status']='IMPORTED'
        except (MissingInput,FileNotFoundError,ImportSideEffect) as e:
            row.update(status='EXTERNAL_INPUT_OR_DRIVER',exception=type(e).__name__,detail=str(e))
        except BaseException as e:
            row.update(status='FAIL',exception=type(e).__name__,detail=str(e))
            if isinstance(e,ModuleNotFoundError):
                row['missing_module']=e.name
                top=(e.name or '').split('.')[0]
                if declared_missing(e.name,entry,contracts):
                    row['status']='EXTERNAL_PACKAGE_OR_CODE'
            if type(e).__name__=='Skipped' and '/tests/' in row['file'] and str(e)=='optional field contact port package':
                row['status']='EXTERNAL_PACKAGE_OR_CODE'
        finally:signal.alarm(0)
        report.append(row)
    return report

def check(entries=None):
    catalogue=json.loads((ROOT/'demos.json').read_text())['demos']
    entries=entries or [e['id'] for e in catalogue if 'document' not in e]
    env={k:v for k,v in os.environ.items() if k in ('PATH','SYSTEMROOT','LD_LIBRARY_PATH')}
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1')
    rows=[]
    for entry in entries:
        expected=[str(p.relative_to(ROOT)) for p in files(entry)]
        with tempfile.TemporaryDirectory(prefix='dental-import-') as temp:
            proc=subprocess.run([sys.executable,'-I','-B',str(Path(__file__).resolve()),'--worker',entry],cwd=temp,env=env,text=True,capture_output=True,timeout=120)
        if proc.returncode:
            rows.extend(dict(file=p,entry=entry,status='FAIL',exception='WorkerExit',detail=proc.stderr[-1500:]) for p in expected)
            continue
        payload=json.loads(proc.stdout)
        if [r['file'] for r in payload]!=expected:raise AssertionError('Import census lost files: '+entry)
        rows.extend(payload)
    errors=[r for r in rows if r['status']=='FAIL']
    return dict(status='FAIL' if errors else 'PASS',entries=len(entries),python_files=len(rows),imported=sum(r['status']=='IMPORTED' for r in rows),external_input_or_driver=sum(r['status']=='EXTERNAL_INPUT_OR_DRIVER' for r in rows),external_package_or_code=sum(r['status']=='EXTERNAL_PACKAGE_OR_CODE' for r in rows),failures=len(errors),rows=rows,scope='Every retained Python file is executed as an import. External-input/driver initialization and explicitly declared external packages or code are reported separately. Every undeclared missing Python module fails; no full source pipeline replay is implied.')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--worker');p.add_argument('--entry',action='append');p.add_argument('--output',type=Path)
    a=p.parse_args()
    result=worker(a.worker) if a.worker else check(a.entry)
    text=json.dumps(result,indent=2)+'\n'
    if a.output:a.output.write_text(text)
    print(text,end='')
    return 0 if a.worker or result['status']=='PASS' else 1

if __name__=='__main__':raise SystemExit(main())
