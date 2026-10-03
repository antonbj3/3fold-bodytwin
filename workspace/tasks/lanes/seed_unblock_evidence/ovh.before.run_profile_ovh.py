"""Explicit account; a separate ephemeral OpenCode database for each cloud job."""
import argparse, collections, fcntl, hashlib, json, os, shutil, subprocess, time
from pathlib import Path
from route_runtime import choose_free_models,apply_route,run_guarded,reserve_slot
ROOT=Path('/opt/agents')
MODELS={'reserve_worker':'opencode-go/reserve_worker-v4.1-flash','swarm':'opencode/space-swarm-free','free_worker':'opencode/free_worker-2.5-preview-free'}
def select(label, requested, directory):
    return reserve_slot(label,requested,directory)
def environment(label,model,directory,web=False):
    provider,model_id=model.split('/',1)
    original=ROOT/'profiles'/label/'data/opencode/auth.json'
    auth=json.loads(original.read_text())
    if provider not in auth:raise SystemExit('Requested provider is not configured.')
    tag=hashlib.sha256((label+str(directory.resolve())).encode()).hexdigest()[:24]
    runtime=ROOT/'runtime'/tag;runtime.mkdir(parents=True,exist_ok=True,mode=0o700)
    env=os.environ.copy()
    for key in list(env):
        if key.endswith('_API_KEY') or key.startswith('OPENCODE_'):env.pop(key)
    for variable,folder in [('XDG_CONFIG_HOME','config'),('XDG_DATA_HOME','data'),('XDG_STATE_HOME','state'),('XDG_CACHE_HOME','cache')]:
        path=runtime/folder;path.mkdir(exist_ok=True,mode=0o700);env[variable]=str(path)
    target=runtime/'data/opencode/auth.json';target.parent.mkdir(exist_ok=True,mode=0o700)
    with os.fdopen(os.open(target,os.O_WRONLY|os.O_CREAT|os.O_TRUNC,0o600),'w') as f:json.dump({provider:auth[provider]},f)
    config={'model':model,'small_model':model,'default_agent':'build',
       'enabled_providers':[provider],'provider':{provider:{'whitelist':[model_id]}},
       'autoupdate':False,'share':'disabled',
       'permission':{'task':'deny','external_directory':'deny','webfetch':'allow' if web else 'deny','websearch':'allow' if web else 'deny'}}
    route=ROOT/'HOME_EGRESS.json'
    if route.is_file() and provider=='opencode':
        apply_route(env,json.loads(route.read_text()),model_id,directory)
    env.update(OPENCODE_DISABLE_PROJECT_CONFIG='true',OPENCODE_CONFIG_CONTENT=json.dumps(config),
       OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONPATH='/opt/bt/vendor')
    return env,runtime
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('profile',choices=('A','B','C','D'));parser.add_argument('model',choices=MODELS)
    parser.add_argument('--dir',required=True);parser.add_argument('--title',default='isolated_profile_job');parser.add_argument('prompt')
    args=parser.parse_args();directory=Path(args.dir).resolve();account,selected,marker=select(args.profile,args.model,directory);model=MODELS[selected]
    (directory/'MODEL_ROUTE.json').write_text(json.dumps({'account':account,'requested':args.model,'selected':selected,'model':model,'started':time.time()}))
    env,runtime=environment(account,model,directory,(directory/'ALLOW_WEB').exists())
    executable=str(ROOT/'.opencode/bin/opencode')
    try:
        completed=run_guarded([executable,'run','--pure','--auto','--agent','build','--model',model,
            '--dir',str(directory),'--format','json','--title',args.title,args.prompt],cwd=directory,env=env,runtime=runtime)
        code=completed.returncode
    finally:
        # Only this launcher's own ephemeral account store; never the configured profiles.
        shutil.rmtree(runtime)
        if marker is not None:marker.unlink(missing_ok=True)
    raise SystemExit(code)
