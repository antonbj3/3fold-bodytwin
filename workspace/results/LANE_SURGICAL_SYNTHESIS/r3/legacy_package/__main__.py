import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import argparse,json
from pathlib import Path
from .scenario import scenario_config,evidence_config

def main():
 p=argparse.ArgumentParser(description='Conditional surgical chain; default evidence mode keeps native outputs null')
 p.add_argument('--config',type=Path);p.add_argument('--mode',choices=['evidence','scenario'],default='evidence');p.add_argument('--out',type=Path,required=True);p.add_argument('--resume',type=Path)
 a=p.parse_args()
 from .chain import run,resume
 if a.resume:
  d=resume(a.resume,a.out);print(json.dumps({'restart_strength90_pct':float(d['strength_pct'][-1]),'review_state':'PENDING_INDEPENDENT_REVIEW'}));return
 cfg=json.loads(a.config.read_text()) if a.config else scenario_config() if a.mode=='scenario' else evidence_config()
 d=run(cfg,mode=a.mode,out=a.out)
 print(json.dumps({'mode':d['mode'],'empirical_joint_gate':d['empirical_joint_gate'],'cells':d.get('numerics',{}).get('cells'),'strength90_pct':d.get('strength_time',{}).get('value',{}).get('pct',[None])[-1]},ensure_ascii=False))
if __name__=='__main__':main()
