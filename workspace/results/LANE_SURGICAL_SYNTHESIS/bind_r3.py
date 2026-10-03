"""Lane-local review dispatch+feedback only; no source graph mutation."""
import os
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import hashlib,json,shutil,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
source=ROOT/'r3/ROUND_SUMMARY.json';digest=hashlib.sha256(source.read_bytes()).hexdigest()
dest=ROOT/'graph_review_workspace/bodytwin/results/ROUND_SUMMARY_R3.json'
if dest.exists():raise FileExistsError(dest)
shutil.copy2(source,dest)
receipts=[]
for suffix in ['INCISION','COLLAGEN','HEMOSTASIS','HEALING']:
    target='BT-CTX-SURG-'+suffix;lane='LANE_SURGICAL_SYNTHESIS_R3_'+suffix
    record=dict(target_id=target,result_file='results/ROUND_SUMMARY_R3.json',sha256=digest,review_state='PENDING_INDEPENDENT_REVIEW',outcome='SHARED_REFERENCE_PASS; AUTONOMOUS8.306_FAIL21.811; CHEMISTRY_FAIL44.340; CONTROL_TIE; NATIVE_UNKNOWN',measured_quantity='Collagen/tensor reference identity, integrated strength RMSE, autonomous precursor capacity and checkpoint equality',units='normalized fixed-reference collagen; mol sites/mol collagen reference; pp; CPU/wall seconds',uncertainty='Synthetic birth/removal and chemistry closures; native assay/traction/posterior and laboratory cost UNKNOWN',population_regime='Synthetic spatial skin transfer; separate mouse HP and rat strength diagnostics, inherited early strength fit',preregistered_gate='PREREG_R3; autonomous8.306+0.1pp failed;90/90 regression tests pass; conservative ledgers and time/space refinements pass',baseline='Same-information conventional FV/six directional cohorts matches exactly; no unique computation-cost gain',negative_result=True,original_lane_result='results/LANE_SURGICAL_SYNTHESIS/r3/ROUND_SUMMARY.json',source_context_import='PENDING_COORDINATOR; lane-local review receipts',scope='Incision/hemostasis closures retained; new explicit collagen reference and downstream inventory law, no empirical admission')
    rp=ROOT/f'GRAPH_FEEDBACK_R3_{suffix}.json'
    with rp.open('x') as f:json.dump(record,f,indent=2)
    dispatch=json.loads(subprocess.check_output([str(ROOT/'graph'),'dispatch','--id',target,'--lane',lane,'--kind','review','--reason','R3 shared collagen reference; preserve autonomous RMSE failure and chemistry capacity obstruction','--task-file',str(ROOT/'PREREG_R3.json')],text=True))
    feedback=json.loads(subprocess.check_output([str(ROOT/'graph'),'feedback','--lane',lane,'--record',str(rp)],text=True))
    receipts.append(dict(target=target,lane=lane,dispatch=dispatch,feedback=feedback))
    record['result_file']=record['original_lane_result']
    with (ROOT/f'GRAPH_COORDINATOR_FEEDBACK_R3_{suffix}.json').open('x') as f:json.dump(record,f,indent=2)
with (ROOT/'GRAPH_BINDING_STATUS_R3.json').open('x') as f:json.dump(dict(review_state='PENDING_INDEPENDENT_REVIEW',scientific_admission=False,local_receipts=receipts,R3_global_import='PENDING_COORDINATOR',write_boundary='lane-local review workspace; no canonical/source graph writes'),f,indent=2)
print('Four lane-local review dispatch+feedback pairs recorded; source graph unchanged')
