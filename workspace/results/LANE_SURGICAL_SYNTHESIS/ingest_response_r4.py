from pathlib import Path
import json,hashlib,shutil,datetime
R=Path(__file__).resolve().parent;P=R.parent/'LANE_SURGICAL_RESPONSE'
files=[]
for name in ['RESPONSE_PORTS_R4.json','VERIFICATION_R4.json','CHECKPOINT_R4.json','NEXT_ROUND_R4_FINAL.md','DERIVATION_R4.md','OBSERVATIONS_R4.json','PREREG_R4_GLYCO_RESERVOIR.json','SOURCE_REVIEW_R4.md','r4/SCORES_R4_v2.json','r4/glyco_reservoir/summary.json','night_rounds/r4.json']:
 d=R/'sources/LANE_SURGICAL_RESPONSE'/name;d.parent.mkdir(parents=True,exist_ok=True)
 if d.exists():raise FileExistsError(d)
 shutil.copyfile(P/name,d);files.append({'source':str(P/name),'snapshot':str(d.relative_to(R)),'sha256':hashlib.sha256(d.read_bytes()).hexdigest()})
for name in ['RESULTS.md','WORK_STATUS.json','NEXT_ROUND.md']:
 d=R/'sources/LANE_SURGICAL_RESPONSE'/(name+'.R4_FINAL');shutil.copyfile(P/name,d);files.append({'source':str(P/name),'snapshot':str(d.relative_to(R)),'sha256':hashlib.sha256(d.read_bytes()).hexdigest()})
with (R/'SOURCE_MANIFEST_R1_RESPONSE_R4.json').open('x') as f:json.dump({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':files},f,indent=2)
pre={'operation':'Consume final RESPONSE R4 as isolated D,H,A chemistry branch with explicit unknown chemistry→native U/I/M rupture map','nominal':'D10=.066,H10=.018,A10=0; independent control first-orderk=.06949679103768976/day,eta=.2833333333333334; pool is synthetic native transfer','prediction':'H42=.026338454833019345 in fixed C10 reference; do not reinterpret as current U/I/M maturity','control':'Same autonomous conservative reaction formula=TIE','gates':'R4 denominator-onlyFAIL, all16glycojointFAIL, nominal23.1798ppFAIL retained; no replacement by best profile','native_dependency':'chemistry_angle_bridge_joint/healing_traction remains null; chemistry branch cannot certify native strength'}
p=R/'PREREG_R1_RESPONSE_R4.json'
with p.open('x') as f:json.dump(pre,f,indent=2)
p.with_suffix('.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest()+'\n')
rows=json.loads((R/'PORT_TABLE_R1.json').read_text())['ports'];rows.extend([
 dict(step='slutlig kemi R4',port='native precursor D10',value=None,unit='mol/molreference collagen',status="UNKNOWN",uncertainty=None,source='LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R4.json',scope='Same reference and absolute glyco pools missing'),
 dict(step='slutlig kemi R4',port='cell-free first-order control k,eta',value=[.06949679103768976,.2833333333333334],unit='1/day,1',status="MEASURED",uncertainty='Apparent endpoint-derived analogue; first/second order nonidentified',source='LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R4.json',scope='External cell-culture typeI assay, native transfer synthetic'),
 dict(step='slutlig kemi R4',port='HP inventory42/10',value=2.220261437908497,unit='1',status="MEASURED",uncertainty=[1.747159090909091,2.8467987804878048],source='LANE_SURGICAL_RESPONSE/r4/SCORES_R4_v2.json',scope='Ratio of group means; fixed assay reference; graphical box, not CI')])
with (R/'PORT_TABLE_R1_FINAL.json').open('x') as f:json.dump({'review_state':'PENDING_INDEPENDENT_REVIEW','ports':rows},f,indent=2,ensure_ascii=False)
with (R/'PORT_TABLE_R1_FINAL.md').open('x') as f:
 f.write("# Final port table, including RESPONSE R4\n\n| Step | Port | Value / unit | Uncertainty | Status | Source / scope |\n|---|---|---|---|---|---|\n")
 for x in rows:f.write('| '+' | '.join([x['step'],x['port'],str(x['value'])+' '+x['unit'],str(x['uncertainty']),x['status'],x['source']+'; '+x['scope']]).replace('\n',' ')+' |\n')
print('R4 frozen last; final table',len(rows),'ports')
