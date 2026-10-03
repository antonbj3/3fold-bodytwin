from pathlib import Path
import json
R=Path(__file__).resolve().parent
rows=[]
def add(step,port,value,unit,status,uncertainty,source,scope):
 rows.append(dict(step=step,port=port,value=value,unit=unit,status=status,uncertainty=uncertainty,source=source,scope=scope))
i=json.loads((R/'sources/LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json').read_text())
b=json.loads((R/'sources/LANE_SURGICAL_BINDINGS/PORTS_R4_FINAL.json').read_text())
g=json.loads((R/'sources/LANE_SKIN_TOUGHNESS_GAP/PORTS_R2_FINAL_V2.json').read_text())
r=json.loads((R/'sources/LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R3.json').read_text())
for a in i['Gamma_pierce_effective']:add('verktyg→mekanik',a['tool']+'.effective_J',a['effective_J_J_m2'],'J/m²',"MEASURED",a['reading_band_J_m2'],'LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json',a['mechanism_scope']+'; fullskin, tool-specific, figure reading not CI')
for name in ['Gamma_cut_empirical_J_m2','Gamma0_intrinsic_J_m2','wound_gap_empirical_m','bleeding_flow_m3_s','vessel_map','healing_strength_law']:
 add('snitt→respons',name,i[name],'m' if name.endswith('_m') else 'm³/s' if 'flow' in name else 'J/m²' if 'Gamma' in name else 'documented map/law',"UNKNOWN",None,'LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json','Cannot fill from another assay')
for name in ['mechanical_strain_zone_halfwidth_m','cell_viability_zone_halfwidth_m','perfusion_loss_zone_halfwidth_m']:
 add('skadezon',name,None,'m',"UNKNOWN",None,'LANE_SURGICAL_INCISION/PORTS_R4_FINAL.json','Distinct widths; no mechanical→biological transfer')
add('snittarbete','Gamma_cut_scenario',[150,380],'J/m²','SYNTETISKT','scenario envelope, not CI','LANE_SURGICAL_BINDINGS/PORTS_R4_FINAL.json','Explicit opt-in; not native cut')
add('mikro→makro','fibril_inventory_upper',6965.365853658535,'J/m²','SYNTETISKT','borrowed interface bound and imposed V/A','LANE_SURGICAL_BINDINGS/INTERFACE_BUDGET_RESULTS_R4.json','Conditional ceiling, not universal native upper bound')
add('bryggning','mode_work_reference',[30380,20600],'J/m²',"MEASURED",'SD 4900/2150; one juvenile porcine fullskin cohort','LANE_SURGICAL_BINDINGS/PORTS_R4_FINAL.json','Tear assay I/III, not cut or healing stress')
add('bryggning','native_terminal_traction',None,'Pa',"UNKNOWN",None,'LANE_SKIN_TOUGHNESS_GAP/PORTS_R2_FINAL_V2.json','Mode recruitment/anchor survival and cut connectivity absent')
add('bryggning','mm_bridge_work',[24500,19137.5],'J/m²','SYNTETISKT','L5mm,eps.5 strength-as-traction; not CI','LANE_SKIN_TOUGHNESS_GAP/MECHANISM_TABLE_R2_FINAL.json','Corrected crack/opening axes; anisotropy FAIL')
add("bleeding",'vessel_radius_histogram',None,'m and 1/m² per layer',"UNKNOWN",None,'BT-FW48-AUTO-55a5902a27c3f6/vessel_map_flow_results.json','Source histogram explicitly assumed; count control OVERpredicts25.3757x')
add('hemostas','shear_adhesion_law',None,'1/min versus 1/s',"UNKNOWN",None,'BT-FW48-AUTO-27a6dac22a7a23/platelet_hemostasis_shear.py','Qualitative anchors support synthetic rate closure; no native shear/closure assay')
for name,unit in [('dry_mass_density','kg/m³'),('vascular_flow_Hb_map','m³/s and mlO2/m³'),('oxygen_boundary','Torr and mlO2/(m² min Torr)'),('cell_damage_law','1/min')]:
 add('syrekant',name,None,unit,"UNKNOWN",None,'LANE_SURGICAL_RESPONSE/RESULTS.md','Same-specimen quantities not identified by resting human skin proxies')
for name,unit in [('physical_collagen_birth_turnover','reference mass/day'),('catalytic_LOX','1/day'),('joint_Q_U_Q_I_Q_M','reference mass tensor'),('chemistry_to_strength_law','Pa or %intact'),('chemistry_61_90','mol/molcollagen')]:
 add("healing→strength",name,None,unit,"UNKNOWN",None,'LANE_SURGICAL_RESPONSE/RESPONSE_PORTS_R3.json','Tracer≠mass, HP≠total maturity; species/cohort joint posterior null')
add("healing→strength",'R3_HP_input_strength_RMSE',8.306367094671637,'percentage points','SYNTETISKT','conditional observation proxy, nominal gate FAIL','LANE_SURGICAL_RESPONSE/r3/crosslink_port/summary.json','HP42 input, not predicted from early chemistry; Levenson heldout diagnostic')
with (R/'PORT_TABLE_R1.json').open('x') as f:json.dump({'review_state':'PENDING_INDEPENDENT_REVIEW','ports':rows},f,indent=2,ensure_ascii=False)
with (R/'PORT_TABLE_R1.md').open('x') as f:
 f.write("# Port table before chain run\n\nMEASURED indicates the source observation within the stated assay; no shared cohorts are assumed. SYNTETISKT indicates the chosen closure/computation. UNKNOWN remains null.\n\n| Step | Port | Value and unit | Uncertainty | Status | Source file / scope |\n|---|---|---|---|---|---|\n")
 for x in rows:f.write('| '+' | '.join(str(x[k]).replace('|','/') for k in ['step','port'])+' | '+str(x['value'])+' '+x['unit']+' | '+str(x['uncertainty'])+' | '+x['status']+' | '+x['source']+'; '+x['scope']+' |\n')
print(len(rows),'dimensioned ports')
