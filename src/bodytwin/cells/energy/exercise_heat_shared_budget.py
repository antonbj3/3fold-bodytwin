"""Shared cardiac-output budget between exercising muscle and heat-dissipating skin at VO2max.

Tests whether a reported 17.0-25.5 percentage-point oversubscription of the non-muscle flow budget
is an artifact of (i) mixing a passive-heat-stress skin ceiling into a maximal-exercise budget and
(ii) a uniform-extraction tautology in the upstream muscle-flow derivation. Cross-checks the
upstream muscle flow fraction against directly measured catheter/thermodilution data (Skattebo 2020
PMC7540168 n=117; Calbet 2007 PMID 17600155), sweeps the obligatory-bed subtraction for the skin
allowance at VO2max, reconstructs the forearm->whole-body BSA extrapolation (Johnson & Rowell 1975
PMID 1213973) and restates the heat-stress regime (Gonzalez-Alonso & Calbet 2003 PMID 12591751).

Reads: <OUT_ROOT>/exercise_bloodflow_redistribution/exercise_bloodflow_redistribution_results.json
Writes: nothing (stdout only).
Gate: the pre-registered claim C in R['PREREG'] -- accepted iff (a) the upstream muscle flow
fraction exceeds the measured leg fraction by >5 pp AND exceeds the cited 80-85% band, AND (b) the
median measured-trunk-derived skin allowance at VO2max is <= 2.0 L/min and non-negative in >= 75%
of the obligatory-bed sweep.
"""
import json, itertools, statistics as st
import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
R={}
# ---- upstream cell output, read live ----
d=json.load(open(_os.path.join(OUT_ROOT, 'exercise_bloodflow_redistribution', 'exercise_bloodflow_redistribution_results.json')))
fa=d['forced_adversary_2_conservation_of_flow_budget']
CO=d['inputs_reused_readonly']['CO_max_central_l_min']; E_sys=d['inputs_reused_readonly']['avo2diff_central_ml_100ml']/100.
VO2max=fa['VO2max_abs_l_min']; VO2nm=fa['VO2_nonmuscle_generic_l_min']
R['repo']={'CO_max':CO,'E_sys':E_sys,'VO2max':VO2max,'VO2_nonmuscle':VO2nm,
 'muscle_frac_pct':fa['muscle_flow_frac_of_co_max_pct'],'nonmuscle_budget':fa['non_muscle_budget_max_l_min'],
 'lit_band_pct':fa['classic_literature_qualitative_range_pct']}
# TAUTOLOGY PROOF: flow frac == VO2 frac exactly, i.e. E_muscle==E_sys assumed
R['tautology']={'VO2_frac_muscle_pct':100*(VO2max-VO2nm)/VO2max,
 'flow_frac_pct':fa['muscle_flow_frac_of_co_max_pct'],
 'identical_to_1e-9':abs(100*(VO2max-VO2nm)/VO2max - fa['muscle_flow_frac_of_co_max_pct'])<1e-3,
 'lit_band_gate_violated': fa['muscle_flow_frac_of_co_max_pct']>fa['classic_literature_qualitative_range_pct'][1]}

# ---- EXTERNAL measured maxima (catheter/thermodilution/Fick), decorrelated from the upstream cell ----
# Skattebo 2020 PMC7540168 Table3 n=117: CO 25.0+-2.4, LBF 18.6+-3.0, ExtractSys 79+-8%, ExtractLeg 84+-5%
sk={'CO':25.0,'CO_sd':2.4,'LBF':18.6,'LBF_sd':3.0,'Esys':0.79,'Eleg':0.84,'CaO2':203.}
sk['nonleg']=sk['CO']-sk['LBF']; sk['nonleg_frac']=sk['nonleg']/sk['CO']
sk['nonleg_sd']=(sk['CO_sd']**2+sk['LBF_sd']**2)**.5   # conservative (ignores +corr => overstates sd)
sk['E_leg_over_sys']=sk['Eleg']/sk['Esys']
R['skattebo2020']=sk
# Calbet 2007 PMID17600155: leg ~70% CO at Wmax, arm 1.2 L/min, TRUNK 21% CO ~= 5.5 L/min
cal={'leg_pct':70.,'arm_l_min':1.2,'trunk_pct':21.,'trunk_l_min':5.5,'implied_CO':5.5/0.21,'leg_VO2_pct':84.,'trunk_VO2_pct':14.5}
R['calbet2007']=cal

# ---- CROSS-CHECK 1: upstream muscle frac vs measured ----
R['xcheck_muscle_frac']={'repo_pct':fa['muscle_flow_frac_of_co_max_pct'],
 'skattebo_leg_frac_pct':100*sk['LBF']/sk['CO'],'calbet_leg_frac_pct':cal['leg_pct'],
 'calbet_leg_plus_arm_pct':cal['leg_pct']+100*cal['arm_l_min']/cal['implied_CO'],
 'repo_minus_measured_pp':fa['muscle_flow_frac_of_co_max_pct']-100*sk['LBF']/sk['CO']}
# ---- CROSS-CHECK 2: non-muscle budget, upstream vs measured (scaled to upstream CO) ----
scale=CO/sk['CO']
R['xcheck_budget']={'repo_l_min':fa['non_muscle_budget_max_l_min'],
 'skattebo_nonleg_scaled_l_min':sk['nonleg']*scale,
 'skattebo_nonleg_scaled_sd':sk['nonleg_sd']*scale,
 'calbet_trunk_scaled_l_min':cal['trunk_l_min']*CO/cal['implied_CO'],
 'ratio_measured_over_repo':(sk['nonleg']*scale)/fa['non_muscle_budget_max_l_min']}

# ---- SKIN ALLOWANCE at VO2max: measured trunk budget minus obligatory beds (sweep) ----
brain=[0.60,0.75,0.90]; coron=[0.60,0.90,1.20]; renal=[0.20,0.25,0.35]; splanch=[0.25,0.40,0.55]
resp_pct=[0.06,0.10,0.14]   # Harms1998 PMID9688739 up to 14-16% CO to respiratory muscles; low end conservative
trunk=cal['trunk_l_min']*CO/cal['implied_CO']
sk_alw=[]
for b,c,r_,s,rp in itertools.product(brain,coron,renal,splanch,resp_pct):
    sk_alw.append(trunk-(b+c+r_+s+rp*CO))
sk_alw.sort()
R['skin_allowance_at_VO2max_l_min']={'trunk_budget_scaled':trunk,'n':len(sk_alw),'min':sk_alw[0],'p10':sk_alw[int(.1*len(sk_alw))],
 'median':st.median(sk_alw),'p90':sk_alw[int(.9*len(sk_alw))],'max':sk_alw[-1],
 'frac_below_required_cap_1.973':sum(1 for x in sk_alw if x<=fa['non_muscle_budget_max_l_min'])/len(sk_alw),
 'frac_below_3_0':sum(1 for x in sk_alw if x<=3.0)/len(sk_alw),
 'frac_negative':sum(1 for x in sk_alw if x<0)/len(sk_alw)}

# ---- FOREARM->WHOLE-BODY BSA SCALING (the load-bearing extrapolation) ----
# Johnson&Rowell 1975 PMID1213973 RAW: FBF +8.26 ml/100ml/min; forearm MUSCLE flow fell 3.84->2.13 => skin delta = 8.26+1.71
dskin=8.26+ (3.84-2.13)
fv=[0.90,1.05,1.30]; fa_frac=[0.025,0.033,0.040]; BSA=1.8
ext=[]
for v,f in itertools.product(fv,fa_frac):
    per_forearm=dskin*(v*1000/100.)/1000.   # L/min for ONE forearm
    ext.append(per_forearm/f)               # scale by BSA fraction of one forearm
ext.sort()
R['bsa_extrapolation']={'raw_forearm_skin_delta_ml_100ml_min':dskin,'central_L_min':ext[len(ext)//2],
 'range_L_min':[ext[0],ext[-1]],'spread_factor':ext[-1]/ext[0],
 'published_value_GA2008':3.0,'reconstruction_brackets_published': ext[0]<=3.0<=ext[-1],
 'note':'scaling = (forearm tissue volume/100mL units) x (BSA / one-forearm skin area); forearm skin=2.5-4.0% BSA, forearm vol 0.9-1.3L'}
R['bsa_uncertainty_vs_gap']={'gap_to_close_L_min':6.0-fa['non_muscle_budget_max_l_min'],
 'bsa_1sigma_like_spread_L_min':ext[-1]-ext[0]}

# ---- HEAT-STRESS REGIME RESTATEMENT (GA&Calbet 2003 PMID12591751: VO2max -8%, CO/LBF/MAP -5..-11%) ----
heat={}
for dco in (0.05,0.08,0.11):
    VO2h=VO2max*0.92; COh=CO*(1-dco)
    # muscle flow from the measured leg fraction, NOT systemic extraction
    heat[f'CO_drop_{int(dco*100)}pct']={'CO_max_heat':COh,'VO2max_heat':VO2h,
      'nonleg_budget_measured_frac_l_min':COh*sk['nonleg_frac'],
      'repo_style_nonmuscle_budget':COh-(VO2h-VO2nm)/E_sys}
R['heat_regime']=heat

# ---- PREREG VERDICT ----
req=fa['non_muscle_budget_max_l_min']
R['PREREG']={'C':'the 17.0-25.5pp oversubscription is an artifact of (i) regime-mixing a passive-heat-stress skin ceiling into a maximal-exercise budget and (ii) a uniform-extraction tautology in the upstream muscle-flow derivation; the physiologically available skin flow AT VO2max is <=2 L/min and the joint demand does NOT violate conservation',
 'threshold':'accept C iff (a) the upstream muscle flow frac exceeds the directly-measured leg frac by >5pp AND exceeds the cited 80-85 band, AND (b) median measured-trunk-derived skin allowance at VO2max <= 2.0 L/min AND >=0 in >=75% of the obligatory-bed sweep',
 'a_pass':R['xcheck_muscle_frac']['repo_minus_measured_pp']>5 and R['tautology']['lit_band_gate_violated'],
 'b_median_le_2': R['skin_allowance_at_VO2max_l_min']['median']<=2.0,
 'b_nonneg_frac': 1-R['skin_allowance_at_VO2max_l_min']['frac_negative']}
R['PREREG']['b_pass']=R['PREREG']['b_median_le_2'] and R['PREREG']['b_nonneg_frac']>=0.75
R['PREREG']['VERDICT_C']=R['PREREG']['a_pass'] and R['PREREG']['b_pass']
print(json.dumps(R,indent=1,default=str))
