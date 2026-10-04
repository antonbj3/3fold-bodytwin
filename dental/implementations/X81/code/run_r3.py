import copy, hashlib, json, math, re
import numpy as np
from scipy.stats import t
from lxml import etree
from common import ROOT, read, save, state, verify_freeze

def added(key):
    r = verify_freeze('INPUT_ADDITION_R3.json')[key]
    p = ROOT / r['snapshot']
    assert hashlib.sha256(p.read_bytes()).hexdigest() == r['sha256']
    return json.loads(p.read_text()) if p.suffix == '.json' else p

def force_intervals():
    curated = added('curated_force')
    res = []
    for (i, q) in enumerate(read('force')['external_force']['rows']):
        row = {'edge_id': f'D-E-X70-FORCE-{i + 1}', 'observed_N': q['observed_N'], 'predicted_N': q['predicted_N'], 'fixed_prediction_acceptance_N': [q['predicted_N'] * math.exp(-0.2), q['predicted_N'] * math.exp(0.2)], 'original_point_gate': q['gate'], 'locator': q['locator'], 'resolution': 'POPULATION'}
        matches = [z for z in curated if z.get('origin') == 'measured' and abs(z['mean'] - q['observed_N']) < 1e-09]
        if matches:
            z = matches[0]
            hw = float(t.ppf(1 - 0.05 / 6, z['n'] - 1) * z['sd'] / math.sqrt(z['n']) + 0.005)
            (lo, hi) = (z['mean'] - hw, z['mean'] + hw)
            (al, ah) = row['fixed_prediction_acceptance_N']
            decision = 'PASS' if al <= lo and hi <= ah else 'FAIL' if hi < al or lo > ah else 'UNRESOLVED'
            row.update(n=z['n'], n_status='INFERRED_110/11_BALANCED_ALLOCATION_NOT_EXPLICIT_PER_GROUP', external_sample_count_admitted=False, sample_SD_N=z['sd'], simultaneous95_response_mean_interval_N=[lo, hi], source_mean_gate=decision, source_mean_gate_interpretation='CONDITIONAL_ON_N10', source_origin='MEASURED_GROUP_MEAN', uncertainty='Conditional balanced n10 (total110/11 groups, no explicit per-cell count); normal iid endpoints and three-comparison Bonferroni; fixed predictor only, endpoint-fit uncertainty UNKNOWN')
        else:
            row.update(source_mean_gate='UNKNOWN', source_origin='WEIBULL_CLOSURE_EXPECTATION', uncertainty='No empirical mean/SD and parameter covariance absent; gap cannot be declared inside measurement noise')
        row['full_physical_gate'] = 'UNKNOWN'
        res.append(row)
    return res

def paired_sufficiency():
    a = np.array([7.0, 9.0, 11.0, 13.0])
    b = np.array([6.0, 8.0, 10.0, 12.0])
    c = b[::-1].copy()
    summaries = lambda x, y: np.array([x.mean(), y.mean(), x.std(ddof=1), y.std(ddof=1), (x - y).mean()])
    err = float(abs(summaries(a, b) - summaries(a, c)).max())
    assert err == 0
    intervals = []
    for y in [b, c]:
        d = a - y
        hw = float(t.ppf(0.975, 3) * d.std(ddof=1) / 2)
        intervals.append([float(d.mean()) - hw, float(d.mean()) + hw])
    return {'identity_error': err, 'marginal_summary_A': summaries(a, b).tolist(), 'marginal_summary_B': summaries(a, c).tolist(), 'state_A': {'X': a.tolist(), 'Y': b.tolist()}, 'state_B': {'X': a.tolist(), 'Y': c.tolist()}, 'paired95_interval_A': intervals[0], 'paired95_interval_B': intervals[1], 'downstream_lower_bound_difference': abs(intervals[0][0] - intervals[1][0]), 'decisions': ['POSITIVE', 'UNRESOLVED'], 'minimum_extension': 'Within-region paired covariance or SD of paired differences, alongside sample/specimen/state keys; moments suffice only under declared Gaussian mean-contrast closure', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/run_r3.py paired permutation', 'compared_quantity': 'Marginal-summary sufficiency for paired contrast interval', 'refutes_us': True}}

def gap_rows(force):
    r1 = json.loads((ROOT / 'rounds/R1_results.json').read_text())
    r2 = json.loads((ROOT / 'rounds/R2_results.json').read_text())
    rows = []
    for e in read('net')['edges']:
        g = e.get('gap', {})
        if g.get('value') is None:
            continue
        row = {'edge_id': e['id'], 'quantity': g['quantity'], 'gap': g['value'], 'unit': g['unit'], 'resolution': g['resolution_level'], 'timescale': e.get('timescale'), 'source_scope': g.get('scope'), 'source_operands': g['operands'], 'consumer_decisions': e.get('decision_ids', []), 'full_physical_decision': 'UNKNOWN', 'measurement_noise': None, 'rigorous_enclosure': 'MISSING_FOR_PHYSICAL_TRANSFER', 'external_referent': e.get('external_referent'), 'route': 'LABB', 'decided_now': "No full physical conclusion from scalar gap", 'own_design': "Request signed consumer margin, right observation support and a enclosure", 'lab_remainder': e.get('knowledge_debt', {}).get('replacement_measurement', g.get('required_measurement'))}
        id = e['id']
        if id == 'D-E-K02':
            row.update(route='EGEN_DESIGN→LABB', decided_now="40 HU_ref-grinden for global calibration falls; error after3 upprepningar134–142 HU_ref", own_design="AIR+LDPE+Teflon minimax on existing baseline; holding errors180.318 HU_ref. No multi-repeat solution against material/protokollbias.", lab_remainder="Under bilinear [1,h,r²,h*r²]-kontrakt:2 HA-levels×2 positioner ger4 kalibrerings-ROI; third held-out level at2 positioner ger2 mismatch-prov. Totalt6 ROI, precision/replikat UNKNOWNParade bone mechanics test if the consumer requires E.", uncertainty_kind='Finite observed repeat range and protocol transfer; no sigma', decision_margin={'threshold_HU_ref': 40, 'excess_HU_ref': g['value'] - 40})
        elif id in ['D-E-K17', 'D-E-K50']:
            row.update(route="AVGJORT_NU_AVGRAPPENDIX", decided_now="At the Hands0.18mm:0.65+0.60N digitiseringsbudget<1.66068N fryst tolerans; diskret gruppkurvegrind klaras", own_design="Fas/lastsignum and all power–the displacement track is preserved; more precision points are not needed for this particular frozen point gate", lab_remainder="Fysisk PDL- Law requires paired jaw/PDL-reply response and passivitet/rate; interpolation between points lacks rigorous containment", decision_margin={'threshold_N': 1.66068, 'gap_plus_two_digitization_bounds_N': g['value'] + 0.6, 'headroom_N': 1.66068 - g['value'] - 0.6}, uncertainty_kind='Conservative two +/-0.3N digitization closures; not biological/sensor uncertainty')
        elif id in ['D-E-K23', 'D-E-PULP-HEAT']:
            row.update(route='EGEN_DESIGN→LABB', decided_now='Proxy-MAE0.602251°C faller mot0.5°C. Lasergruppkontrast:2/9 positiva med simultant95%+driftbudget, villkorat n4.', own_design="9 kontraster:139 projected independent samples without drift against450 Even.2Hz/350mJ cannot be determined by:0.40°C drift oavsett n.", lab_remainder="Samma tand/dentin/tid: synkron pulpsensor + bracketing referens/blank as well as tool tracks; separate operation, location and heat source.", uncertainty_kind='Group SD, differential drift and unmatched proxy error are distinct; sensor SD UNKNOWN', decision_margin={'engineering_MAE_threshold_degC': 0.5, 'excess_degC': g['value'] - 0.5})
        elif id == 'D-E-CANAL-RELEASE':
            new = read('canal_new')
            row.update(route='EGEN_DESIGN→LABB', decided_now="Releasemasker changes virtual2mm-beslut; fruset X73 visar12/748, not independent annotator noise.", own_design="Reuse local full-cylinder/union-query on142 Image pairs; no new mask copies are counted as independent measurements.", lab_remainder="A19: local anatomy reference and independent reader on the same specimen/protokoll; releasespread does not replace true wall error.", source_dropout=new['measurement']['dropout'], uncertainty_kind='Located release revisions, not physical sigma')
        elif id == 'D-E-REPLICA-CT':
            q = next((x for x in r2['paired']['replica'] if x['system'] == 'Amann' and x['region'] == 'occlusal'))
            row.update(route='EGEN_DESIGN→LABB', decided_now='Amann MO-kontrast60.82µm; konservativt simultant95%[-19.95,141.59]µm korsar0. Ingen generell korrigeringsfaktor.', own_design="Projicerat15 paired copings for this contrast in8-contrast family at worst covariance; requesting regional paired raw values first.", lab_remainder="Separate method and conditions with four observation states; X59 kalibrering10+held10 kept for new lab systems.", conditional_interval=q['simultaneous95_worst_covariance_interval_um'], uncertainty_kind='Paired covariance UNKNOWN; pooled Pearson cannot replace region covariance')
        elif id == 'D-E-PRELOAD-HISTORY':
            row.update(route="AVGJORT_NU_AVGRAPPENDIX", decided_now="Samma dry25Ncm protocol: first–a tenth of the mean preload76.2N,95% Worst Covariance[47.264,105.136]N; minskning avgjord i samplingmodellen.", own_design="For a replicated first–tionde kontrast: projicerat6 par med stabil publicerad SD; keep sample-ID and order.", lab_remainder="Ny bro/skruv: local support reactions, actual preload history and clamping requirements; not insertionstorque/chewingtransfer.", conditional_interval=r1['preload_group_answer']['worst_covariance_pointwise95_interval_N'], uncertainty_kind='POPULATION normal mean-contrast interval, unknown covariance worst-case')
        elif id.startswith('D-E-X69'):
            row.update(route='EGEN_DESIGN→LABB', decided_now="Befintligt IOS–CBCT-pair gives surface residue per tooth, no measured physical TRE eller sigma; inga lokala2mm-decisions may be absolved from this p95.", own_design="Use whole dental regions and fixed registration dates; the same paid data already exists in Hao2023 Demo1.", lab_remainder="A01/M11: independent fiducial/ytfacit in the same region and framework, gingivareference if the consumer needs it.", uncertainty_kind='Observed p95 surface residual is not sigma, TRE or signed local clearance')
        elif id.startswith('D-E-X70-FORCE'):
            f = next((x for x in force if x['edge_id'] == id))
            row.update(route="AVGJORT_NU_AVGRAPPENDIX" if f['source_mean_gate'] == 'PASS' else 'EGEN_DESIGN→LABB', decided_now="Villkorat n10-balans (not verified per cell), fixed prediction and published mean response: " + f['source_mean_gate'] + " inom fryst0.20-logband; fit/full physical uncertainty UNKNOWN.", own_design="Measure SD/n for the exact force point; no Weibull derived mean response is recorded as an empirical agent.", lab_remainder="Test-matched support/cement/as-built/brott and separately held batch if it is a new design issue.", conditional_gate=f, uncertainty_kind=f['uncertainty'])
        elif id == 'D-E-X72-FRICTION':
            row.update(route='EGEN_DESIGN', decided_now="70.2139 – 78.2302 MPa is conditional condensing, not measured model gap. The consumer's voltage limit is missing.", own_design="For peak <=tau:tau<70.21389941372837MPa gives FAIL and tau> = 78.23015713238145 MPa, PASS is within the declared cone model; intermediate positions require local friction information.", lab_remainder="Contact force /friction on the same surface if query crosses the interval; FE / Floating Numbers/Physical enclosure still UNKNOWN .", uncertainty_kind='Conditional optimizer sandwich; rigorous physical/rounding bound absent')
        elif id == 'D-E-TANDLAST-PARTIAL-REFERENCE':
            row.update(route='EGEN_DESIGN', decided_now="7percentage points are missing coordinate mass, not measurement noise.14-dot facit can't fall16-koordinatsnormalisering.", own_design="Keep missing two FDI as unknown with sum7pp; question only jointly14-support or retrieve missing reference values.", lab_remainder="Absolut/per-patch power is needed first if an N-consumer requires it; new lab list for the normalisation gap is unjustified.", uncertainty_kind='Support mismatch, no random noise estimate')
        rows.append(row)
    assert len(rows) == 17
    return rows
GROUP_DESIGNS = {'M01': ("Oberoende CMM/lfield reference and fixed date on the same12 pilotkronor; bevara signerad geometri.", "As-built scan lacks manufacturing/scanner separation; physical limit and variance unknown."), 'M02': ("Samma12 in dry and cemented condition, local points and height/load; M01 First.", "Cunali data are: PVS/dry and cannot replace cement film; no universal60.82µm-korrektion."), 'M04': ("Samma12 without measurement-induced damage; force trace, stiffness, two blinded fracture-origin readers.", "Chen protocol groups are already measured but are not our new lab; the whole72 retains12/cell."), 'M07': ("Co-registered patch power and DIC in at least two separate load modes; select excitation modes by rank of the current compliance contract.", "Kraftskala/patchlthe number of pilots must be distinguished;4 is X71-antagande, precision UNKNOWN."), 'M08': ("Four simultaneous support reactions on the same bridge + real preload race; for reuse-sign6 projected new couples, separate question.", "Publicerat25-parse response determines preload reduction; clamping force/threshold for new bridge is missing."), 'M09': ("Measure local slip and total system displacement simultaneously at the same force; more RFARepeats don't identify slip.", "BIC and compliance separately; biological decision threshold/measuring variance UNKNOWN."), 'M10': ("Under a bilinear level/position model:4 kalibrerings-ROI vid2 levels×2 positioner, samt2 ROI on hold level; then paired the same-region bone mechanics.", "Totalt6 ROI for this contract; minimum ranking4 and2 No minimum precision certified scan/test size is known."), 'M11': ("Independent high resolution anatomy reference + realized pose; three maximum/three mandibulary scouts according to X71.", "Hao pair and release masks already measure comparisons but no absolute local truth/TRE."), 'M14': ("Keep2 Spacer modes20+held spacer10+held batch10The n can be reduced only when local film responses and variance exist.", "Samma single-spacer72- Series can't identify spacer sensitivity."), 'A01': ("Reusing Hao2023 co-registrable pair for digital surface test; add independent fiducial/gingivafacit On the same jaw.", "Existing couples have no target scanners.TRE and soft tissue facitis."), 'A02': ("Samma individ/tid: Power N, pose and EMG synchronously; excite at least as many independent loading modes as active unknown forces.", "T-scan percent or datasets with separate individuals do not identify absolute load."), 'A03': ("Exponering/cykler/lastsvans is logged on the same outcome individual; sequential design against a concrete life margin.", "Short static load missing parafunction tail; sampling length UNKNOWN."), 'A04': ("TMJ-pose and reaction synchronously on the same jaw; separate several directions for torque scale.", "Condyl track without force does not measure reaction."), 'A05': ("Two steps with the same end geometry; moment/matning/temperaturhistoria and test-matched histology.", "Temperaturmedel/sluttemperatur does not replace damage history; laser pulp data is the wrong tissue/recipient."), 'A06': ("Same hole/leg region: profile, fabric/material and moment–feeding; hold out a densitet/geometri.", "Drill torque without paired benfacit does not identify the deposit law."), 'A07': ("Same sample before load→snitt→new load, time stamp/rest time; compare signed response.", "Final equilibrium without history does not determine irreversibility."), 'A08': ('Longitudinellt sammaimplantat: ISQ/MBL/BIC+last; separat heldkohort.', "Cell response is no ISQ-measurement; no matchedhistological data in these snapshots."), 'A09': ("Keramikfraktur med material/support/exponering and censoring, same outcome individual.", "Total survival/implantatstatus does not determine material specific fracture."), 'A10': ("Samma implantat: pretestµCT+preload and documented fatigue/runout, heldbatch.", "Break load on kronor/skruv-reuse is not implant fatigue."), 'A11': ("Samma materialpar/yta/miljisland over time; local loss and response.", "Group wear without rumslig/tidsmatchning does not determine the contact consumer."), 'A12': ("Two known pressure fields + paired intakelio scans and independent load reference on the same mucosa.", "Two unloaded scans provide no stiffness; load fields must be measured, count precision UNKNOWN."), 'A13': ("The same crown process log and signed before/after fields; M01/M02 sammaID.", "More random end scans do not identify sintering/set cause."), 'A14': ("Korsad orientering×fully oriented layer height, frozen process/lot.", "A modified CAD surface is not a manufactured Accuracy measurement."), 'A15': ("Samma Ti-prov XCT pre-test; yta/process and mechanics; full-blown.", "Geometridata/fra tight fatigue table cannot be paired without a test ID."), 'A16': ("A hash-locked export sample through real CAM to tillverkad/metrologiskt measured portion; then holding test.", "Importerparity without real manufacturing does not measure process compatibility."), 'A17': ("Separate specimen/scanner/lot holdout med samma frystaθ and the same local portobservables.", "Refit on the holder does not provide generalisation facit."), 'A18': ("Hela72-protokollets batch/support/provtemperatur, full metrology and fracture testing; new quote for extra60 setup measurements.", "12-scouter is not enough for the wholeQ or originalets12/cell."), 'A19': ("Independent reader+repeated scan+anatomy reference in fixed frame; identify bias and regional covariance separately.", "0 identifierade independentreaderpar iX73; release audits are not sigma."), 'A20': ("Kraft–displacement loop per jaw/tand/PDL vid flera rates, paratgeometri.", "0.65N group curves are under the point gate but full PDL-materialmodell/passivitet remains."), 'A21': ("The same bone implant test: local grinding and BIC vid samma tid; histologi efter synkron last.", "Total system movement and contact surface percentage is not local slip.")}

def group_rows():
    c = read('contracts')
    rows = []
    physical = {k: v for (k, v) in c['edges'].items() if v.get('status') == 'BLOCKED_ON_PHYSICAL_MEASUREMENT'}
    for (gid, g) in c['groups'].items():
        es = [k for (k, v) in physical.items() if gid in v.get('measurement_ids', [])]
        if not es:
            continue
        (design, reason) = GROUP_DESIGNS[gid]
        rows.append({'id': gid, 'name': g['name'], 'route': 'LABB_VILLKORAT', 'gap': 'UNKNOWN', 'decided_now': "Full target relationship not yet settled; see limited source answers in the gap table", 'own_design': design, 'remaining_reason': reason, 'edge_ids': es, 'touch_count': len(es), 'alone_type_ready_edges': [k for k in es if physical[k]['measurement_ids'] == [gid]], 'observed_decision_flips_per_specimen': None, 'actual_cost': None, 'operator_hour_planning_interval': g.get('operator_hours_interval'), 'cost_resolution': 'PHENOMENOLOGICAL' if g.get('operator_hours_interval') else 'UNKNOWN', 'minimum_n_for_physical_precision': None, 'replacement_measurement': g['measurement'], 'specimen_contract': g['specimen'], 'prerequisites': g.get('prerequisites', []), 'resolution': g.get('resolution_level'), 'timescales': sorted(set((physical[k]['timescale'] for k in es))), 'lab_activation_gate': 'Named decision margin + matched observable + quoted pilot cost required; use sequential refinement only if its uncertainty still straddles threshold', 'additional_model_and_holdout_requirements': {k: physical[k]['additional_prerequisites'] for k in es}})
    assert len(rows) == 30
    return (rows, physical)

def close_prerequisites(ids, groups):
    selected = set(ids)
    while True:
        expanded = selected | {p for g in selected for p in groups[g].get('prerequisites', [])}
        if expanded == selected:
            return selected
        selected = expanded

def ready_edges(ids, edges):
    return sorted((k for (k, v) in edges.items() if set(v.get('measurement_ids', [])) <= set(ids)))

def lab_frontier(physical):
    groups = read('contracts')['groups']
    x71 = {x['id']: x for x in read('x71_measurements')}
    plans = []
    candidates = [[gid] for gid in groups if gid in x71 and any((gid in v.get('measurement_ids', []) for v in physical.values()))]
    candidates += [['M01', 'M02', 'M04'], ['M01', 'M02', 'M07'], ['M01', 'M02', 'M08']]
    seen = set()
    for candidate in candidates:
        ids = close_prerequisites(candidate, groups)
        identity = tuple(sorted(ids))
        if identity in seen:
            continue
        seen.add(identity)
        hours = [sum((groups[i]['operator_hours_interval'][j] for i in ids)) for j in [0, 1]]
        count = (12 if set(ids) & {'M01', 'M02', 'M03', 'M04'} else 0) + sum((max(x71[i]['specimens'], x71[i]['new_specimens']) for i in ids if i not in {'M01', 'M02', 'M03', 'M04'}))
        ready = ready_edges(ids, physical)
        ctrl = []
        for (eid, v) in physical.items():
            missing = False
            for mid in v.get('measurement_ids', []):
                if mid not in ids:
                    missing = True
            if not missing:
                ctrl.append(eid)
        assert sorted(ctrl) == ready
        plans.append({'ids': sorted(ids), 'operator_h_assumed': hours, 'specimens_including_shared_cohorts': count, 'type_ready_edges': ready, 'type_ready_count': len(ready), 'rate_upper_proxy_per_specimen': len(ready) / count, 'potential_queries_per_upper_hour': len(ready) / hours[1], 'real_decision_flips_per_specimen': None, 'real_decision_flips_per_hour': None, 'possible_reversal_count_interval': [0, len(ready)], 'uncertainty': 'Structural upper proxy; actual reversals and variance UNKNOWN', 'physical_edges_closed': 0, 'cost_status': 'PHENOMENOLOGICAL_X71_NOT_QUOTE'})
    plans.sort(key=lambda x: (-x['potential_queries_per_upper_hour'], x['operator_h_assumed'][1], x['ids']))
    pilot = next((p for p in plans if p['ids'] == ['M01', 'M02', 'M04']))
    original = read('coverage')['plans']['X71_crown_pilot']['all_listed_acquisition_types_present']
    assert pilot['type_ready_edges'] == sorted(original)
    return {'priced_finite_plans': plans, 'pilot': pilot, 'rank_rule': 'Descending structural type-ready queries/upper assumed operatorhour; tie cost then IDs. Not actual flips or calibrated value of information.', 'unpriced_groups': [k for k in groups if k.startswith('A')], 'conditional_week_policy': 'For original72 crowns, retain M01→M02→M04 on same12. No full72 cost available. If scanner/absolute length truth missing, acquire M01 first.', 'observed_lab_flips': None, 'global_optimum': 'UNKNOWN', 'new_lab_acquisitions': 'NOT_RUN', 'comparison_to_X71': 'X71 five D questions differ from R3 two AND-type-ready edges. No empirically justified first-week reordering from these source data.', 'comparison_to_R3': 'R3 M01 touches12 but alone satisfies one listed edge; counted correctly. M11 supplies sole type on4 but specimen/precision/model requirements remain.'}

def main():
    verify_freeze('PREREG_R3.json')
    verify_freeze('FROZEN_PREDICTIONS.json')
    f = force_intervals()
    (gr, physical) = group_rows()
    gaps = gap_rows(f)
    out = {'claim_type': 'capability', 'force_intervals': f, 'paired_sufficiency': paired_sufficiency(), 'numeric_gap_rows': gaps, 'physical_group_rows': gr, 'physical_edge_obligations': physical, 'lab_frontier': lab_frontier(physical), 'coverage': {'numeric_gap_edges': len(gaps), 'physical_groups': len(gr), 'physical_edges': len(physical), 'total_open_selected_edges': len(gaps) + len(physical)}, 'dropout': {'net_edges': len(read('net')['edges']), 'selected_open_edges': len(gaps) + len(physical), 'excluded_open_nonphysical_edges': [k for (k, v) in read('contracts')['edges'].items() if v.get('status') != 'BLOCKED_ON_PHYSICAL_MEASUREMENT'], 'nonphysical_exclusion_reason': 'K10/K11/K18 computation/scope; K51 installation; K52 licence. Never convert to laboratory work.', 'independent_canal_reader_pairs': 0, 'independent_canal_pair_rejection_fraction': 1.0, 'calibration_candidates': 35, 'calibration_nonmonotonic_rejections': sum((not x['monotonic'] for x in json.loads((ROOT / 'rounds/R2_results.json').read_text())['calibration']['all_candidates'])), 'laser_sample_count': 'n4 inferred balance; all18 group readings retain that conditional status'}, 'full_physical_decisions_certified': 0, 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    save('rounds/R3_results.json', out)
    (ROOT / 'rounds/HANDOFF_R3.md').write_text('R3: complete17 gap/30 group/51 physical-edge inventory, force source mean uncertainty, paired covariance exact sufficiency witness, and conditional cost frontier. None of the unknown physical chains closes. Next construction is a matched drift/reference or HA-position test with genuine lab metadata; do not buy repeats before ruling out differential bias.\n')
    state('R3_COMPLETE', 'COMPLETE_DECISION_TABLE_PHYSICAL_UNKNOWN_RETAINED', 'Run independent source/raw checks, actual rejecting mutations, package one-command demo and graph feedback')
    print('R3 coverage', out['coverage'], 'force mean gates', [q['source_mean_gate'] for q in f])
if __name__ == '__main__':
    main()
