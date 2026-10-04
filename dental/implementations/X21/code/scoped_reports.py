"""A region-indexed report observable resolves legitimate mixed anterior/lateral findings."""
from dental_release.paths import expand as _release_expand
import json, re, datetime, time
from pathlib import Path
from contact import P, sha, write

def prereg():
    r = dict(round='R6', claim_type='capability', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), capability='Represent anterior, lateral/posterior and unspecified crossbite statements separately so one report can say lateral absent and anterior present without a false contradiction.', obstacle='R5 scalar negation injection cannot reject global presence in five reports that explicitly combine anterior presence with lateral absence. Prefix word bags also confuse lateral incisors with lateral crossbite.', changed_operation='Use the directly adjacent modifier of each crossbite mention to key the observation by region. Preserve each source clause, memberhash, polarity and unspecified scope. Region-specific injected opposite polarity must be rejected against its own source clause. Map named anterior teeth to corresponding signed geometric candidate margins; leave anatomical validityUNKNOWN.', consumer='Researcher comparing tooth contact/opposition with the actual clinical finding in the same region.', gates={'literal_no_lateral_controls': 125, 'all_region_positive_injections_rejected': True, 'mixed_scope_records_preserved': True, 'no_false_conflict_for_anterior_present_lateral_absent': True}, scope_rule='Adjacent anterior -> anterior; lateral/posterior -> lateral; otherwise unspecified. Nonadjacent lateral-incisor noun does not define crossbite region.', external_referent=dict(kind='independent_measurement', locator=_release_expand('@DENTAL_DATA_ROOT@/geometry/Bite2Text/Bite2Text.zip::*/reports_ios_en/*.txt, raw/CORRECTED_REPORTS_R5.json source SHA256; five dual-region reports @DENTAL_CASE_ID@/18_1855,@DENTAL_CASE_ID@/18_821,@DENTAL_CASE_ID@/18_1801,@DENTAL_CASE_ID@/18_1811,@DENTAL_CASE_ID@/18_1780'), compared_quantity='Explicit source polarity of region-indexed crossbite statements, not physical tooth diagnosis', refutes_us=True), strongest_equally_informed_control='Literal source clauses; no algorithm superiority claim', falsifiers=['A no-lateral clause maps to lateral presence', 'An anterior crossbite adjacent modifier is miskeyed by lateral-incisor phrase', 'An opposite region-specific injection is accepted'], cost=dict(preparation='Reuse hashed R5 source clauses; no data duplication', fit='none', discovery='Changed observation representation after mixed-region negative control revealed scalar insufficiency', validation='All125 literal negative region statements checked; five dual-region geometries retained', queries='Region-indexed lookup only', fallback='Unspecified scopeUNKNOWN, independent anatomical truthUNKNOWN'))
    write(P / 'PREREG_R6.json', r)
    (P / 'PREREG_R6.json.sha256').write_text(sha(P / 'PREREG_R6.json') + '\n')
    write(P / 'DECOMPOSITION_R6.json', dict(idea=r['capability'], leaves=[dict(name='region-specific statements', status='EXTERNALLY_MEASURED', relation='Literal no lateral crossbites and anterior crossbite present in same source report.', stopping_argument='Independent published text, not anatomy truth.'), dict(name='statement identity', status='DERIVED_UNDER_ASSUMPTIONS', relation='Negative(region lateral) and positive(region anterior) have different proposition keys; no logical contradiction.', stopping_argument='Propositional identity depends on explicit adjacent modifier, not a statistical model.'), dict(name='geometry-to-region linkage', status='UNKNOWN', relation='Named anterior tooth IDs map to X11 opposing point order.', stopping_argument='FDI/point validity on target remains unmeasured.')]))

def run():
    st = time.perf_counter()
    reports = json.load(open(P / 'raw/CORRECTED_REPORTS_R5.json'))
    op = json.load(open(P / 'raw/OPPOSITION_FEATURES_R2.json'))
    out = {}
    controls = []
    mixed = []
    for (c, rs) in reports.items():
        out[c] = []
        for r in rs:
            propositions = []
            occurrences = {}
            for clause in r['crossbite_scope_clauses']:
                s = clause['clause']
                mentions = list(re.finditer('cross[ -]?bites?', s, re.I))
                index = occurrences.get(s, 0)
                occurrences[s] = index + 1
                mention = mentions[min(index, len(mentions) - 1)]
                modifier = re.search('\\b(anterior|posterior|lateral)\\s*$', s[:mention.start()], re.I)
                region = modifier.group(1).lower() if modifier else 'unspecified'
                region = 'lateral' if region == 'posterior' else region
                polarity = 'absent' if clause['explicitly_negated'] else 'present'
                propositions.append(dict(region=region, polarity=polarity, clause=s))
            states = {}
            for region in ['anterior', 'lateral', 'unspecified']:
                values = {p['polarity'] for p in propositions if p['region'] == region}
                states[region] = next(iter(values)) if len(values) == 1 else 'CONFLICT' if len(values) > 1 else 'UNKNOWN'
            record = dict(member=r['member'], source_sha256=r['sha256'], states=states, propositions=propositions, physical_status='UNKNOWN')
            out[c].append(record)
            literal = [p for p in propositions if re.search('\\bno lateral cross[ -]?bites?', p['clause'], re.I)]
            for p in literal:
                controls.append(dict(case_id=c, member=r['member'], region=p['region'], reference=p['polarity'], injected='present', rejected=p['region'] == 'lateral' and p['polarity'] == 'absent'))
            if states['lateral'] == 'absent' and states['anterior'] == 'present':
                named = sorted({int(k) for p in propositions if p['region'] == 'anterior' for k in re.findall('\\b[1-4][1-8]\\b', p['clause'])})
                geo = [t for t in op.get(c, {}).get('teeth', []) if t['upper_fdi'] in named or t['lower_corresponding_fdi'] in named]
                mixed.append(dict(case_id=c, source_report=record, named_teeth=named, geometry_opposition_candidates=geo, contradiction=False, explanation='Same report describes different spatial regions, so global present is compatible with lateral absence.'))
    write(P / 'raw/SCOPED_REPORTS_R6.json', out)
    write(P / 'raw/SCOPED_NEGATION_CONTROLS_R6.json', controls)
    write(P / 'raw/MIXED_SCOPE_GEOMETRY_R6.json', mixed)
    gates = dict(n_literal_controls=len(controls) == 125, all_injected_region_values_rejected=all((c['rejected'] for c in controls)), mixed_scope_preserved=len(mixed) >= 5, no_mixed_report_contradiction=all((not r['contradiction'] for r in mixed)))
    result = dict(round='R6', claim_type='capability', external_referent=json.load(open(P / 'PREREG_R6.json'))['external_referent'], reports=sum(map(len, out.values())), literal_no_lateral_controls=len(controls), mixed_anterior_positive_lateral_negative_reports=len(mixed), gates=gates, primary_gate='PASS' if all(gates.values()) else 'FAIL', all_wrong_region_labels_rejected=all((c['rejected'] for c in controls)), wall_s=time.perf_counter() - st, remaining_unknown='Unspecified report scope, target point/tooth validity, physical contact/pressure. R5 global classifier metrics descriptive; R6 does not create clinical diagnostic accuracy.', original_failed_R5_control_preserved=True)
    write(P / 'raw/RESULTS_R6.json', result)
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    import sys
    prereg() if sys.argv[1] == 'prereg' else run()
