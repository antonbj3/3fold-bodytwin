"""Reproduce the frozen X43 status proposal. Standard library only; no simulations."""
from dental_release.paths import expand as _release_expand
import copy, datetime, hashlib, json, re, time, resource
from pathlib import Path
from collections import Counter
from xml.sax.saxutils import escape
HERE = Path(__file__).resolve().parent
ROOT = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))
CATEGORIES = ['PRESENT', 'PARTIAL', 'REFUTED', 'IN_PROGRESS', 'RESEARCH', 'MISSING']

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    return json.loads(path.read_text())

def write(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')

def baseline():
    text = (HERE / 'raw/source_snapshot/notes/NORTHSTAR_SPEC.md').read_text()
    rows = {}
    for (line_no, line) in enumerate(text.splitlines(), 1):
        if not re.match('^\\| K\\d{2} \\|', line):
            continue
        fields = line.split('|')
        k = fields[1].strip()
        raw = fields[3].strip()
        status = next((s for s in CATEGORIES if raw.startswith(s)), None)
        if status is None:
            raise ValueError(f'Unknown baseline status: {k} {raw}')
        rows[k] = {'k': k, 'chain': fields[2].strip(), 'status': status, 'original_status_text': raw, 'source_line': line_no, 'mixed_refutation': raw.startswith('PARTIAL/REFUTED')}
    if set(rows) != {f'K{x:02}' for x in range(1, 53)}:
        raise ValueError('Baseline must contain exactly K01-K52')
    return rows

def validate(rows, refs, manifest, requirements):
    errors = []
    if len(rows) != 52 or {r['k'] for r in rows} != {f'K{x:02}' for x in range(1, 53)}:
        errors.append('K_COVERAGE')
    known = {r['relative_path']: r for r in manifest['sources']}
    for r in rows:
        if r['proposed_status'] not in CATEGORIES:
            errors.append('INVALID_STATUS')
        if not r.get('evidence_keys'):
            errors.append('MISSING_EVIDENCE')
        for key in r.get('evidence_keys', []):
            e = refs.get(key)
            if not e or e['path'] not in known:
                errors.append('UNLOCATED_EVIDENCE')
                continue
            snap = HERE / known[e['path']]['snapshot']
            if not 0 < e['line'] <= len(snap.read_text().splitlines()):
                errors.append('BAD_LINE')
        if r['proposed_status'] == 'PRESENT' and r['old_status'] != 'PRESENT' and (not r['full_chain_reviewed']):
            errors.append('UNSUPPORTED_PRESENT')
        if r['proposed_status'] != r['old_status'] and (not any((refs[k]['kind'] == 'independent_review' for k in r['evidence_keys'] if k in refs))):
            errors.append('UNREVIEWED_STATUS_CHANGE')
        if not r.get('resolution') or not r.get('remaining_requirement'):
            errors.append('MISSING_SCOPE_OR_RESOLUTION')
        if 'PHENOMENOLOGICAL' in r['resolution'] and (not r.get('phenomenological_debt')):
            errors.append('UNPAID_CLOSURE')
    if set(requirements) != {r['k'] for r in rows if r['proposed_status'] == 'MISSING'}:
        errors.append('MISSING_ACQUISITION_PLAN')
    counts = Counter((r['proposed_status'] for r in rows))
    if sum(counts.values()) != 52:
        errors.append('COUNT_SUM')
    return sorted(set(errors))

def validate_summary(summary, rows):
    actual = Counter((r['proposed_status'] for r in rows))
    errors = []
    if sum(summary.values()) != 52:
        errors.append('COUNT_SUM')
    if any((summary.get(k) != actual[k] for k in CATEGORIES)):
        errors.append('COUNT_BY_ROWS')
    return errors

def source_check(manifest):
    rows = []
    for r in manifest['sources']:
        snap = HERE / r['snapshot']
        original = Path(r['path'])
        rows.append({'path': r['relative_path'], 'snapshot_matches': snap.exists() and sha(snap) == r['sha256'], 'live_matches': original.exists() and sha(original) == r['sha256']})
    return rows

def links(keys, refs):
    return ', '.join((f"[{k}]({ROOT / refs[k]['path']}:{refs[k]['line']})" for k in keys))

def render_chart(old, new):
    width = 850
    height = 410
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="850" height="410" fill="#fafbfc"/>', '<g font-family="sans-serif" font-size="15" fill="#17212b">', '<text x="24" y="30" font-size="21">Northstar: 52 administrative chain statuses</text>', '<text x="24" y="57">Gray: frozen specification. Blue: X43 proposal. No added physical validation.</text>']
    maxv = max(max(old.values()), max(new.values()), 1)
    for (i, s) in enumerate(CATEGORIES):
        y = 90 + i * 43
        scale = 590 / maxv
        parts += [f'<text x="24" y="{y + 12}">{escape(s)}</text>', f'<rect x="150" y="{y}" width="{old[s] * scale}" height="13" fill="#9da7af"/>', f'<text x="{160 + old[s] * scale}" y="{y + 12}">{old[s]}</text>', f'<rect x="150" y="{y + 16}" width="{new[s] * scale}" height="13" fill="#2776a8"/>', f'<text x="{160 + new[s] * scale}" y="{y + 28}">{new[s]}</text>']
    parts += ['<text x="24" y="374">K33 counts as PARTIAL; its refuted strength component is retained separately.</text>', '<text x="24" y="397">Inherited scoped PRESENT components do not approve the full patient chain.</text>', '</g></svg>']
    (HERE / 'status_counts.svg').write_text('\n'.join(parts) + '\n')

def main():
    start = time.monotonic()
    stamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    reg = read(HERE / 'PREREG_X43_v1.json')
    cur = read(HERE / 'raw/CURATED_EVIDENCE.json')
    manifest = read(HERE / 'SOURCE_MANIFEST.json')
    base = baseline()
    rows = copy.deepcopy(cur['chains'])
    refs = cur['references']
    for r in rows:
        b = base[r['k']]
        r.update(old_status=b['status'], original_status_text=b['original_status_text'], chain=b['chain'], status_changed=r['proposed_status'] != b['status'], source_line=b['source_line'], review_state='PENDING_INDEPENDENT_REVIEW', mixed_refutation_preserved=b['mixed_refutation'])
        r['evidence'] = [dict(refs[k]) for k in r['evidence_keys']]
        r['external_referent'] = {'kind': 'external_review', 'locator': [str(ROOT / refs[k]['path']) + f":{refs[k]['line']}" for k in r['evidence_keys'] if refs[k]['kind'] == 'independent_review'], 'compared_quantity': 'Scoped reviewed delivery versus ' + b['chain'], 'refutes_us': False} if any((refs[k]['kind'] == 'independent_review' for k in r['evidence_keys'])) else None
        r['uncertainty_kind'] = 'scope limitation, not a newly estimated probability'
    errors = validate(rows, refs, manifest, cur['missing_requirements'])
    source_checks = source_check(manifest)
    frozen = read(HERE / 'FROZEN_PREDICTIONS.json')
    if sha(HERE / 'raw/CURATED_EVIDENCE.json') != frozen['curated_evidence_sha256']:
        errors.append('CURATED_MAP_DRIFT')
    if sha(HERE / 'SOURCE_MANIFEST.json') != frozen['source_manifest_sha256']:
        errors.append('SOURCE_MANIFEST_DRIFT')
    lock = read(HERE / 'DEMO_CODE_LOCK.json')
    if any((sha(HERE / f) != h for (f, h) in lock['files'].items())):
        errors.append('DEMO_CODE_DRIFT')
    if any((not x['snapshot_matches'] or not x['live_matches'] for x in source_checks)):
        errors.append('SOURCE_DRIFT')
    hash_checks = read(HERE / 'raw/REVIEW_HASH_CHECKS.json')
    live_bindings = []
    for r in hash_checks:
        f = Path(r['result_file'])
        actual = sha(f) if f.is_file() else None
        live_bindings.append(dict(r, actual_sha256=actual, status='MATCH' if actual == r['expected_sha256'] else 'DRIFT'))
    if any((x['status'] != 'MATCH' for x in live_bindings)):
        errors.append('REVIEW_RESULT_DRIFT')
    old = Counter((b['status'] for b in base.values()))
    new = Counter((r['proposed_status'] for r in rows))
    old = {s: old[s] for s in CATEGORIES}
    new = {s: new[s] for s in CATEGORIES}
    errors.extend(validate_summary(new, rows))
    if new != frozen['expected_counts']:
        errors.append('FROZEN_COUNTS_DRIFT')
    if {r['k']: r['proposed_status'] for r in rows} != frozen['expected_status_per_K']:
        errors.append('FROZEN_STATUS_DRIFT')
    mutations = []
    bad = copy.deepcopy(rows)
    bad.append(copy.deepcopy(rows[0]))
    mutations.append({'fault': '53 rows instead of 52', 'errors': validate(bad, refs, manifest, cur['missing_requirements'])})
    bad = copy.deepcopy(rows)
    bad[7]['evidence_keys'] = ['NOT_A_SOURCE']
    mutations.append({'fault': 'K08 missing review locator', 'errors': validate(bad, refs, manifest, cur['missing_requirements'])})
    bad = copy.deepcopy(rows)
    bad[46]['proposed_status'] = 'PRESENT'
    mutations.append({'fault': 'K47 full-chain PRESENT unsupported by digital export only', 'errors': validate(bad, refs, manifest, cur['missing_requirements'])})
    bad = copy.deepcopy(rows)
    bad[7]['evidence_keys'] = ['SPEC']
    mutations.append({'fault': 'K08 promotion with no new reviewed evidence', 'errors': validate(bad, refs, manifest, cur['missing_requirements'])})
    bad = copy.deepcopy(rows)
    bad[15]['proposed_status'] = 'PARTIAL'
    mutations.append({'fault': 'K16 status changed but missing-chain plan not updated', 'errors': validate(bad, refs, manifest, cur['missing_requirements'])})
    forged = copy.deepcopy(new)
    forged['MISSING'] += 1
    mutations.append({'fault': 'summary MISSING forged from 7 to 8 and total from 52 to 53', 'errors': validate_summary(forged, rows)})
    altered = copy.deepcopy(manifest)
    altered['sources'][0]['sha256'] = '0' * 64
    mutations.append({'fault': 'source hash substituted with zero hash', 'errors': ['SOURCE_DRIFT'] if not source_check(altered)[0]['snapshot_matches'] else []})
    if not all((m['errors'] for m in mutations)):
        errors.append('MUTATION_SURVIVED')
    validation = {'validated_at': stamp, 'baseline_errors': errors, 'mutation_probes': mutations, 'all_bad_inputs_rejected': all((m['errors'] for m in mutations)), 'source_checks': source_checks, 'review_result_checks': live_bindings, 'note': 'Checks establish bookkeeping/version invariants. They do not automate semantic interpretation of reviewer scope or constitute independent scientific review.'}
    write(HERE / 'VALIDATION.json', validation)
    if errors:
        raise SystemExit('FROZEN GATE FAIL: ' + ','.join(errors))
    active = read(HERE / 'raw/source_snapshot/results/DEMO48_PACKAGE/PACKAGE_ACTIVE_INDEX.json')
    demos = read(HERE / 'raw/source_snapshot/results/DEMO48_PACKAGE/demos.json')
    mismatched = set(active['active_demo_ids']) - {x['id'] for x in demos['demos']}
    identities = sorted({(r['result_file'], r['expected_sha256']) for r in live_bindings})
    coverage = []
    for (path, h) in identities:
        job = Path(path).parent.name
        review_paths = {r['review'] for r in live_bindings if r['result_file'] == path and r['expected_sha256'] == h}
        ks = [r['k'] for r in rows if review_paths & {refs[k]['path'] for k in r['evidence_keys']}]
        coverage.append({'job': job, 'result_file': path, 'result_sha256': h, 'review_paths': sorted(review_paths), 'exact_review_reference_Ks': ks, 'disposition': 'DIRECT_K_BINDING' if ks else 'CONTEXT_ONLY_OR_OUTSIDE_K; see coverage notes; not counted as an independent promotion'})
    write(HERE / 'DELIVERY_COVERAGE.json', {'result_identities': coverage, 'note': 'Direct bindings use exact REVIEW JSON. Aggregate XREVIEW documents support other rows; unlisted K is not rejection of the scientific result.'})
    result = {'schema': 'northstar-status-proposal/v1', 'lane': 'X43-northstar-status', 'claim_type': 'information_link', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'created_at': stamp, 'prereg_sha256': sha(HERE / 'PREREG_X43_v1.json'), 'source_manifest_sha256': sha(HERE / 'SOURCE_MANIFEST.json'), 'frozen_curated_map_sha256': sha(HERE / 'raw/CURATED_EVIDENCE.json'), 'external_referent': {'kind': 'external_review', 'locator': [str(ROOT / 'results/LANE_AUDIT_3B_SOL/RESULTS.md'), *[str(ROOT / refs[k]['path']) for k in ['REST', 'CLIN', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'G1', 'G2']]], 'compared_quantity': 'Scope, review class and principal quantity of delivered suboperations against K01–K52 capability contracts', 'refutes_us': False}, 'old_counts': old, 'proposed_counts': new, 'status_changes': [{'k': r['k'], 'old': r['old_status'], 'proposed': r['proposed_status']} for r in rows if r['status_changed']], 'chains': rows, 'missing_chains': cur['missing_requirements'], 'evidence_census': {'frozen_documents': len(manifest['sources']), 'review_hash_bindings': len(live_bindings), 'matching_review_bindings': sum((x['status'] == 'MATCH' for x in live_bindings)), 'unique_reviewed_result_identities': len(identities), 'review_report_glob_count': len(list((ROOT / 'results').glob('*XREVIEW*/XREVIEW.md'))), 'legacy_demos_json_rows': len(demos['demos']), 'active_demo_count': active['active_count'], 'publication_candidates': active['publication_candidate_ids'], 'optional_audit_tools': active['optional_audit_tools'], 'active_ids_not_in_legacy_json': sorted(mismatched)}, 'attrition': {'K_candidates': 52, 'K_dropped': 0, 'review_bindings_missing_or_drift': 0, 'unsupported_PRESENT_promotions': 0, 'note': 'Suboperation evidence is retained with its limitations; it is not deleted because whole-chain promotion is denied. Dataset-specific losses appear in README_DEMO and NORTHSTAR_STATUS_PROPOSAL.'}, 'validation': {'positive_gate': 'PASS', 'bad_inputs_rejected': len(mutations), 'bad_inputs_tested': len(mutations), 'canonical_spec_unchanged': True}, 'full_cost': {'fit': 'None', 'physical_measurements': 0, 'new_datasets': 0, 'questions_to_user': 0, 'numerical_search': 'None; semantic status mapping is manually curated', 'inherited_review_cost': 'Not rerun or charged as new independent experiments', 'preparation_and_reading_wall_seconds': 'UNKNOWN; tool/model reasoning not precisely timed', 'elapsed_from_prereg_to_current_run_s': (datetime.datetime.now(datetime.timezone.utc) - datetime.datetime.fromisoformat(reg['frozen_at'])).total_seconds(), 'census_validation_render_wall_s': time.monotonic() - start, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'fallback': 'Frozen source snapshot, inherited status for unsupported changes; coordinator resolves remaining semantic decisions'}, 'not_established': ['New physical prediction accuracy', 'Clinical validity', 'Manufactured specimen validation', 'Automatic semantic status inference', 'Scientific admission or source graph mutation']}
    write(HERE / 'results.json', result)
    render_chart(old, new)
    summary = ['# Proposed Northstar status from reviewed deliveries', '', f"""All 52 K chains now have traceable status proposals. {len(result['status_changes'])} status labels change; no new row is promoted to PRESENT. External reference for this status revision consists of other lanes version-bound reviews, not a new physical measurement.""", '', f"""Snapshot: {stamp}. Proposals only, **PENDING_INDEPENDENT_REVIEW**. NORTHSTAR_SPEC.md unchanged; SHA256 `{sha(ROOT / 'notes/NORTHSTAR_SPEC.md')}`.""", '', 'PRESENT retained for seven previously scoped contracts. This does not retrospectively review every old run or claim closure of each broad patient chain. PARTIAL denotes an executed/reviewed suboperation with remaining requirements. REFUTED retains a controlling negative prediction. RESEARCH denotes a construction without a complete shared dental consumer. MISSING means this exact chain is unsupported by this frozen selection, even when adjacent tools exist.', '', '## Counted totals', '', '| Status | Frozen K table | Proposal | Difference |', '|---|---:|---:|---:|']
    for s in CATEGORIES:
        summary.append(f'| {s} | {old[s]} | {new[s]} | {new[s] - old[s]:+d} |')
    summary += ['| **Total** | **52** | **52** | **0** |', '', 'Old K33 PARTIAL/REFUTED counts unambiguously as PARTIAL while its refuted strength component remains. Specification Table 1b predates sections 4/5 and must be recalculated from K rows rather than copied.', '', '![Counted statuses](status_counts.svg)', '', '## Changed rows', '', '| K | From -> proposal | Reviewed support and limitation |', '|---|---|---|']
    for r in rows:
        if r['status_changed']:
            summary.append(f"| {r['k']} | {r['old_status']} → **{r['proposed_status']}** | {r['rationale']} {links(r['evidence_keys'], refs)} |")
    summary += ['', '## Proposal for every K chain', '', 'Resolution belongs to the evidence quantity. K, demo and review record counts are administrative. PHENOMENOLOGICAL denotes explicit model debt, replaced by the measurement in the last column. A locally signed gap or canal wall must not be replaced by a population mean.', '', '| K and contract | Current -> proposed status | Delivery/review class and principal value | Rationale / missing requirements |', '|---|---|---|---|']
    for r in rows:
        classes = '; '.join(dict.fromkeys((refs[k]['review_class'] for k in r['evidence_keys'])))
        src = [f"""[{r['k']}-requirement]({ROOT / 'notes/NORTHSTAR_SPEC.md'}:{r['source_line']})""", links(r['evidence_keys'], refs)]
        summary.append(f"""| **{r['k']}** — {r['chain'].replace('|', '/')} | {r['old_status']} -> **{r['proposed_status']}** | {classes}. {r['main_number']} [{', '.join(r['resolution'])}]. {'; '.join(src)} | {r['rationale']} **Remaining:** {r['remaining_requirement']} |""")
    summary += ['', '## Still MISSING: concrete data, code and measurement', '', '| Chain | Data | Code | Decisive measurement / control |', '|---|---|---|']
    for (k, v) in cur['missing_requirements'].items():
        summary.append(f"| {k} | {v['data']} | {v['code']} | {v['measurement']} |")
    summary += ['', '## Source selection, attrition and comparison', '', f"""Deterministic census of {len(manifest['sources'])} small frozen input files: all {result['evidence_census']['review_report_glob_count']} immediate `results/*XREVIEW*/XREVIEW.md` files, package README/demos/active index, coordinator log and AUD3b older Sol review. {len(live_bindings)}/{len(live_bindings)} review/result hash bindings match, representing {len(identities)} unique result identities. Duplicate decision formats, copied package reviews and multiple K consumers do not count as new experiments. No K records excluded (0/52). No result bindings excluded for hash drift (0/{len(live_bindings)}).""", '', f"""`demos.json` has {len(demos['demos'])} records from earlier packaging. `PACKAGE_ACTIVE_INDEX.json` has {active['active_count']} active demos; X39 is separate as a private publication candidate and X35_STATISTICS is an optional audit tool. {len(mismatched)} active demo IDs are absent from older demos.json; GENCAD_V2 is also archived and aliased to V3. Proposed reference for active package count is therefore the later active index. Package count must not count physically validated capabilities.""", '', 'Important data attrition (different denominators retained):', '', '- X38 rejects 4/15 source meshes (26.67%) under numeric triangle policy; 55/55 covers only eleven accepted meshes. Region labels synthetic.', '- X25: 13/36 censored, including four full runouts and nine earlier stops; Brier compares only 27 known endpoints. Three survivors at 100 N are a separate configuration, not all four runouts.', '- X36: 64/148 cement cells and 19/49 crown groups outside primary models. Only 16/30 crown groups in two held-out studies scoreable; other 14 UNKNOWN. Not a whole-corpus design rule.', '- X30: 140/423 candidate cases rejected; 283 accepted for slot selection, while 3449 untruncated tooth rows cover 282 CT scans. Another 989/4528 slots and 90/3539 measured tooth rows fail. These denominators must not be equated.', '- IL3: 11/22 link proposals rejected; six duplicates, four excessive interpretations, one missing modulus unit. After duplicates: 5/16. Two of 35 pore files lack units/correct quantity. L05/L06 external CTF deposition unresolved.', '- GenCAD v3: 22/24 larger volumetric domains export, but all 22 fail exterior gate. R6 exports 8/8, but 0/8 pass 0.35 mm; digital roof benchmark and anatomical full crown track are separate.', '', 'Control outcome: exact 52 frozen K rows recalculated and compared with curated rows; all total, locator and version gates PASS. Seven actually changed bad bookkeeping inputs rejected in `VALIDATION.json`. Machine checks cannot decide correctness of a semantic review mapping; coordinator independent reading remains. No algorithmic advantage claimed.', '', '## Deliveries requiring their own contracts or narrower language', '', 'X17 sleep apnea, X19/X27/X37 aligners and X22 synthetic caries broaden the project but lack exact complete K contracts here. They require scoped coverage proposals; positive geometry/benchmark results do not upgrade patient load, biological healing or contact force. X29 bridge connector group measurement provides K36/K49 context without a MEASURED connector law. X35 statistics improve uncertainty reporting without new specimens or physical transfer. X39 is a private manuscript artifact reusing a cohort with publication_ready=false. See `DELIVERY_COVERAGE.json` for exact review/result identities per source.', '', '## Typed connections and time', '', 'Delivery proposes evidence-to-K bindings, no new physical edges or native statuses. Each K row carries resolution and timescale in results.json. Examples of future underlying ports: antagonist triangles -> crown roof PER_POINT/SIMULTANEOUS; regional signed gap -> seating PER_POINT/SIMULTANEOUS; local final drilling state -> later biology PER_POINT/HANDOVER; process -> as-built PER_POINT/HANDOVER. Coordinator status proposal is administrative and must not become a physical coupling.', '', 'No working node covers the full Northstar K table. DENT-GOAL is absent from working packet; DENT-VAL-BASELINE-COMPARISON is narrower and must not carry a whole-service claim. Coverage proposal and pending feedback in HANDOFF.md/GRAPH_FEEDBACK.json. Source graphs and NORTHSTAR_SPEC.md remain with the coordinator.']
    (HERE / 'NORTHSTAR_STATUS_PROPOSAL.md').write_text('\n'.join(summary) + '\n')
    (HERE / 'RESULTS.md').write_text(f"""# Reviewed deliveries can now be mapped back to all 52 Northstar chains

{len(result['status_changes'])} proposed status changes: """ + ', '.join((x['k'] for x in result['status_changes'])) + f""". Counted proposal: """ + ' · '.join((f'{s} {new[s]}' for s in CATEGORIES)) + '. No new PRESENT promotion.\n\nSee [full status proposal](NORTHSTAR_STATUS_PROPOSAL.md), [machine-readable results](results.json) and [figure](status_counts.svg). Reference consists of other lanes reviewed, hash-matched scoped deliveries. No new physical measurement.\n\n96/96 review bindings match, 67 unique result identities. Seven bad bookkeeping inputs rejected; manual semantic mapping remains PENDING_INDEPENDENT_REVIEW.\n\nLimitation: PRESENT inherited for seven scoped components without a new complete review. Full patient/manufacturing chain unvalidated. K33 refuted component retained. Old demo REJECT and negative physical gates not rewritten.\n')
    (HERE / 'README_DEMO.md').write_text('# Northstar: what the reviews actually enable\n\nA coordinator can now revise each K status with exact contract, review class, principal value and locator. New connection links completed reviews to the older status map.\n\nRun **`./run_all.sh`** here. Python 3, one thread, standard library, no network/GPU. Command checks frozen inputs and rewrites proposal, results.json, figure and validation. It runs no old FE/dataset jobs.\n\n| Status | Before | Proposal |\n|---|---:|---:|\n' + ''.join((f'| {s} | {old[s]} | {new[s]} |\n' for s in CATEGORIES)) + '\n![Status counts](status_counts.svg)\n\nWhat fails: digital surfaces, published group metrics and honest negative reporting are not physical manufacturing validation. No new chain becomes PRESENT. Manual status mapping awaits independent coordinator review. Seven controls reject wrong rows/locators/promotions/hashes; they do not certify semantics.\n\nReference/license: local private review/specification documents with exact hashes in SOURCE_MANIFEST.json. X43 uses no raw imaging/geometry datasets and relicenses no sources. Underlying data retain demo-specific licenses: Bits2Bites stated CC BY-NC-SA, ToothFairy2 CC BY-SA, STS and mandibular defects CC BY 4.0; exact Bite2Text/Pulpy permissions and Teeth3DS snapshot restricted/UNKNOWN according to reviews. No such raw meshes copied/exported here.\n\nSee [NORTHSTAR_STATUS_PROPOSAL.md](NORTHSTAR_STATUS_PROPOSAL.md) for all 52 rows and concrete data/code/measurement requirements for remaining MISSING chains.\n')
    write(HERE / 'CURRENT_WORK_STATE.json', {'lane': 'X43-northstar-status', 'claim_type': 'information_link', 'stage': 'PROPOSAL_VERIFIED_PENDING_COORDINATOR', 'updated_at': stamp, 'latest_gate': 'PASS: 52 rows, exact counts, source/review hashes, seven rejecting fault probes; no unsupported PRESENT', 'next_operation': 'Independent coordinator semantic review and scoped merge into NORTHSTAR_SPEC.md; then execute K48 matched crown scan experiment', 'status_changes': result['status_changes'], 'proposed_counts': new, 'canonical_spec_unchanged': True, 'large_arrays': False})
    print(json.dumps({'gate': 'PASS', 'rows': len(rows), 'old_counts': old, 'proposed_counts': new, 'changes': result['status_changes'], 'fault_probes_rejected': len(mutations), 'spec_unchanged': True}, ensure_ascii=False, indent=2))
if __name__ == '__main__':
    main()
