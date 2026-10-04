'Item-complete reporting ledgers with explicit scope and conservative evidence status.\n\nShort English labels paraphrase official items. No checklist score is a quality score.\nCRIS 2014 has NO official numbered checklist: local thematic prompts are labelled as such.\n'
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parent
CLAIM = 'AI identity|Summary|Usage background|Hypotheses|Temporal design|Study task|Data sources|Selection rules|Preprocessing|Subset selection|Deidentification|Missing data|Imaging protocol|Reference method|Reference rationale|Annotation source|Test annotation|Rater variation|Data partitions|Independent split|Sample rationale|Model description|Software libraries|Parameter initialization|Training method|Model selection|Ensemble|Performance metrics|Statistical uncertainty|Robustness|Explainability|Internal test|External test|Trial registration|Participant flow|Cohort characteristics|Result uncertainty|Diagnostic accuracy|Error analysis|Limitations|Usage implications|Protocol access|Resource access|Funding'.split('|')
TRIPOD_IDS = '1 2 3a 3b 3c 4 5a 5b 6a 6b 6c 7 8a 8b 8c 9a 9b 9c 10 11 12a 12b 12c 12d 12e 12f 12g 13 14 15 16 17 18a 18b 18c 18d 18e 18f 19 20a 20b 20c 21 22 23a 23b 24 25 26 27a 27b 27c'.split()
TRIPOD = 'Study identity|Abstract|Care context|Target users|Inequality|Purpose|Data sources|Collection dates|Study centers|Inclusion|Treatments|Data preparation|Outcomes|Outcome assessors|Outcome blinding|Predictor selection|Predictor measurement|Predictor assessors|Sample rationale|Data attrition|Data division|Predictor handling|Model building|Cluster heterogeneity|Performance metrics|Recalibration|Prediction computation|Class imbalance|Fairness|Model output|Cohort differences|Ethics|Funding|Conflicts of interest|Protocol|Registration|Data access|Code access|Patient involvement|Participant flow|Cohort characteristics|Cohort comparison|Analysis denominator|Model specification|Result intervals|Cluster results|Update results|Interpretation|Limitations|Input errors|User requirements|Further research'.split('|')
STARD_IDS = '1 2 3 4 5 6 7 8 9 10a 10b 11 12a 12b 13a 13b 14 15 16 17 18 19 20 21a 21b 22 23 24 25 26 27 28 29 30'.split()
STARD = 'Study identity|Summary|Test usage|Hypotheses|Temporal design|Inclusion|Recruitment basis|Recruitment setting|Sample type|Index test|Reference test|Reference rationale|Test thresholds|Reference thresholds|Test blinding|Reference blinding|Statistics|Indeterminate answers|Data attrition|Variation analysis|Sample rationale|Flowchart|Cohort characteristics|Disease grades|Alternative diagnoses|Test timing relation|Cross-tabulation|Precision intervals|Adverse effects|Limitations|Usage implications|Registration|Protocol|Funding'.split('|')
CRIS = 'Sample rationale|Meaningful effect|Specimen preparation|Provenance clusters|Losses|Random allocation|Blinding|Primary/secondary statistics'.split('|')
SOURCES = {'CLAIM_2024': dict(doi='10.1148/ryai.240300', url='https://pubs.rsna.org/doi/pdf/10.1148/ryai.240300', items=44, verified=True, locator='Table, printed page3, items1–44', note='Imaging AI only; not a quality score. Non-AI geometry uses explicitly adapted questions.'), 'TRIPOD_AI': dict(doi='10.1136/bmj-2023-078378', url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11019967/#tbl2', items=52, verified=True, locator='Table2:27 main items,52 subitems', note='Clinical prediction models. Engineering material forecasts and generative CAD fall outside formal scope.'), 'STARD_2015': dict(doi='10.1136/bmj.h5527', url='https://www.equator-network.org/wp-content/uploads/2015/03/STARD-2015-checklist.pdf', items=34, verified=True, locator='2015 checklist,30 main items plus a/b subitems', note='Diagnostic accuracy; geometric category and synthetic comparators are scoped adaptations. 10.1136/bmjopen-2016-012799 is explanation/elaboration, not the2015 statement.'), 'CRIS_CONCEPTS': dict(doi='10.4103/0972-0707.136338', url='https://pmc.ncbi.nlm.nih.gov/articles/PMC4127685/', items=8, verified=True, locator='CURRENT LACUNA IN REPORTING IN VITRO STUDIES; Need for CRIS Guidelines', note='2014 concept note has no numbered finalized checklist. Eight LOCAL thematic prompts, not official CRIS items or certification.')}
CLAIM_DEMOS = {'PROOF_LANE', 'GENCAD_V2', 'GENCAD_V3', 'X3', 'X5', 'X9', 'X12', 'X15', 'X20', 'X22', 'X24', 'X30'}
TRIPOD_DEMOS = {'X7', 'X21', 'X1B', 'X13', 'X14', 'X16', 'X17', 'X19', 'X25', 'X27', 'X29'}
STARD_DEMOS = {'X7', 'X21', 'X22'}
CRIS_DEMOS = {'PROOF_LANE', 'X1B', 'X2', 'X10', 'X13', 'X14', 'X16', 'X19', 'X23', 'X25', 'X27', 'X28', 'X29', 'X33'}

def items(guideline):
    assert len(CLAIM) == 44 and len(TRIPOD_IDS) == len(TRIPOD) == 52 and (len(STARD_IDS) == len(STARD) == 34)
    if guideline == 'CLAIM_2024':
        return list(zip(map(str, range(1, 45)), CLAIM))
    if guideline == 'TRIPOD_AI':
        return list(zip(TRIPOD_IDS, TRIPOD))
    if guideline == 'STARD_2015':
        return list(zip(STARD_IDS, STARD))
    return list(zip(['LOCAL-' + str(i) for i in range(1, 9)], CRIS))

def status(id, guideline, number, label, demo):
    """Every YES/PARTIAL has scoped artifact evidence; absent reports stay not located.

 This is a demo-artifact review, not an audit of each source study's investigators.
 """
    p = Path(demo['path'])
    readme = p / 'README_DEMO.md'
    files = {x.name for x in p.iterdir() if x.is_file()}
    common = 'NOT_LOCATED'
    reason = 'Not established item by item in reviewed demo artifacts; additional reporting or the primary study owner is required.'
    evidence = []
    good = {'Hypotheses', 'Purpose', 'Study task', 'Usage background', 'Target users', 'Reference method', 'Reference rationale', 'Reference test', 'Index test', 'Model description', 'Model output', 'Limitations', 'Interpretation', 'Usage implications', 'Further research', 'User requirements'}
    if label in good:
        common = 'PARTIAL'
        reason = 'README describes the demo object; the complete guideline item, including clinical/primary empirical requirements, is not established.'
        evidence = [str(readme)]
    if label in {'Data sources', 'Data access', 'Resource access', 'Code access', 'Prediction computation', 'Model specification', 'Protocol', 'Protocol access'}:
        common = 'PARTIAL'
        reason = 'Local sources/code/protocols exist; private machine-specific access, licensing and export conditions must follow the source.'
        evidence = [str(readme)] + [str(p / f) for f in sorted(files) if f.startswith('PREREG') or f == 'run_all.sh']
    if label in {'Participant flow', 'Data attrition', 'Indeterminate answers', 'Analysis denominator', 'Provenance clusters', 'Sample rationale'}:
        common = 'PARTIAL'
        reason = 'Counts and attrition are available, but representative independent sampling, a planned relevant effect and/or raw specimen lineage are incomplete.'
        evidence = [str(readme), str(ROOT / 'STATISTICS.md') + '#' + id]
    if label in {'Statistical uncertainty', 'Result uncertainty', 'Precision intervals', 'Result intervals', 'Cluster results', 'Primary/secondary statistics'}:
        have = any((m.get('cluster_ci95') for m in demo['metrics']))
        common = 'PARTIAL' if have else 'MISSING'
        reason = 'X35 conditional case/study intervals; no refit/measurement/transport uncertainty.' if have else 'Missing: raw independent specimen lineage/variance, or statistics do not apply to a deterministic value.'
        evidence = [str(ROOT / 'STATISTICS.md') + '#' + id, str(ROOT / 'results.json')]
    if label in {'Sample type', 'Temporal design', 'Test timing relation', 'Data division', 'Independent split', 'Data partitions'}:
        common = 'PARTIAL'
        reason = 'Retrospective local analysis. Split contracts are inspected as case IDs, not guaranteed population/center independence.'
        evidence = [str(readme)]
    if label in {'External test', 'Trial registration', 'Registration'}:
        common = 'NOT_APPLICABLE' if label == 'Trial registration' else 'MISSING'
        reason = 'No new clinical trial.' if common == 'NOT_APPLICABLE' else 'Missing: external untouched test cohort for the same prediction target; external paper/reference is not external testing.' if label == 'External test' else 'Missing: registration number. Local PREREG is not registry registration.'
        evidence = [str(readme)]
    if label in {'Disease grades', 'Alternative diagnoses', 'Adverse effects', 'Treatments', 'Patient involvement'}:
        common = 'NOT_APPLICABLE'
        reason = 'Scoped retrospective geometry/material/synthetic study without new treatment or patient recruitment; clinical endpoints absent.'
        evidence = [str(readme)]
    if label in {'Meaningful effect', 'Random allocation', 'Blinding', 'Outcome blinding', 'Test blinding', 'Reference blinding', 'Specimen preparation', 'Rater variation', 'Outcome assessors', 'Predictor assessors', 'Collection dates', 'Funding', 'Conflicts of interest', 'Ethics', 'Imaging protocol', 'Deidentification', 'Recruitment setting', 'Cohort characteristics', 'Cohort comparison', 'Inequality', 'Fairness'}:
        common = 'NOT_LOCATED'
        reason = 'Missing in this demo reanalysis or requires the original data/study owner. Not inferred from README keywords.'
        evidence = []
    if label in {'Training method', 'Parameter initialization', 'Model selection', 'Ensemble', 'Model building', 'Class imbalance', 'Predictor handling', 'Predictor selection', 'Recalibration', 'Update results'} and id not in {'X7', 'X21', 'X22', 'X13', 'X14', 'X25'}:
        common = 'NOT_APPLICABLE'
        reason = 'No corresponding trained AI model in the selected main query; deterministic geometry or literature reanalysis.'
        evidence = [str(readme)]
    evidence = [e + ':1' if e == str(readme) else e for e in evidence]
    return dict(item=number, label=label, status=common, reason=reason, evidence=evidence)

def build(demos):
    out = {}
    dest = ROOT / 'CHECKLISTS'
    dest.mkdir(exist_ok=True)
    for demo in demos:
        id = demo['id']
        gs = []
        if id in CLAIM_DEMOS:
            gs.append('CLAIM_2024')
        if id in TRIPOD_DEMOS:
            gs.append('TRIPOD_AI')
        if id in STARD_DEMOS:
            gs.append('STARD_2015')
        if id in CRIS_DEMOS:
            gs.append('CRIS_CONCEPTS')
        matrices = []
        lines = [f'# {id}: reporting matrix', '', f"Scope: demo package inspected at {demo['path']}. Status concerns local artifacts, not original study quality. PRESENT means reported, not true. NOT_LOCATED means absent from the reviewed view.", '']
        if not gs:
            lines += ['The four requested guidelines do not formally apply to this deterministic design/geometry query. STATISTICS.md records provenance, denominators and open empirical requirements. No empty checklist counts as PASS.', '']
        for g in gs:
            scope = 'DIRECT_IMAGING_AI_SYNTHETIC' if g == 'CLAIM_2024' and id == 'X22' else 'LOCAL_CONCEPT_AUDIT' if g == 'CRIS_CONCEPTS' else 'ADAPTED_SCOPE_NOT_FORMAL_COMPLIANCE'
            rows = [status(id, g, n, l, demo) for (n, l) in items(g)]
            matrices.append(dict(guideline=g, scope=scope, source=SOURCES[g], rows=rows))
            lines += [f'## {g}', '', f"Source: [{SOURCES[g]['doi']}]({SOURCES[g]['url']}); {SOURCES[g]['locator']}. {SOURCES[g]['note']}", f'Application: {scope}.', '', '| Item | Short label | Status | Reason and locator |', '|---|---|---|---|']
            for r in rows:
                refs = '; '.join((f"[{Path(e.split('#')[0]).name}]({e})" for e in r['evidence']))
                lines.append(f"| {r['item']} | {r['label']} | {r['status']} | {r['reason']} {refs} |")
            lines += ['']
        (dest / (id + '.md')).write_text('\n'.join(lines) + '\n')
        out[id] = dict(demo=id, checklists=matrices, status='INCOMPLETE_REPORTING' if gs else 'NOT_APPLICABLE_SCOPE', count_is_not_quality_score=True)
    (dest / 'checklists.json').write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n')
    (dest / 'SOURCES.json').write_text(json.dumps(SOURCES, ensure_ascii=False, indent=2) + '\n')
    return out
