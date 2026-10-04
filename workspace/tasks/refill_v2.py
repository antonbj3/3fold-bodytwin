"""Graph-guided The swarm refill. Argument 0 previews the next refill without writes."""
from __future__ import annotations

import glob
import fcntl
import hashlib
import json
import os
import random
import re
import shutil
import sys
import time
from pathlib import Path

def _private_input_pattern(public_pattern):
    """Combine public exclusions with required operator-maintained private terms."""
    import json
    from pathlib import Path
    import re
    config = Path.home() / ".bodytwin" / "private_source_terms.json"
    terms = json.loads(config.read_text())
    if not isinstance(terms, list) or not terms or any(not isinstance(t, str) or not t.strip() for t in terms):
        raise ValueError("A nonempty private-source exclusion list is required")
    return "(?:" + public_pattern + ")|(?:" + "|".join(re.escape(t) for t in terms) + ")"


W = Path(__file__).resolve().parents[1]
R = W / 'results'
MAP = R / 'CX-QGRAPH/Q_GRAPH_MAP.json'
PAIRS = R / 'CX-QGRAPH/PAIRS.json'
# Manual triage: a graph path is a lead, not evidence of shared measured data.
# The excluded mitochondria cluster and two-hop energy paths repeatedly generated
# combinations with no common subject, observable, or clock.
ACTIONABLE_PAIRS = {
    ('Q009', 'Q018'), ('Q012', 'Q115'), ('Q029', 'Q127'),
    ('Q037', 'Q038'), ('Q043', 'Q090'), ('Q115', 'Q148'),
    ('Q148', 'Q168'),
}
VALIDATION_TARGETS = {
    'Q005': 'renal filtration and creatinine: independent held-out patient stays, time-aligned laboratory and urine observations',
    'Q012': 'TMA/TMAO: independent cohort with measured precursor and product, never infer conversion fraction from a single cross-sectional ratio',
    'Q013': 'local versus plasma drug PK: independently recalculate molar units (1 ng/mL = 1 microgram/L; at MW 528.9 g/mol, 1 ng/mL = 1/528.9 micromol/L) and recheck all AUC and ratio claims against the raw matched subjects',
    'Q052': 'bone porosity and permeability: paired geometry and flow from the same specimen with explicit voxel length unit',
    'Q146': 'nasal nicotine PK: same formulation and dose with time-resolved concentration; keep NHANES serum cotinine as an exposure context, not nasal-route validation',
}
# These findings require fresh source checks, not a silent edit of a producer's result.
QUALITY_TARGETS = {
    'Q005': 'patient.csv has 29 observed columns, while the report says 31; reconcile the downloaded bytes, parser and row/column counts',
    'Q052': 'NIfTI xyzt_units=0; establish voxel length from an independent source before any metric geometry comparison',
    'Q146': 'population serum cotinine is exposure context, not nasal-route PK; locate dose-matched time-series measurements',
    'Q013': 'audit the DATX Q012-Q013 molar conversion from raw source values: 1 ng/mL = 1/528.9 micromol/L at MW 528.9 g/mol',
}
# Each item is a named gap in the completed dental packets, not an inferred
# graph dependency. Keep the existing negative findings in the new context.
DENTAL_FOLLOWUPS = {
    'LOAD-CONTACT': ('BT-DENT-LOADIF', 'Test the lever-only tooth-load distribution against a public simultaneous multi-tooth contact map; distinguish capacity ceiling from service load.'),
    'LOAD-ANATOMY': ('BT-DENT-LOADIF', 'Replace assumed condyle distance, PDL support area and tooth inclination with measured, source-traced anatomy; propagate uncertainty.'),
    'KIN-ARCH': ('BT-DENT-DATX-KINxGEOM', 'Find a registered mandibular arch transform and opposing dentition for the same subject; freeze frame and contact definitions before testing.'),
    'PDL-COMP': ('BT-DENT-DATX-PDLxLOAD', 'Find a public PDL compression measurement with matched load mode and units; test the existing tension-derived parameter as a baseline.'),
    'TMJ-FORCE': ('BT-DENT-DATX-TMJxLOAD', 'Seek direct joint-force measurements or a clearly labelled validated force estimate and test the existing load-ratio interval.'),
    'CERAM-AREA': ('BT-DENT-DATX-LOADxFAT', 'Measure effective regional contact area and moisture-dependent ceramic fatigue parameters before updating pf5.'),
}
PRIVATE_MARKERS = re.compile(_private_input_pattern('(?<![A-Za-z0-9])(collaborator|collaborator|gc)(?![A-Za-z0-9])'), re.I)
NON_DATA_PREFIXES = ('AUTO-', 'AUDIT-', 'BUILD-', 'FIX-', 'GRAPH-', 'SOLVE-')
HEAD = ("Read inputs/NIGHT_PREAMBLE.md. Just work here, ≤ 45 min, 1 thread. Web allowed for published literature and public datasets. PREREG.md + sha256 before calculation; RESULTS.md begins with the line \"{j}\"; results.json. No judgment words, do not find on the measurement data.\n\n")


def done(d: Path) -> bool:
    return (d / 'RESULTS.md').is_file()


def criterion_met(d: Path) -> bool:
    """Conservative explicit criterion check, with UNKNOWN and failure vetoes."""
    if not done(d):
        return False
    text = (d / 'RESULTS.md').read_text(errors='replace').lower()
    if re.search("\\b(ej uppfyllt|not fulfilled|not met|criterion fail|kriterium: unknown)\\b", text):
        return False
    try:
        data = json.loads((d / 'results.json').read_text())
    except (OSError, ValueError):
        data = {}
    for key in ('frozen_criteria', 'criteria', 'criterion', 'criterion_check'):
        part = data.get(key)
        if isinstance(part, dict):
            for flag in ('all_pass', 'all_passed', 'all_met', 'criterion_met'):
                if isinstance(part.get(flag), bool):
                    return part[flag]
    for key in ('criterion_met', 'all_passed'):
        if isinstance(data.get(key), bool):
            return data[key]
    status = data.get('status')
    if isinstance(status, str):
        if status.lower() in ('pass', 'passed', 'criterion_met'):
            return True
        if any(x in status.lower() for x in ('fail', 'unknown', 'not met')):
            return False
    return bool(re.search(r'\b(kriterium uppfyllt|criterion met|all criteria passed)\b', text))


def model_files(q: str):
    d = R / f'BT-HX-{q}'
    return [(d / f, f'{q}_{f}') for f in ('model.py', 'PREREG.md', 'RESULTS.md', 'results.json')] + [
        (d / 'inputs/QUESTION.md', f'{q}_QUESTION.md')]


def load_graph():
    mapping = json.loads(MAP.read_text())['questions']
    pairs = json.loads(PAIRS.read_text())['pairs']
    score = {}
    priority = {a['id']: a.get('priority', 0) or 0 for a in json.loads((W / 'PRIORITY.json').read_text())['actions']}
    next_score = {a['id']: a.get('score', 0) or 0 for a in json.loads((W / 'NEXT_ACTIONS.json').read_text())['actions']}
    for item in mapping:
        q = item['question_id'].replace('BT-HX-', '')
        nodes = [m['node_id'] for m in item['matches']]
        score[q] = max([float(priority.get(n, 0)) + .01 * float(next_score.get(n, 0)) for n in nodes] + [0])
    return score, pairs


def robust_models():
    """Read the R1 table only; R2 and prose mentions never qualify."""
    summary = R / 'BT-HX-SYNTH-2021/SUMMARY.md'
    if not summary.is_file():
        return []
    section = summary.read_text(errors='replace').split('### R1 — Robusta (alla tre villkor uppfyllda)', 1)
    if len(section) < 2:
        return []
    section = section[1].split('### R2', 1)[0]
    return sorted(set(re.findall(r'BT-HX-(Q\d{3})', section)))


def graph_data_candidates():
    """Unmapped, actionable empirical graph nodes; source priority is advisory."""
    mapped = {m['node_id'] for q in json.loads(MAP.read_text())['questions'] for m in q['matches']}
    priorities = {a['id']: a for a in json.loads((W / 'PRIORITY.json').read_text())['actions']}
    actions = json.loads((W / 'NEXT_ACTIONS.json').read_text())['actions']
    candidates = []
    for action in actions:
        node = action['id']
        p = priorities.get(node, {})
        if (node in mapped or not action.get('actionable') or not p.get('actionable')
                or action.get('status') not in ('OPEN', 'ASSUMED')
                or action.get('type') not in ('EMPIRICAL', 'DATA-ANCHOR')
                or node.split(':', 1)[-1].startswith(NON_DATA_PREFIXES)
                or PRIVATE_MARKERS.search(node + ' ' + str(action.get('claim', '')))):
            continue
        candidates.append((node, float(p.get('priority') or 0), float(action.get('score') or 0)))
    return sorted(candidates, key=lambda x: (-x[1], -x[2], x[0]))


def graph_data_id(node):
    return 'BT-DATG-' + hashlib.sha256(node.encode()).hexdigest()[:16].upper()


def graph_data_brief(node, priority, next_score):
    return (f'Nodstyrd PUBLIK data search for graphnode {node}. Grafprioritet {priority:g}, next-action-score {next_score:g} are only sample signals. population/regim, unit of measurement, time base and pre-search decision criteria. Find an actual independent raw data sample; specify URL/DOI, license, filhash, parsad form, units and measurement uncertainty in DATA_SOURCES.json. Compare with the exact question in the package. Read no private graphs from the cloud.ID: not unambiguously indicate observable and intended comparison: account NEEDS_TASK_DEFINITION and the missing fields before wide search; do not guess internal meaning. This is a task gap, not a negative research result. Do not copy internal graph text, the collaborator/GC/restricted model data or private raw data to a ALLOW_WEB-paket or any web service. If no matching source is available: UNKNOWN and an explicit coverage gap. Upgrade no graph status.')


def integrated_models():
    return sorted(d.name.removeprefix('BT-INT-') for d in R.glob('BT-INT-Q*') if done(d))


def plan(limit: int = 80, now: float | None = None):
    now = time.time() if now is None else now
    score, pairs = load_graph()
    selected = []
    finished = [d for pattern in ('BT-DAT-Q*', 'BT-DATG-*', 'BT-DATX-Q*', 'BT-INT-Q*')
                for d in R.glob(pattern) if d.is_dir() and done(d) and (d/'ALLOW_WEB').exists()]
    auditable=sorted((d for d in finished if not (R/('BT-AUD-V2-'+d.name)).exists()),key=lambda d:(d/'RESULTS.md').stat().st_mtime,reverse=True)
    audit_budget=min(len(auditable), max(1,limit//10)) if limit else 0
    def add(j, brief, files=()):
        cap = limit if j.startswith('BT-AUD-') else max(0, limit - audit_budget)
        if len(selected) < cap and not (R / j).exists() and j not in {x[0] for x in selected}:
            selected.append((j, brief, list(files)))
            return True
        return False

    qs = [f'Q{i:03}' for i in range(1, 169) if done(R / f'BT-HX-Q{i:03}')]
    for q in sorted(qs, key=lambda x: (-score.get(x, 0), x)):
        if len(selected) >= int(limit * .45): break
        add(f'BT-DAT-{q}', "Find PUBLIC measured datasets constraining the model in inputs/. For each: URL, licence, size, format, subjects, parameter/output. Download a sample ≤ 50 MB, check the actual parsed shape, source-based units and coordinate frame; save in samples/ and DATA_SOURCES.json with sha256. A missing unit must not be replaced with an assumption in a measurement comparison. An interval with only a measured lower bound must not PASS against a two-sided band. Distinguish measured value, derivation and model assumption.", model_files(q))

    for pair in pairs:
        if len(selected) >= int(limit * .65): break
        qa, qb = [x.replace('BT-HX-', '') for x in pair['questions']]
        if (qa, qb) not in ACTIONABLE_PAIRS:
            continue
        a, b = R / f'BT-DAT-{qa}', R / f'BT-DAT-{qb}'
        if not all(done(d) and (d / 'DATA_SOURCES.json').is_file() for d in (a, b)):
            continue
        files = [(a / 'DATA_SOURCES.json', f'{qa}_DATA_SOURCES.json'), (b / 'DATA_SOURCES.json', f'{qb}_DATA_SOURCES.json')]
        detail = json.dumps(pair['evidence'][0], ensure_ascii=False)
        add(f'BT-DATX-{qa}-{qb}', f'KOMBINERA dataset. Graflink: {detail}. Load yourself and parse at least one actual data sample per source; check sha256, form, units, coordinate frame and time base. Write an integral SI-conversion for each quantity of measurement: for molecular mass 528,9 g/mol applies to: 1 ng/mL = 1/528,9 µmol/L, not 1/0,5289. Visa explicit om individ-ID, states and observables can be matched. If there is no common measurement: report unmatched and no parameter improvement, without synthetic joy. Calculate uncertainty before/after only for identifiable quantity. A coupling annotation is hypothesis, not evidence.', files + model_files(qa) + model_files(qb))

    for q, finding in QUALITY_TARGETS.items():
        if len(selected) >= int(limit * .65): break
        dat = R / f'BT-DAT-{q}'
        if not done(dat): continue
        files = [(dat / 'RESULTS.md', f'{q}_DAT_RESULTS.md'),
                 (dat / 'DATA_SOURCES.json', f'{q}_DATA_SOURCES.json'),
                 (R / 'CX-SWARMPLANNER/QUALITY_LOG.md', 'QUALITY_LOG.md')]
        add(f'BT-DATQ-{q}', f'Oberoende datakvalitetsomtag: {finding}. Reopen actual source file via samples/ or SAMPLES_MOVED.txt; kontrollera SHA256, parses, shape, units and measurement definition with the sharpened DAT-mallen. Frys det specifika felet i PREREG. Rapportera korrigerat eller UNKNOWN, with an impact on previous conclusion. collaborator/GC/restricted model data shall not be included.', files)

    for q in robust_models():
        if len(selected) >= int(limit * .85): break
        d = R / f'BT-HX-{q}'
        files = model_files(q) + [(R / f'BT-HX-{q}-R/RESULTS.md', f'{q}_R_RESULTS.md'),
                                  (W / 'bodytwin_core/__init__.py', 'bodytwin_core__init__.py'),
                                  (W / 'bodytwin_core/contact_band_fast.py', 'TEMPLATE_module.py')]
        if (R / f'BT-DAT-{q}/DATA_SOURCES.json').is_file():
            files.append((R / f'BT-DAT-{q}/DATA_SOURCES.json', f'{q}_DATA_SOURCES.json'))
        if done(d):
            add(f'BT-INT-{q}', "Package the model robust in the synthesis packet’s R1 list as a bodytwin_core module. Common interface with units, source docstring and pytest for limiting cases and verified reference. Preserve the synthesis corrections and limitations.", files)

    for q in integrated_models():
        if len(selected) >= int(limit * .85): break
        dat, integrated = R / f'BT-DAT-{q}', R / f'BT-INT-{q}'
        if not done(dat): continue
        files = [(integrated / 'RESULTS.md', f'{q}_INT_RESULTS.md'),
                 (dat / 'RESULTS.md', f'{q}_DAT_RESULTS.md'),
                 (dat / 'DATA_SOURCES.json', f'{q}_DATA_SOURCES.json')]
        add(f'BT-VAL2-{q}', "Validate the integrated module against a SECOND independent public dataset, not the source/individuals/endpoints used for fitting or the first DAT packet. Freeze source separation, baseline, units, population/regime and outcome in PREREG. Open and parse the sample, check licence and sha256 and run the module’s public interface. Report a negative finding or UNKNOWN if independent matching is missing; simulation is not external validation. Internal the collaborator/GC/restricted model data must not be included.", files)

    # Validation after the data and integration passes: these known gaps
    # have concrete measured endpoints, unlike another broad ontology pair.
    for q, target in VALIDATION_TARGETS.items():
        dat = R / f'BT-DAT-{q}'
        if not done(dat):
            continue
        files = model_files(q) + [(dat / 'RESULTS.md', f'{q}_DAT_RESULTS.md')]
        if (dat / 'DATA_SOURCES.json').is_file():
            files.append((dat / 'DATA_SOURCES.json', f'{q}_DATA_SOURCES.json'))
        add(f'BT-VAL-{q}', f'Independently measured validation of {target}. Freezer outputs, units of measurement, baseline and full-out separation in PREREG before data analysis. Read and check an actual public sample; enter DOI/URL, licens, sha256, number of individuals, sample format and uncertainty. Do not reuse the same individuals or endpoint for both fit and test. If the matched open sample is missing, the outcome is UNKNOWN; if the model loses to baseline, it is a valid negative result. collaborator/GC/restricted model data may not be loaded or described here.', files)

    for pair in pairs:
        if len(selected) >= limit: break
        qa, qb = [x.replace('BT-HX-', '') for x in pair['questions']]
        if (qa, qb) not in ACTIONABLE_PAIRS:
            continue
        if not (criterion_met(R / f'BT-HX-{qa}') and criterion_met(R / f'BT-HX-{qb}')):
            continue
        detail = json.dumps(pair['evidence'][0], ensure_ascii=False)
        add(f'BT-CPL-{qa}-{qb}', f'Koppla modellerna via den grafrelaterade storheten. Relation: {detail}. Specify the port, units and frozen test against a measured value that no model alone predicts. Validate the physics of the relationship before any dependence.', model_files(qa) + model_files(qb))

    for name, (parent, objective) in DENTAL_FOLLOWUPS.items():
        if len(selected) >= limit - audit_budget: break
        d = R / parent
        if done(d):
            add(f'BT-DENT2-{name}', f'Continue From {parent}: {objective} Preserve past negative results and uncertainty. individ/regim, SI-units and a frozen counterexample; no synthetic joy or graph status change. Use only public data in this ALLOW_WEB-paket.', [(d / 'RESULTS.md', parent + '_RESULTS.md'), (d / 'results.json', parent + '_results.json')])

    # Fill spare capacity only after the finite validation and coupling work.
    # This also keeps refills useful after those finite families are exhausted.
    for node, priority, next_score in graph_data_candidates():
        if len(selected) >= limit - audit_budget: break
        add(graph_data_id(node), graph_data_brief(node, priority, next_score))

    # One focused, stable audit per target; include the dominant DATG stream.
    for d in auditable[:audit_budget]:
        files=[(d/'RESULTS.md',d.name+'_RESULTS.md')]
        files += [(d/f,d.name+'_'+f) for f in ('results.json','DATA_SOURCES.json','SAMPLES_MOVED.txt','BRIEF.md') if (d/f).is_file()]
        sample=next((f for f in sorted((d/'samples').rglob('*')) if f.is_file() and f.stat().st_size<20e6),None)
        if sample:files.append((sample,d.name+'_'+sample.name))
        add('BT-AUD-V2-'+d.name, 'Oberoende granskning av '+d.name+". Check exact task definition, primary source, actual parsing, units, sample shape and independence between fit/test. Choose the conclusion that most affects the next model decision and recompute it. Report separately: verified finding, model assumption, missing material or wrong task definition. An honest UNKNOWN can be correct but does not mean increased model capability. State concretely which next experiment or implementation decision the result supports. No graph status may be upgraded.",files)
    selected.sort(key=lambda item: 0 if item[0].startswith('BT-AUD-') else 1)
    return selected


def write_packets(selected):
    for j, brief, files in selected:
        d = R / j
        storage=Path('external_mount')
        if storage.is_dir() and not d.exists():
            dest=storage/j;dest.mkdir(exist_ok=False);d.symlink_to(dest,target_is_directory=True)
        (d / 'inputs').mkdir(parents=True)
        shutil.copy2(W / 'tasks/NIGHT_PREAMBLE.md', d / 'inputs/NIGHT_PREAMBLE.md')
        (d / 'ALLOW_WEB').touch()
        for src, name in files:
            if src.is_file() and src.stat().st_size < 50e6:
                shutil.copy2(src, d / 'inputs' / name)
        (d / 'BRIEF.md').write_text(f'# {j}\n' + HEAD.format(j=j) + brief + '\n')
    if selected:
        with (W / 'tasks/lanes/bt_queue.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            with (W / 'tasks/lanes/bt_queue.txt').open('a') as f:
                f.write(''.join(f"{'AB'[i % 2]} swarm {j}\n" for i, (j, _, _) in enumerate(selected)))


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    requested = int(args[0]) if args else 80
    if requested < 0:
        raise ValueError('limit must be nonnegative')
    preview = requested == 0
    selected = plan(80 if preview else requested)
    if not preview:
        write_packets(selected)
    print(('WOULD create' if preview else 'Created') + f' {len(selected)} packets:')
    for j, _, _ in selected:
        print(j)


if __name__ == '__main__':
    main()
