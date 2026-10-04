"""Freeze the executable union and create a small, clean writable package instance."""
from dental_release.paths import expand as _release_expand
import datetime, hashlib, json, os, shutil
from pathlib import Path
HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent / 'DEMO48_PACKAGE'
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X86_sunday_refresh'))
ALIAS = DATA / 'readonly_package'

def execution_manifest():
    revision = HERE / 'FROZEN_EXECUTIONS_R3.json'
    return revision if revision.exists() else HERE / 'FROZEN_EXECUTIONS.json'

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda : f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def write(p, value):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')

def main():
    DATA.mkdir(parents=True, exist_ok=True)
    ALIAS.mkdir(exist_ok=True)
    catalog = json.loads((PACKAGE / 'demos.json').read_text())['demos']
    rows = []
    for r in catalog:
        if r['id'] == 'GENCAD_V2':
            continue
        batch = r['copy_dir'].split('/')[0]
        if batch == 'demos':
            batch = 'base'
        rows.append({'id': r['id'], 'batch': batch, 'route': ['run_demos.sh', r['id']], 'catalog': r, 'scope': r.get('execution_scope', 'Read bundled runner and independent review')})
    for (batch, key) in [('batch3', 'demos'), ('batch4', 'demos'), ('batch5', 'demos'), ('batch6', 'demos'), ('batch7', 'executables'), ('batch9', 'executables')]:
        for r in json.loads((PACKAGE / batch / 'INDEX.json').read_text())[key]:
            if r['id'] not in {q['id'] for q in rows}:
                rows.append({'id': r['id'], 'batch': batch, 'route': ['/usr/bin/python3', batch + '_runner.py', r['id']], 'catalog': r, 'scope': r.get('scope', 'Reviewed batch pipeline with frozen controls')})
    assert len(rows) == 93 and len({r['id'] for r in rows}) == 93
    priority = ['X82', 'NET_R4', 'X84', 'X85', 'PROOF_LANE_FULL_CROWN_R3']
    rows.sort(key=lambda r: (priority.index(r['id']) if r['id'] in priority else len(priority),))
    write(HERE / 'FROZEN_EXECUTIONS.json', {'frozen_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'source_catalog_sha256': sha(PACKAGE / 'demos.json'), 'executions': rows, 'excluded': [{'id': 'GENCAD_V2', 'reason': 'Explicitly superseded by GENCAD_V3 in batch4/INDEX.json'}, {'id': 'X35_STATISTICS', 'reason': 'Optional audit tool, not a demonstration'}, {'id': 'LANE_XBREAK_HUNT_2', 'reason': 'Independent review: UNDECIDABLE, not packaged'}, {'id': 'LANE_XBREAK_HUNT_3', 'reason': 'Independent review: REJECT, not packaged'}]})
    prereg = {'claim_type': 'capability', 'frozen_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'capability': 'Run every active packaged demonstration with one command and connect reviewed results to a lab-facing front', 'obstacle': 'Current dispatch catalog omits older batch executables and stops after the first failing route; prior front ran nine selections', 'changed_operation': 'Frozen union of active executable identities; serial isolated copied-package replay; per-demo receipts and failure continuation', 'consumer': 'Prosthetics researcher with CAD/CAM and measurement lab; independent coordinator', 'metrics': [{'name': 'attempted_unique_demos', 'expected': 93, 'tolerance': 0}, {'name': 'unsupported_scientific_promotions', 'expected': 0, 'tolerance': 0}, {'name': 'preserved_source_mutations', 'expected': 0, 'tolerance': 0}, {'name': 'unreported_execution_failures', 'expected': 0, 'tolerance': 0}], 'decision_criterion': 'Delivery complete only after all 93 commands have an individual outcome, timing and original frozen validator result; report failures without changing scientific tolerances', 'numerical_tolerances': 'Each original package validator retains its original tolerance and expected outputs, hash-bound before replay', 'strongest_equally_informed_control': 'Existing independent batch validators and manual exact-review reading; no algorithm superiority claim', 'falsifiers': ['A wrong numeric result must fail its existing comparison', 'A wrong executable identity or missing receipt fails coverage', 'Missing review locator prevents a scientific claim'], 'external_referent': {'kind': 'external_review', 'locator': str(HERE.parent / 'LANE_XREVIEW_BATCH23') + '; ' + str(HERE.parent / 'LANE_XREVIEW_BATCH24') + '; ' + str(HERE.parent / 'LANE_XREVIEW_BATCH25'), 'compared_quantity': 'Version-bound numerical claims, scope and demo readiness', 'refutes_us': True}, 'summary_sufficiency_test': 'Equal exact executable counts with one substituted identity; compare target executable membership. Minimal extension: identity+review+scope per demo.', 'cost': {'preparation': 'Elapsed measured; installed packages and local immutable data retained', 'fit': 0, 'discovery': 'Read source dispatch chain and frozen review files', 'validation': 'All 93 serial executions plus existing fault injections; measured wall time', 'questions': 0, 'fallback': 'Continue after failed demo; preserve error and use direct existing batch runner for omitted legacy IDs'}, 'constraints': {'threads_max': 4, 'intermediate_limit_bytes': 3000000000, 'subagents': 0, 'physical_measurements': 0, 'graph_mutations': False}}
    write(HERE / 'PREREG_X86.json', prereg)
    write(HERE / 'DECOMPOSITION.json', {'idea': 'Make every scoped result executable and lab-readable', 'mechanism': 'Join immutable executable identity, validator, reviewed claim and lab observation contract', 'equation': 'coverage = |attempted executable identities intersect frozen active identities|; PASS_i = original_validator_i(output_i)', 'operation': 'Read batch indices, freeze union, isolate writable replay state, record validator and run costs', 'representation': 'Per-demo record with identity, scope, reviewed hash, command, receipt, duration and scientific unknowns', 'leaves': [{'name': 'Executable union and exact coverage', 'status': 'DERIVED_UNDER_ASSUMPTIONS', 'stop': 'Set membership suffices for administrative coverage; it says nothing about physics'}, {'name': 'Independent batch verdicts', 'status': 'EXTERNALLY_MEASURED', 'stop': 'External review files are the administrative referent; their model scope is retained'}, {'name': 'Scientific model predictions', 'status': 'CONSTITUTIVE_CLOSURE', 'stop': 'No new physical derivation; source model assumptions remain'}, {'name': 'Same-patient force, fabricated crown accuracy, generalization and cold installation', 'status': 'UNKNOWN', 'stop': 'No corresponding measurement or OS-installation trial is authorized or available'}], 'timescales': {'digital_geometry_to_static_model': 'SIMULTANEOUS', 'CAD_to_fabricated_measurement': 'HANDOVER'}, 'resolution': 'Administrative counts/timings PHENOMENOLOGICAL; original physical values retain their finest declared support'})
    for name in ['PREREG_X86.json', 'FROZEN_EXECUTIONS.json', 'DECOMPOSITION.json']:
        (HERE / (name + '.sha256')).write_text(sha(HERE / name) + '  ' + name + '\n')
    clone = DATA / 'clean_package'
    if clone.exists():
        raise RuntimeError('Refuse to overwrite an existing clean instance')
    clone.mkdir()
    bindings = []
    copied = []

    def copy_level(src, dest, depth=0):
        for f in src.iterdir():
            out = dest / f.name
            if f.name == '__pycache__':
                continue
            if f.is_dir() and (not f.is_symlink()):
                if f.name in ['runs', 'work', 'logs']:
                    out.mkdir()
                    continue
                if depth == 0 and f.name.startswith('batch'):
                    out.mkdir()
                    copy_level(f, out, depth + 1)
                    continue
                out.symlink_to(ALIAS / f.relative_to(PACKAGE), target_is_directory=True)
                bindings.append({'source': str(f), 'alias': str(ALIAS / f.relative_to(PACKAGE)), 'copy_link': str(out), 'mode': 'read-only outer namespace'})
            elif f.is_symlink():
                target = os.readlink(f)
                out.symlink_to(target)
            else:
                shutil.copy2(f, out)
                copied.append({'relative_path': str(f.relative_to(PACKAGE)), 'sha256': sha(f), 'bytes': f.stat().st_size, 'mode': f.stat().st_mode & 511})
    copy_level(PACKAGE, clone)
    for q in clone.rglob('LATEST*VERIFICATION*.json'):
        if not q.is_symlink():
            q.unlink()
    write(HERE / 'CLEAN_COPY_MANIFEST.json', {'source': str(PACKAGE), 'clean_instance': str(clone), 'immutable_input_bindings': bindings, 'copied_files': copied, 'scope': 'Fresh writable command/configuration state; unchanged large local inputs shared read-only. Existing host Python environments remain. Not a clean OS installation.'})
    for name in ['PKG_demo48', 'DEMO48_PACKAGE_BATCH4']:
        (DATA / 'mapped_data' / name).mkdir(parents=True)
    print('Frozen', len(rows), 'executions; copied metadata bytes', sum((r['bytes'] for r in copied)), flush=True)
if __name__ == '__main__':
    main()
