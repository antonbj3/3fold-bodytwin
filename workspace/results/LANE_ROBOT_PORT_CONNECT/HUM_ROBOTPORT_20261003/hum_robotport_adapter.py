"""BodyTwin sidecar-v2 extension: strict rational intervals, no evidence admission.

The original202 are referenced, never overwritten or silently rebound.
Legacy ROBOT_PORT_CONNECT validation checks node structure only. This explicit
opt-in extension adds interval semantics and mathematical certificate replay.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from hum_grasp_certificate import rational, verify

FIELDS = {'port_id', 'node_id', 'quantity', 'original_domain', 'direction',
          'source_value', 'source_unit', 'source_regime', 'source_validity_stage',
          'source_provenance_class', 'source_provenance_detail', 'source_locator',
          'source_scope', 'binding_status', 'computed_value', 'computed_unit',
          'computed_provenance_class', 'computed_origin', 'priority',
          'required_computation', 'required_output_unit', 'lower_inputs_required',
          'candidate', 'notes', 'bench_assay_requirements_counted_as_robot_gaps',
          'review_state'}
IDS = {'HUM-CONTACT-SIM::normal_force_per_finger_N',
       'HUM-TOOLGRASP::dry_scalpel_tangential_capacity_N'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_interval(value, *, unit=None, sources=None):
    if not isinstance(value, dict) or not {'lower', 'upper', 'unit', 'uncertainty_kind',
                                         'source_refs', 'certificate_status'} <= value.keys():
        raise ValueError('explicit interval with unit, provenance and status required')
    # Wire format is rational strings: prevents binary-float or bool coercion.
    if any(type(value[k]) is not str for k in ['lower', 'upper']):
        raise ValueError('interval endpoints must be exact rational strings')
    lo, hi = rational(value['lower']), rational(value['upper'])
    if lo > hi or str(lo) != value['lower'] or str(hi) != value['upper']:
        raise ValueError('ordered canonical rational endpoints required')
    if not value['unit'] or (unit is not None and value['unit'] != unit):
        raise ValueError('quantity/unit mismatch')
    if not value['uncertainty_kind'] or not value['source_refs']:
        raise ValueError('uncertainty and source required')
    if sources is not None and any(s not in sources for s in value['source_refs']):
        raise ValueError('unresolved source')
    if value['certificate_status'] not in {'EXACT_CONDITIONAL_MODEL_BOUND',
            'PUBLISHED_INPUT_NOT_DETERMINISTIC_BOUND', 'FIXED_DESIGN_INPUT'}:
        raise ValueError('unsupported certificate status')
    return lo, hi


def validate_catalog(catalog):
    if catalog.get('schema') != 'LANE_ROBOT_PORT_CONNECT/sidecar-v2':
        raise ValueError('BodyTwin sidecar-v2 required')
    if catalog.get('extension_contract') != 'HUM_ROBOTPORT/interval-v1':
        raise ValueError('interval extension required')
    if catalog.get('review_state') != 'PENDING_INDEPENDENT_REVIEW' or catalog.get('scientific_admission') is not False:
        raise ValueError('independent review/admission must be preserved')
    ports = catalog['ports']
    if len(ports) != 2 or {p['port_id'] for p in ports} != IDS:
        raise ValueError('two distinct new HUM records required')
    sources = catalog['sources']
    for source in sources.values():
        if len(source.get('sha256', '')) != 64 or any(c not in '0123456789abcdef' for c in source['sha256']):
            raise ValueError('full source hash required')
        if not source.get('file') or not source.get('scope'):
            raise ValueError('source file and scope required')
    for p in ports:
        if not FIELDS <= p.keys():
            raise ValueError('BodyTwin native port fields missing')
        if p['node_id'] != 'manufacturing:'+p['port_id'].split('::')[0]:
            raise ValueError('unrelated node/quantity binding refused')
        if p['computed_unit'] != 'N' or p['required_output_unit'] != 'N':
            raise ValueError('force units required')
        if p['review_state'] != catalog['review_state'] or p['bench_assay_requirements_counted_as_robot_gaps'] is not False:
            raise ValueError('review/ownership mismatch')
        if p['binding_status'] != 'COMPUTED_CONDITIONAL_MODEL_INTERVAL':
            raise ValueError('conditional binding required')
        if p['computed_provenance_class'] != 'UNSOURCED':
            raise ValueError('internal derived output cannot inherit measurement class')
        validate_interval(p['computed_value'], unit='N', sources=sources)
        if p['computed_value']['certificate_status'] != 'EXACT_CONDITIONAL_MODEL_BOUND':
            raise ValueError('computed output must retain conditional model status')
        for group in ['quantities', 'inputs']:
            for value in p[group].values():
                validate_interval(value, sources=sources)
        origin = p['computed_origin']
        if not origin['commit'] or not origin['file'] or any(s not in sources for s in origin['source_refs']):
            raise ValueError('commit/file provenance required')
        cert = p['certificate']['payload']
        if not verify(cert) or p['certificate']['model_status'] != cert['model_status']:
            raise ValueError('invalid mathematical certificate')
        if p['certificate']['physical_status'] != 'UNKNOWN':
            raise ValueError('no empirical physical certification available')
        for name, pair in [('gap_finger0_m', cert['request']['gaps_m'][0]),
                           ('gap_finger1_m', cert['request']['gaps_m'][1]),
                           ('mu', cert['request']['mu'][0]),
                           ('step_s', [cert['request']['step_s']]*2)]:
            b = p['inputs'][name]
            if (rational(b['lower']), rational(b['upper'])) != tuple(map(rational, pair)):
                raise ValueError('input/certificate drift')
        for name, value in p['inputs'].items():
            expected_unit = 's' if name == 'step_s' else '1' if name == 'mu' else 'm'
            if value['unit'] != expected_unit:
                raise ValueError('input unit drift')
            if name.startswith(('surface_', 'glove_')) and value['certificate_status'] != 'PUBLISHED_INPUT_NOT_DETERMINISTIC_BOUND':
                raise ValueError('geometry dispersion cannot become a deterministic guarantee')
        expected = cert['bounds']['normal_force_per_finger'][0] if 'CONTACT-SIM' in p['port_id'] else cert['bounds']['tangential_capacity']
        for key in ['lower', 'upper', 'unit', 'uncertainty_kind']:
            if p['computed_value'][key] != expected[key]:
                raise ValueError('computed value not bound to certificate')
        # Check every exported readout, not just the headline force.
        for name, value in p['quantities'].items():
            if value['certificate_status'] != 'EXACT_CONDITIONAL_MODEL_BOUND':
                raise ValueError('readout assurance drift')
            expected = cert['bounds'].get(name)
            if expected is None:
                raise ValueError('unlicensed readout')
            if any(value[k] != expected[k] for k in ['lower', 'upper', 'unit', 'uncertainty_kind']):
                raise ValueError('readout not bound to certificate')
        if p['certificate']['empirical_validation'] != 'UNKNOWN' or not p['certificate']['missing_observables']:
            raise ValueError('missing physical observables must remain explicit')
    if set(catalog['connected_ports']) != IDS or len(catalog['connected_ports']) != 2:
        raise ValueError('connected port inventory drift')
    return {'ok': True, 'computed_conditional_extension_ports': 2,
            'empirically_certified_ports': 0, 'scope': 'interval contract and exact reduced-model replay'}


def validate_base(catalog, native_dir):
    native_dir = Path(native_dir)
    base = json.loads((native_dir/'PORT.json').read_text())
    if sha(native_dir/'PORT.json') != catalog['base_catalog']['sha256']:
        raise ValueError('BodyTwin base changed; rebase required')
    if len(base['ports']) != 202 or any(p['computed_value'] is not None for p in base['ports']):
        raise ValueError('base inventory changed')
    if any(p['port_id'] in IDS for p in base['ports']):
        raise ValueError('new IDs collide with original catalog')
    graph_file = Path(base['source_branch'])
    if sha(graph_file) != base['source_branch_sha256']:
        raise ValueError('native branch hash drift')
    engine = native_dir/'engine_r1/anchor_graph_tools.py'
    if sha(engine) != catalog['native_validator_sha256']:
        raise ValueError('native validator drift')
    sys.path.insert(0, str(engine.parent))
    spec = importlib.util.spec_from_file_location('bodytwin_native_robotport', engine)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    result = mod.validate(json.loads(graph_file.read_text()),
                          ledger_path=native_dir/'NOT_PROVIDED_FOLD_LEDGER.jsonl', repo_root=native_dir)
    if not result['ok']:
        raise ValueError(result)
    return {'ok': True, 'original_records_unchanged': 202,
            'original_strict_computed_bindings': 0, 'native_validator': result,
            'scope': 'unchanged ROBOT_PORT_CONNECT node schema, dependencies and load; no empirical or interval validation'}


def validate_contact_port(text, certificate):
    from field_engine.contact_port_v1 import Port, Status, replay_witnesses
    p = Port.from_json(text)
    if p.physical_status != Status.UNKNOWN or p.branch_set.complete:
        raise ValueError('continuous state family must not become a complete singleton')
    if len(p.contacts) != 2:
        raise ValueError('two contact frames required')
    for c, bounds in zip(p.contacts, certificate['request']['gaps_m']):
        if (c.gap.lower, c.gap.upper) != tuple(map(rational, bounds)) or c.gap.unit.value != 'm':
            raise ValueError('gap transport drift')
    if not replay_witnesses(p.branch_set, {'hum_grasp': lambda d: d == certificate and verify(d)}):
        raise ValueError('transport witness replay failed')
    return {'ok': True, 'physical_status': p.physical_status.value,
            'branch_set_status': p.branch_set.status.value}


def compose_catalog(catalog, native_dir):
    """Read-only composed view; retain all202 records and their historical facts."""
    import copy
    validate_catalog(catalog)
    validate_base(catalog, native_dir)
    base = json.loads((Path(native_dir)/'PORT.json').read_text())
    result = copy.deepcopy(base)
    result['ports'] += copy.deepcopy(catalog['ports'])
    result['connected_ports'] += catalog['connected_ports']
    result['composition'] = {'base_port_count': 202, 'extension_port_count': 2,
        'total_port_count': 204, 'extension_contract': catalog['extension_contract'],
        'extension_conditional_bindings': 2, 'extension_empirically_certified_bindings': 0,
        'base_summary_scope': 'unchanged historical202; does not describe composed204',
        'consumer_node_resolution': 'external CAD manufacturing namespace; BodyTwin graph admission pending'}
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('catalog', type=Path)
    ap.add_argument('--native-dir', type=Path)
    ap.add_argument('--compose', action='store_true', help='emit composed204-port JSON to stdout without writing')
    args = ap.parse_args()
    catalog = json.loads(args.catalog.read_text())
    if args.compose:
        if args.native_dir is None: ap.error('--compose requires --native-dir')
        print(json.dumps(compose_catalog(catalog, args.native_dir), ensure_ascii=False, indent=2))
        return
    result = {'extension': validate_catalog(catalog)}
    if args.native_dir:
        result['base'] = validate_base(catalog, args.native_dir)
    print(json.dumps(result, ensure_ascii=False))

if __name__ == '__main__':
    main()
