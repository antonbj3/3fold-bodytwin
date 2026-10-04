"""Freeze scope and vendor existing reviewed operations before running geometry."""
import ast
import datetime as dt
import hashlib
import json
from pathlib import Path
H = Path(__file__).resolve().parents[1]
DENT = H.parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def main():
    sources = []
    vendor = H / 'code/designgate/vendor'
    vendor.mkdir(parents=True, exist_ok=True)
    (vendor / '__init__.py').touch()
    for name in ['checks.py', 'geometry.py']:
        p = DENT / 'results/PROOF_LANE_GENCAD_V3/code/legacy' / name
        out = vendor / name
        out.write_bytes(p.read_bytes())
        sources.append(dict(source=str(p), source_sha256=sha(p), snapshot=str(out.relative_to(H)), snapshot_sha256=sha(out), operation='verbatim reviewed graph checks/geometry'))
    p = DENT / 'results/LANE_NEXT_P_OCCLUSION_VALIDATION/code/continuous_gap.py'
    s = p.read_text()
    names = {'cross', 'plane', 'inside', 'exact_batch', 'broad_phase', 'extremum'}
    parts = [ast.get_source_segment(s, n) for n in ast.parse(s).body if isinstance(n, ast.FunctionDef) and n.name in names]
    out = vendor / 'continuous.py'
    out.write_text('"""Verbatim function snapshots; experiment runner/imports omitted. See SOURCE_REUSE.json."""\nimport time\nimport numpy as np\nfrom scipy.spatial import cKDTree\n\n' + '\n\n'.join(parts) + '\n')
    sources.append(dict(source=str(p), source_sha256=sha(p), snapshot=str(out.relative_to(H)), snapshot_sha256=sha(out), operation='unaltered function bodies', functions=sorted(names)))
    dump(H / 'SOURCE_REUSE.json', sources)
    prereg = dict(round='R1', claim_type='capability', frozen_utc=dt.datetime.now(dt.timezone.utc).isoformat(), capability='Independent command line decision per STL design rule with coordinate witnesses, explicit UNKNOWN and source-linked thresholds.', obstacle='STL discards units, region identity, material product and acquisition pose; reviewed benchmark accepts only matched piecewise-affine height graphs.', changed_operation='Preserve STL face identity, bind region sidecar to SHA256, enforce graph preconditions, expose inherited checks and continuous triangle-overlap minima without resampling.', consumer='Prosthetics researcher inspecting exported crown/preparation/antagonist before manufacture.', metric={'closed_form_error_mm': 1e-05, 'continuous_full_control_parity_mm': 1e-07, 'transformed_witness_error_mm': 0.0001, 'required_injected_error_rejections': 'all specified attacks', 'false_PASS_count': 0}, decision_criterion='Every declared scope produces a dimensioned answer/witness or UNKNOWN; all false-pass and injected-error controls must reject. Real geometry may legitimately FAIL or UNKNOWN. No claim of complete clinical/manufacturing validation.', strongest_equally_informed_control='Execute inherited PROOF_LANE graph checks directly; continuous unpruned full projected-pair control; analytic parallel-plane and crossing-triangle solutions. These are verification, not a method contest.', external_referent={'kind': 'published_dataset', 'locator': 'https://ditto.ing.unimore.it/bits2bites/; local zip and X18 frozen predictions, SHA256 manifest', 'compared_quantity': 'Continuous signed projected separation of original registered antagonist from the X18 exported conditional crown roof', 'refutes_us': False}, falsifiers=['Any unmarked unit/frame ambiguity yields a dimensional PASS', 'A subtriangle penetration is missed', 'Wrong sidecar hash is accepted', '3Y without product IFU gets a thickness PASS', 'Non-graph geometry inherits horizontal insertion/milling PASS'], full_cost={'preparation': 'source inspection, STL parse/dedup and metadata validation timed', 'fit': 'none', 'discovery': 'implementation wall time not independently metered; UNKNOWN', 'validation': 'injection/control wall time and RSS measured', 'queries': 'each complete CLI wall/CPU/RSS measured', 'fallback': 'UNKNOWN carries exact missing input; manual adjustment minutes UNKNOWN absent measured lab record'}, resource_contract={'threads': 1, 'max_intermediates_bytes': 3000000000, 'GPU': False, 'expected_RAM_GB': '<1; steps >=4GB require heavy_run.sh'}, resolution='PER_POINT witnesses; PER_SURFACE_REGION extrema; PER_TOOTH verdict', timescale='SIMULTANEOUS for geometric checks; HANDOVER if output becomes manufacturing input', sources=sources)
    if (H / 'PREREG_R1.json').exists():
        raise RuntimeError('Refusing to overwrite preregistration')
    dump(H / 'PREREG_R1.json', prereg)
    (H / 'PREREG_R1.sha256').write_text(sha(H / 'PREREG_R1.json') + '\n')
    leaves = [('STL triangles', 'DERIVED_UNDER_ASSUMPTIONS', 'Piecewise planar interpolant is the geometric object; no unknown smooth surface certified.', 'Stop at observed facets; scanner error needs independent scan calibration.'), ('Units and shared pose', 'UNKNOWN', 'STL has no enforced physical length/frame metadata.', 'Explicit declaration with SHA-bound sidecar replaces missing semantics; registration uncertainty remains separate.'), ('Wall graph bound', 'DERIVED_UNDER_ASSUMPTIONS', 'min dz/sqrt(1+L_inner^2) <= minimum Euclidean graph separation <= min dz.', 'Requires matched domain, continuous height graphs, no folded projections; rim/axial walls excluded.'), ('Continuous projected contact', 'DERIVED_UNDER_ASSUMPTIONS', 'Affine gap extrema lie at vertices of each projected triangle intersection.', 'Stops at rigid geometry along a declared axis; force/pressure not identified.'), ('Product wall limit', 'EXTERNALLY_MEASURED', 'Manufacturer IFU specification, not measured fracture probability.', 'Use exact product/indication; class 3Y alone not sufficient.'), ('Cement film limits', 'CONSTITUTIVE_CLOSURE', 'Caller-provided lab acceptance interval in mm; axial nominal gap in R1.', 'No universal cement threshold inferred; actual film after seating requires triple scan.'), ('Insertion and milling', 'DERIVED_UNDER_ASSUMPTIONS', 'Inherited proof only for horizontal prep/intaglio patch and unbounded open halfspace.', 'Full crown path, shaft and fixture absent -> UNKNOWN.'), ('Manual adjustment minutes', 'UNKNOWN', 'No paired STL defect/lab time records read.', 'Need per-rule timed adjustment records; never infer minutes from mm.')]
    dump(H / 'DECOMPOSITION_R1.json', dict(idea='CAD-independent geometric decisions', mechanism='Explicit information contract plus scoped inherited operators', equation='rule = compare observable interval with source/caller threshold; unknown contract => UNKNOWN', operation='STL -> hash-bound regions -> rotated common frame -> continuous facets -> witnesses -> report', representation='triangles + sidecar + interval/status + provenance', leaves=[dict(leaf=a, status=b, relation=c, stop_argument=d) for (a, b, c, d) in leaves]))
    dump(H / 'CURRENT_WORK_STATE.json', dict(lane='X34-stl-design-gate', status='R1_PREREG_FROZEN', latest_gate='No numerical experiment yet', next_operation='Implement metadata and graph STL adapter, then freeze fixture predictions before validation', updated_utc=prereg['frozen_utc']))
if __name__ == '__main__':
    main()
