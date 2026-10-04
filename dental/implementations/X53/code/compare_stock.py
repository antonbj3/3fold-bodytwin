"""Prospective signed-stock comparison with explicit measurement intervals.

This scores the conditional scenario transport, not physical model validity
or clinical fit. The actual scene must be measured and frozen in a new round.
"""
import argparse, csv, math
from common import *

def compare(file):
    freeze_file = ROOT / 'FROZEN_PREDICTIONS_R5.json'
    frozen = read(freeze_file)
    assert digest(frozen['payload']) == frozen['payload_sha256']
    cases = {c['case_id']: c for c in read(ROOT / 'LAB_CASES_20.json')['cases']}
    rows = []
    with open(file, newline='') as f:
        for obs in csv.DictReader(f):
            out = dict(observation=obs, resolution='PER_POINT', timescale='HANDOVER')
            required = ['case_id', 'point_index', 'mesh_sha256', 'tool_library_sha256', 'frozen_predictions_sha256', 'observed_stock_um', 'measurement_bound_um', 'scene_conditioning_verified', 'measurement_locator', 'measurement_kind']
            if any((not obs.get(k) for k in required)):
                out['status'] = 'UNKNOWN_MISSING_FIELDS'
                rows.append(out)
                continue
            try:
                c = cases[obs['case_id']]
                path = DATA / 'R5' / (c['design_id'] + '.json')
                assert sha(path) == frozen['payload']['point_files'][str(path)]
                if obs['mesh_sha256'] != c['mesh_sha256'] or obs['tool_library_sha256'] != sha(ROOT / 'TOOL_LIBRARY.json') or obs['frozen_predictions_sha256'] != sha(freeze_file) or (obs['scene_conditioning_verified'] != 'true'):
                    out['status'] = 'UNKNOWN_UNMATCHED_SCENE'
                    rows.append(out)
                    continue
                a = read(path)
                index = int(obs['point_index'])
                if index < 0 or index >= len(a['records']):
                    raise ValueError('point index outside frozen map')
                pred = a['records'][index]['local_component_stock_um']
                actual = float(obs['observed_stock_um'])
                bound = float(obs['measurement_bound_um'])
                if not math.isfinite(actual) or not math.isfinite(bound) or bound < 0:
                    raise ValueError('finite observation and nonnegative bound required')
                if pred is None:
                    out['status'] = 'UNKNOWN_NO_LOCAL_PREDICTION'
                    rows.append(out)
                    continue
                delta = abs(pred - actual)
                lo = max(0.0, delta - bound)
                hi = delta + bound
                threshold = read(ROOT / 'PREREG_R5.json')['metrics']['overcut_reporting_um']
                status = 'FAIL_SCENARIO_STOCK' if lo > threshold else 'CONSISTENT_SCENARIO_ONLY' if hi <= threshold else 'UNKNOWN_THRESHOLD_UNRESOLVED'
                out.update(status=status, predicted_stock_um=pred, error_interval_um=[lo, hi], threshold_um=threshold, empirical=obs['measurement_kind'] == 'independent_measurement')
            except (KeyError, ValueError, IndexError) as e:
                out.update(status='UNKNOWN_INVALID_OBSERVATION', reason=str(e))
            rows.append(out)
    counts = {s: sum((r['status'] == s for r in rows)) for s in sorted({r['status'] for r in rows})}
    return dict(claim_type='capability', rows=rows, counts=counts, fit_performed=False, prediction_sha256=sha(freeze_file), scenario_transport='FAIL' if counts.get('FAIL_SCENARIO_STOCK', 0) else 'UNKNOWN' if not rows or any((r['status'].startswith('UNKNOWN') for r in rows)) else 'CONSISTENT_PENDING_REVIEW', physical_model_validity='UNKNOWN', limitation='Measurement bound is supplied, not estimated here; scenario input geometry and complete CAM stock have no established rigorous model error enclosure')
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('observations')
    ap.add_argument('--output', default='STOCK_COMPARISON.json')
    args = ap.parse_args()
    dump(args.output, compare(args.observations))
