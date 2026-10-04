from common import *
import collections

def run():
    rows = read(ROOT / 'raw/X60_PARTICIPANT_ROWS.json')
    out = []
    for name in sorted({r['participant'] for r in rows}):
        rr = [r for r in rows if r['participant'] == name]
        good = [r for r in rr if r['L1'] == 'PASS']
        ratios = [r['optimization']['model_force_quantile_ratio'] for r in good]
        out.append(dict(participant='X60_' + name, requested=len(rr), case_clusters=len({r['case_key'] for r in rr}), families=sorted({r['family'] for r in rr}), PASS=len(good), status_counts=dict(collections.Counter((r['L1'] for r in rr))), conditional_quantile_ratio_median=float(np.median(ratios)), conditional_quantile_ratio_range=[min(ratios), max(ratios)], calibrated_absolute_force='UNKNOWN', rigorous_transfer_enclosure='MISSING', resolution='POPULATION', conditional_model_resolution='PER_TOOTH'))
    dump(ROOT / 'raw/X60_TABLE.json', out)
    return out
if __name__ == '__main__':
    run()
