"""Exercise the actual CLI with preserved fresh output paths on every replay."""
import csv, datetime, json, subprocess, sys
from pathlib import Path
from lab_alarm import R, csv_rows, write, sha

def run():
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    dst = R / 'raw/cli_replays' / stamp
    dst.mkdir(parents=True)
    rows = csv_rows(R / 'raw/SYNTHETIC_nominal_SIDECAR.csv')
    ex = json.loads((R / 'raw/EXAMPLE_R3_scanner_bias.json').read_text())['timeline']
    fields = list(rows[0]) + ['gauge_sinter_ratio', 'gauge_sinter_sd']
    for (row, v) in zip(rows, ex):
        row.update(post_sinter_ref_length_mm=str(10 * v['scan']), sinter_ratio_sd='.001', gauge_sinter_ratio=str(v['gauge']), gauge_sinter_sd='.0005')
    side = R / 'raw/SYNTHETIC_PAIRED_SCANNER_BIAS_SIDECAR.csv'
    with side.open('w') as h:
        w = csv.DictWriter(h, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    jobs = [('shifted', ['raw/SYNTHETIC_shifted_weibull.csv', '--profile', 'inputs/DEMO_PROFILE.json', '--sidecar', 'raw/SYNTHETIC_shifted_weibull_SIDECAR.csv']), ('unknown', ['raw/SYNTHETIC_nominal.csv']), ('paired', ['raw/SYNTHETIC_nominal.csv', '--profile', 'inputs/DEMO_PAIRED_PROFILE.json', '--sidecar', str(side)])]
    commands = []
    outs = {}
    for (name, args) in jobs:
        path = dst / (name + '.json')
        cmd = [sys.executable, str(R / 'code/lab_alarm.py'), *args, '--output', str(path)]
        res = subprocess.run(cmd, cwd=R, text=True, capture_output=True)
        commands.append(dict(argv=cmd, returncode=res.returncode, stdout=res.stdout, stderr=res.stderr))
        assert res.returncode == 0, res.stderr
        outs[name] = json.loads(path.read_text())
    (a, b, c) = (outs['shifted'], outs['unknown'], outs['paired'])
    q = c['posterior']
    g = dict(shifted_Weibull_link_detected='material_scale' in a['monitor']['first_alarm_at_specimen'], no_profile_stays_UNKNOWN=b['status'] == 'UNKNOWN_NO_FROZEN_LAB_LIKELIHOOD', paired_scanner_link_detected='scanner_bias' in c['monitor']['first_alarm_at_specimen'], paired_gauge_does_not_label_sinter='sinter' not in c['monitor']['first_alarm_at_specimen'], sinter_truth_covered=q['sinter_factor']['credible_95'][0] <= 0.8 <= q['sinter_factor']['credible_95'][1], scanner_bias_truth_covered=q['scanner_bias']['credible_95'][0] <= 0.006 <= q['scanner_bias']['credible_95'][1])
    out = dict(gates=g, paired_first_hits=c['monitor']['first_alarm_at_specimen'], paired_joint_posterior=q['paired_sinter_metrology'], csv_interface_test=True, synthetic=True, replay_path=str(dst), commands=commands)
    write(dst / 'VERIFY.json', out)
    write(R / 'raw/INTEGRATED_CLI_CHECKS.json', out)
    assert all(g.values()), g
    print(json.dumps(dict(gates=g, paired_first_hits=out['paired_first_hits'])))
if __name__ == '__main__':
    run()
