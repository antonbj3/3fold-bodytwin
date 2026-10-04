"""New input-information factorial, with independent all-face retraction fault."""
from geometry import *
from experiment import local_tri

def run():
    fr = json.loads((H / 'FROZEN_PREDICTIONS_R4.json').read_text())
    pred = json.loads((H / 'raw/PREDICTIONS_R4.json').read_text())
    rows = json.loads((H / 'raw/RESULTS_R4_ROWS.json').read_text())
    assert sha(H / 'PREREG_R4.json') == fr['prereg_sha256']
    assert sha(H / 'raw/PREDICTIONS_R4.json') == fr['predictions_sha256']
    assert all((sha(p['path']) == p['sha256'] for p in fr['files']))
    good = [r for r in rows if r['status'] == 'SCORED']
    assert good
    for row in good:
        for arm in ['generic_with', 'personal_with']:
            assert row['continuous'][arm]['maximum_penetration_mm'] <= 1e-07
    r = good[0]
    p = next((p for p in pred if p['case'] == r['case'] and p['fdi'] == r['fdi']))
    z = np.load(p['file'])
    (data, _) = pair(r['case'])
    U = local_tri(data['upper']['tri'], z['xy'])
    fault = signed_gap(U, z['xy'], z['z_personal_with'] + 1, z['faces'])
    rejected = fault['maximum_penetration_mm'] > 1e-07
    assert rejected
    dump(H / 'raw/VERIFICATION_R4.json', dict(status='PASS', scored_sites=len(good), all_frozen_source_prediction_hashes_match=True, every_informed_roof_nonpenetrating=True, injected1mm_height_rejected=rejected, external_referent=dict(kind='closed_form', locator='raw/VERIFICATION.json independent halfplane LP; raw/VERIFICATION_R2.json full polygon inequality equivalence', compared_quantity='Affine minimum over every projected roof/antagonist triangle intersection', refutes_us=True), dental_information_gate=json.loads((H / 'rounds/R4.json').read_text())['summary']['decision']))
if __name__ == '__main__':
    run()
