from common import *

def run():
    key = split()['dev_round1'][0]
    ts = tasks(key)
    refs = dev_reference(key)
    rows = []
    for orig in ts[::4]:
        if orig['status'] != 'READY':
            continue
        t = scene(orig)
        ok = np.isfinite(t['ceiling']) & np.isfinite(refs[t['family']])
        if not ok.any():
            continue
        near = t['prior'].copy()
        far = t['prior'].copy()
        near[ok] = t['ceiling'][ok] - 0.05
        far[ok] = t['ceiling'][ok] - 0.15
        m1 = metric(t, near, refs[t['family']])
        m2 = metric(t, far, refs[t['family']])
        passcheck = m1['pred_contact_mm2'] > 0 and m2['pred_contact_mm2'] == 0
        rows.append(dict(family=t['family'], injected_near_area_mm2=m1['pred_contact_mm2'], injected_far_area_mm2=m2['pred_contact_mm2'], rejected_band_confusion=passcheck, resolution='PER_SURFACE_REGION', facit='native registered dev ceiling/reference; computational band injection not a new measurement'))
    if not rows or not all((r['rejected_band_confusion'] for r in rows)):
        raise RuntimeError('contact-band control failed')
    dump(ROOT / 'raw/CONTACT_CONTROLS.json', dict(controls=rows, PASS=True))
    print('contact-band faults', len(rows), 'PASS')
if __name__ == '__main__':
    run()
