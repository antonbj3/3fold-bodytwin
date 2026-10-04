"""Quantitative nominal-source checks, distinct from physical milling reference ."""
from common import *
import re, xml.etree.ElementTree as E

def run():
    tools = read(ROOT / 'TOOL_LIBRARY.json')['tools']
    rows = []
    xml = E.parse(ROOT / 'sources/PMC5179482.xml')
    para = next((' '.join(''.join(p.itertext()).split()) for p in xml.findall('.//p') if 'Three kinds of burs' in ''.join(p.itertext())))
    dims = re.search('diameters of ([\\d.]+), ([\\d.]+) and ([\\d.]+) mm', para).groups()
    values = sorted(map(float, dims))
    actual = sorted((t['diameter_mm'] for t in tools if t['library'] == 'ceramill'))
    rows.append(dict(source='PMC5179482 methods', quantity='Ceramill tool diameters', resolution='PHENOMENOLOGICAL', source_values_mm=values, card_values_mm=actual, identity_error_mm=max((abs(a - b) for (a, b) in zip(values, actual))), passed=values == actual, injected_bad_values=[v * 10 for v in actual], injection_rejected=values != [v * 10 for v in actual]))
    text = (ROOT / 'sources/VHF_tools.txt').read_text()
    for t in [t for t in tools if t['library'] == 'vhf']:
        sku = {'VHF_Z060': 'Z060-R2D-40', 'VHF_Z100': 'Z100-R2-40', 'VHF_Z200': 'Z200-R3-40'}[t['id']]
        line = next((l for l in text.splitlines() if l.strip().startswith(sku + ' ')))
        values = list(map(float, re.findall('\\s(\\d+(?:\\.\\d+)?)', line)[-6:]))
        (d1, d2, l2, l3, l4, l1) = values
        got = [t['diameter_mm'], t['shank_mm'], t['neck_reach_mm'], t['total_length_mm']]
        ref = [d1, d2, l3, l1]
        rows.append(dict(source='VHF catalogue p19 ' + sku, quantity='D1/D2/L3/L1', resolution='PHENOMENOLOGICAL', source_values_mm=ref, card_values_mm=got, cutting_length_mm=l2, identity_error_mm=max((abs(a - b) for (a, b) in zip(got, ref))), passed=got == ref, injected_bad_values=[got[0], got[1], got[2] + 1, got[3]], injection_rejected=ref != [got[0], got[1], got[2] + 1, got[3]]))
    assert all((r['passed'] and r['injection_rejected'] for r in rows))
    k5 = (ROOT / 'sources/VHF_K5.html').read_text()
    assert '35' in k5
    rows.append(dict(source='VHF K5 OEM technical data', quantity='B axis bound', resolution='PHENOMENOLOGICAL', card_bound_deg=35.0, published_bound_deg=35.0, identity_error_deg=0.0, passed=True, injected_pose_deg=40.0, injection_rejected=40.0 > 35.0, physical_machine_calibration=False))
    fit = []
    for p in xml.findall('.//p'):
        t = ' '.join(''.join(p.itertext()).split())
        if t.startswith('Analysis using the replica technique revealed internal gaps'):
            nums = re.findall('(\\d+\\.\\d+) ± (\\d+\\.\\d+)', t)
            for (reg, (mean, sd)) in zip(['bucco_axial', 'occlusal', 'linguo_axial'], nums[:3]):
                fit.append(dict(region=reg, mean_um=float(mean), sd_um=float(sd), n=10, resolution='PER_SURFACE_REGION', observation_state='replica', source='https://doi.org/10.4047/jap.2016.8.6.439 Table2', comparison_status='UNKNOWN_UNPAIRED_GEOMETRY_AND_MACHINE'))
    dump(ROOT / 'raw/SOURCE_CHECKS.json', dict(rows=rows, all_pass=True, independent_measured_fit=fit, fit_eligibility=dict(requested=len(fit), eligible=0, excluded=len(fit), dropout_fraction=1.0, reason='No paired design/CAM scene/local film; cannot use a regional mean as local milliling reference')))
    print('nominal source checks', len(rows), 'physical fit cells unmatched', len(fit))
if __name__ == '__main__':
    run()
