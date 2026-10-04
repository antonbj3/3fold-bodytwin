from common import *
import numpy as np

def main():
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64)
    a = points.copy()
    b = points.copy()
    surface_error = float(np.max(np.abs(a - b)))
    state_a = {'14': 'natural_unrestored', '16': 'restored'}
    state_b = {'14': 'restored', '16': 'natural_unrestored'}
    sa = sum((v == 'restored' for v in state_a.values()))
    sb = sum((v == 'restored' for v in state_b.values()))
    allowa = int(state_a['14'] == 'natural_unrestored')
    allowb = int(state_b['14'] == 'natural_unrestored')
    count_a = np.array([100.0, 0.0, 100.0])
    count_b = np.array([0.0, 100.0, 100.0])
    tests = [dict(summary='Bit-identical outer surface coordinates', identity_error=surface_error, downstream='Native-only material applicability at tooth14', downstream_A=allowa, downstream_B=allowb, downstream_difference=abs(allowa - allowb), minimum_extension='Observed restorative-treatment state on the queried tooth; physical response additionally needs material/boundary/volume measurements'), dict(summary='Total reported restorations', A=sa, B=sb, identity_error=abs(sa - sb), downstream='Whether queried tooth14 is restored', downstream_A=1 - allowa, downstream_B=1 - allowb, downstream_difference=abs(allowa - allowb), minimum_extension='FDI membership; for arbitrary downstream tooth queries preserve the full typed per-tooth observation map'), dict(summary='Total represented tooth vertices', A=float(count_a.sum()), B=float(count_b.sum()), identity_error=float(abs(count_a.sum() - count_b.sum())), downstream='Presence at first queried FDI under100-vertex rule', downstream_A=int(count_a[0] >= 100), downstream_B=int(count_b[0] >= 100), downstream_difference=1, minimum_extension='Count at queried FDI and anatomical identity; total support is insufficient'), dict(summary='Balanced accuracy', A=0.75, B=0.75, identity_error=0.0, state_A={'TN': 50, 'FP': 50, 'FN': 0, 'TP': 100}, state_B={'TN': 100, 'FP': 0, 'FN': 50, 'TP': 50}, downstream='Sensitivity', downstream_A=1.0, downstream_B=0.5, downstream_difference=0.5, minimum_extension='Full class-conditional confusion matrix and patient-cluster support; neither BA nor tooth count alone validates screening')]
    write('raw/SUMMARY_SUFFICIENCY.json', dict(external_referent=dict(kind='our_own_fixture', locator=str(ROOT / 'code/sufficiency.py'), compared_quantity='Representation sufficiency counterexamples only; not clinical detector reference', refutes_us=True), tests=tests, affine_physical_sensitivity='Not reported. Rigorous physical uncertainty enclosure is missing.', pass_gate=all((t['identity_error'] == 0 and t['downstream_difference'] > 0 for t in tests))))
if __name__ == '__main__':
    main()
