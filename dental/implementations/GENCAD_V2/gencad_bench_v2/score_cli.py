from .scorer import score
from .measurements import score as measurements
from .common import state
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--round', default='R1')
    a = ap.parse_args()
    r = score(a.round)
    if a.round == 'R1':
        measurements()
        state('R1_SCORED', str(r['totals']), 'Preserve R1 and construct changed operation')
    print(r['totals'], flush=True)
