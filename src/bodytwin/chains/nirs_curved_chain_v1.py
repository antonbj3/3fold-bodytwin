"""Conditional curved-fixture chain using frozen propagation and declared bands."""
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.chains.nirs_chain_v1 import leg,SEED,DRAWS,NOMINAL,BANDS


def main():
    paths=np.load(ROOT/'reports/photon_curved_detector_paths_l4.npz')['paths_0']
    legs=[leg(paths),leg(paths)]
    gates=dict(full_chain_repeat_exact=legs[0]==legs[1],finite_physical=all(r['finite_physical'] for r in legs),
               zero_band=all(r['zero_band_exact'] for r in legs),zero_blood=all(r['zero_blood_ratio']==1 for r in legs),
               oxygen_direction=all(r['oxygen_direction'] for r in legs),
               sample_composition=all(r['joint_span']<=r['sum_isolated_spans'] for r in legs),
               grid_composition=all(r['grid_span']<=r['grid_isolated_span_sum'] for r in legs),
               samples_inside_measured_grid=all(r['random_within_grid'] for r in legs))
    report=dict(legs=legs,gates=gates,seed=SEED,draws=DRAWS,nominal_inputs=NOMINAL.tolist(),uniform_half_bands=BANDS.tolist(),
                clinical_calibration_certified=False,scope='Conditional synthetic curved fixture; excludes anatomy, extinction, finite photon and clinical uncertainty.')
    (ROOT/'reports/nirs_curved_chain_v1.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
