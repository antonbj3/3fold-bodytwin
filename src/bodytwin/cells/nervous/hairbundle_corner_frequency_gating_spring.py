"""Hair-bundle mechanotransduction: corner frequency of the gating spring.

Bullfrog saccular hair-bundle gating-spring stiffness kappa=1.35 mN/m and bundle drag
lambda=2.8 uN.s/m (Nadrowski, Martin & Julicher 2004, fit to STATIC force-displacement +
response-spectrum data) predict a corner frequency f_c = kappa/(2*pi*lambda). The claim under
test says this lands inside the SEPARATELY-measured 5-50 Hz spontaneous-oscillation band, while
its own pre-registered falsifier is weaker: "the corner-frequency-vs-oscillation-band check
failing to replicate (>1 order-of-magnitude mismatch) for any independently-published
same-species kappa/lambda pair" -- i.e. the bar is ORDER-OF-MAGNITUDE agreement, not literal band
membership. This cell checks BOTH and reports them as separate, distinguished numbers rather than
silently substituting the easier one.

GEOMETRY: an overdamped (low-Reynolds-number, hair-bundle scale) linear spring-dashpot has a
single corner (half-power) frequency f_c = kappa/(2*pi*lambda) where kappa is the restoring
stiffness and lambda the drag coefficient -- direct 1-pole low-pass transfer function algebra,
not curve-fitting.

PRE-REGISTERED GATES:
  G1 (arithmetic replication): f_c computed from the claim's kappa,lambda must match the
     claim's stated 76.7 Hz to <1%.
  G2 (claim's falsifier bar, order-of-magnitude vs the 5-50Hz spontaneous-oscillation
     band): f_c / 50 Hz (distance above the band) must be < 10x -- i.e. NOT a >1-order-of-
     magnitude mismatch. This is the actual pre-registered bar, used here as the gate (not the
     looser prose "landing inside").
  INFORMATIONAL (not a gate, reported honestly): whether f_c is LITERALLY inside [5,50] Hz --
     pre-registered to be checked and reported as FALSE/overclaim-flag if it fails, since the
     claim's prose says "inside" while its own falsifier language says "order of magnitude".
  G3 (VOID FLOOR, pre-registered to FAIL): swap kappa and lambda into the WRONG slots of the
     same formula (f_c' = lambda/(2*pi*kappa)) -- a scrambled-parameter control using the
     identical two numbers and the identical formula shape. This must land >1 order of
     magnitude from the 5-50Hz band.

Reads: nothing (both parameters embedded).
Writes: hairbundle_corner_frequency_gating_spring.json (raw numeric results).
"""

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

import os
import json
import math

OUT_PATH = _os.path.join(OUT_ROOT, "hairbundle_corner_frequency_gating_spring",
                         "hairbundle_corner_frequency_gating_spring.json")

KAPPA_N_PER_M = 1.35e-3   # Nadrowski-Martin-Julicher 2004, bullfrog saccular hair bundle, N/m
LAMBDA_NS_PER_M = 2.8e-6  # same source, bundle drag coefficient, N.s/m
CLAIM_FC_HZ = 76.7
OSC_BAND_HZ = (5.0, 50.0)  # independently-measured spontaneous bundle-oscillation band
ORDER_OF_MAG_FACTOR = 10.0


def corner_freq(kappa, lam):
    return kappa / (2.0 * math.pi * lam)


def order_of_mag_distance(f, band):
    lo, hi = band
    if lo <= f <= hi:
        return 1.0  # inside band = 0 log-distance
    if f < lo:
        return lo / f
    return f / hi


def main():
    fc_live = corner_freq(KAPPA_N_PER_M, LAMBDA_NS_PER_M)
    g1_pass = abs(fc_live - CLAIM_FC_HZ) / CLAIM_FC_HZ < 0.01

    dist_live = order_of_mag_distance(fc_live, OSC_BAND_HZ)
    g2_pass = dist_live < ORDER_OF_MAG_FACTOR  # the claim's pre-registered falsifier bar

    literally_inside_band = OSC_BAND_HZ[0] <= fc_live <= OSC_BAND_HZ[1]

    # VOID FLOOR: scrambled parameter slots, same formula shape, same two numbers
    fc_void = corner_freq(LAMBDA_NS_PER_M, KAPPA_N_PER_M)  # swapped: lambda/(2*pi*kappa)
    dist_void = order_of_mag_distance(fc_void, OSC_BAND_HZ)
    void_should_fail = dist_void >= ORDER_OF_MAG_FACTOR
    g3_pass = void_should_fail

    result = {
        "node_id": "hair-bundle mechanotransduction",
        "distinct_from": "the hearing thermal limit (mammalian equipartition floor, different species/params, not rebuilt here); the cochlear Hopf amplifier (active-amplifier compressive nonlinearity, different mechanism, not rebuilt here)",
        "kappa_N_per_m": KAPPA_N_PER_M,
        "lambda_Ns_per_m": LAMBDA_NS_PER_M,
        "f_c_hz_live": round(fc_live, 3),
        "claim_stated_f_c_hz": CLAIM_FC_HZ,
        "osc_band_hz": list(OSC_BAND_HZ),
        "gates": {
            "G1_arithmetic_match_lt1pct": g1_pass,
            "G2_order_of_magnitude_vs_band_own_falsifier_bar": g2_pass,
        },
        "informational_literal_inside_band_claim_prose": literally_inside_band,
        "overclaim_flag": (not literally_inside_band) and g2_pass,
        "void_floor": {
            "description": "kappa/lambda swapped into wrong slots of the same formula (f_c'=lambda/(2*pi*kappa))",
            "f_c_hz_void": round(fc_void, 6),
            "order_of_mag_distance_void": round(dist_void, 2),
            "void_correctly_failed_order_of_mag_bar": g3_pass,
        },
        "verdict": "PASS" if (g1_pass and g2_pass and g3_pass) else "DISAGREE",
        "symmetric_qc_note": (
            "The claim's prose ('landing inside the 5-50Hz band') is LITERALLY FALSE -- "
            "76.7 Hz is not in [5,50] Hz. The claim's pre-registered falsifier (b) in "
            "honest_gaps uses an order-of-magnitude bar instead, which DOES pass (76.7/50 = "
            "1.53x, well under 10x). Reported as an overclaim in the prose, not a refutation "
            "of the underlying order-of-magnitude convergence argument, which holds."
        ),
    }

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
