"""CIEDE2000 for positive-b* colours (the selected study's regime)."""
import math
import numpy as np
import mpmath as mp

def _formula(first, second, lib):
    (L1, a1, b1) = first
    (L2, a2, b2) = second
    C1 = lib.sqrt(a1 * a1 + b1 * b1)
    C2 = lib.sqrt(a2 * a2 + b2 * b2)
    cm = (C1 + C2) / 2
    G = (1 - lib.sqrt(cm ** 7 / (cm ** 7 + 25 ** 7))) / 2
    ap1 = (1 + G) * a1
    ap2 = (1 + G) * a2
    cp1 = lib.sqrt(ap1 * ap1 + b1 * b1)
    cp2 = lib.sqrt(ap2 * ap2 + b2 * b2)
    h1 = lib.atan2(b1, ap1)
    h2 = lib.atan2(b2, ap2)
    dh = h2 - h1
    hm = (h1 + h2) / 2
    cpm = (cp1 + cp2) / 2
    lm = (L1 + L2) / 2
    r = lib.pi / 180
    T = 1 - 0.17 * lib.cos(hm - 30 * r) + 0.24 * lib.cos(2 * hm) + 0.32 * lib.cos(3 * hm + 6 * r) - 0.2 * lib.cos(4 * hm - 63 * r)
    SL = 1 + 0.015 * (lm - 50) ** 2 / lib.sqrt(20 + (lm - 50) ** 2)
    SC = 1 + 0.045 * cpm
    SH = 1 + 0.015 * cpm * T
    RT = -2 * lib.sqrt(cpm ** 7 / (cpm ** 7 + 25 ** 7)) * lib.sin(60 * r * lib.exp(-((hm / r - 275) / 25) ** 2))
    dl = (L2 - L1) / SL
    dc = (cp2 - cp1) / SC
    dH = 2 * lib.sqrt(cp1 * cp2) * lib.sin(dh / 2) / SH
    return dl * dl + dc * dc + dH * dH + RT * dc * dH

def de00(first, second):
    assert first[2] > 0 and second[2] > 0
    return math.sqrt(max(0.0, _formula(first, second, math)))

def interval_de00(box, target):
    """Outward enclosure of this exact formula, conditional on rectangular Lab box.

    Boxes are descriptive source mean +/-2SD scenarios, NOT confidence intervals.
    mpmath interval arithmetic encloses dependencies conservatively. No linearization.
    """
    mp.iv.dps = 35
    assert all((lo <= hi for (lo, hi) in box)) and box[2][0] > 0 and (target[2] > 0)
    iv = mp.iv
    lab = [iv.mpf([str(lo), str(hi)]) for (lo, hi) in box]
    ref = [iv.mpf(str(x)) for x in target]
    v = _formula(lab, ref, iv)
    z = iv.sqrt(v if v.a >= 0 else iv.mpf([0, v.b]))
    return [max(0.0, float(np.nextafter(float(z.a), -np.inf))), float(np.nextafter(float(z.b), np.inf))]
