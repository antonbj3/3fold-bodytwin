# Styrning LANE_LENS_RESOLUTION — efter r30

## Verified, from the raw integers
I read `HISTORY_STATIC_V1.json` and `HISTORY_MOVEMENT_*_V1.json` and recalculated the span from
integer fields with your own scaling: **0,08875955352806436 fs** for the static mail, exactly your
speech; and 0,0151–0,0203 fs for the ten movement lines, all with the same 26 000 rader. Kvoten 4,3653 mot
the movement MAXIMUM is thus, and the static position is completely outside the distribution of movement.

I tried two explanations and both fell:
* **Drift.** A Straight Trend Explains **1,0 %** of the static span. It is not operation and
  cannot reference away with a sliding window.
* **Enstaka spikar.** Robust measurements show a wider distribution throughout its length, not just in the tail:
  MAD 0,003355 mot 0,001361/0,001333 (faktor 2,5), IQR 0,006887 mot 0,002721/0,002668 (faktor 2,5),
  p99−p1 0,043100 mot 0,009096/0,008831 (faktor 4,8), spann faktor 4,4–5,7.

## What's worth more than your headline
The ratio grows with the quantum: **2,5× at the core and 4,8× i svansarna**The static reference is therefore:
worse in TWO ways — wider core and heavier tails — and a single span number hides that there are two
different defects. Report the decomposition instead of the span. And because it is the static
the record naively chosen as the calibration reference makes this a consequential result: the choice of reference
kostar en faktor 2,5 in the core even if you cut the tails.

## Hindret
Ditt eget: "worst-case calibration remains unacquired". It now has a MEASURED form instead of being
ett tomrum — a worst case to cover a heavier tail cannot be estimated
ur en enda statisk post; `static_runs = 1` mot `motion_runs = 10`.

## Changed operation
Get more static entries from the same archive (`whole_archive_downloaded = False`) and measure if the factor 2,5
in the core holds over several static items. A factor estimated from n = 1 mot n = 10 is directional safe
men saknar bredd; med tre–fyra statiska poster blir den ett tal med intervall.

## Starkaste kontrollen
Samma 26 000 lines, the same exact integer arithmetic, but the movement records treated as if they were
statiska (same estimate, same window)The factor only counts if it survives that the estimater is
identical in both arms.

## Falsifieraren
If the core factor falls below 1,5 when more static entries are received, the static reference is not:
Less in the core and the whole conclusion is reduced to a tail story. Print it before the pickup.
