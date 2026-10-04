# Styrning LANE_SETTING_TO_DIFFUSIVITY — efter r2

## You found why the setting cannot be calibrated straight
`primary_delivered_power_at_10_25W_time_means_W = [9,43; 25,08]` — vid nominellt 10 and 25 W is the
faktiskt levererade tidsmedeleffekten 9,43 and 25,08 W, that is, −5,7 % and +0,32 %Close to nominal.

Men vid nominellt 50 W is the reported transient `[50; 13; 36] W` and
`primary_50W_time_mean_W = None`I counted the ratio: **50/13 = 3,8462**. Den levererade effekten
thus falls to a quarter within the same nominal setting, and the source reports no
time average effect there.

**That is a second, independent reason for the same refusal.** My Decision Refuses From Attitude 60
because the correlation is losing its grip. Now there is a mechanism: over about 25 W is the setting
not the delivered power. Two approaches, the same conclusion.

## Hindret, exakt
You have two layer temperatures, `71,75293210216537` and `56,91036995856369 °C`, med gapet
**14,842562143601683 K**, and a support band `[24,056905302588405; 72,0019017909824] °C`. Den
critical support half-width is **0,19567759147916347 mm**So, toler.ансone in depth is two tenths
millimetres while the temperature is indefinite within 14,8 K. A calibration cannot be narrower than the
Worst of them.

## Changed operation
Kalibrera mot **levererad tidsmedeleffekt**, never against the setting, and enter for each source entry
which of the two it reports, where the time-average effect is missing — som vid 50 W — is the entry is not
usable for calibration and shall be excluded for that reason, not approximated by the transient.

Dela sedan `D` per tissue as previously directed requested. Liver, muscle, skin, spread within each
sig, med uteslutningarna redovisade.

## Starkaste kontrollen
Calibration against the nominal setting, that is, the naive choice. If the tissue-specific calibration
against delivered power does not narrow the spread during **14,897233653406495** you measured in r1 har
The trade didn't buy anything.

## Falsifierare
If the mean time effect is missing in more than half of the items above 25 W, so goes the high end of
setting scale not to calibrate out of published data at all. Say it with the number of entries, and let
the refusal of the decision remains the only management of that area.
