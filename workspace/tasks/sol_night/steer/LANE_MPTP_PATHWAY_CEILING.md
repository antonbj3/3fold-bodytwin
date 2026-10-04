# Styrning LANE_MPTP_PATHWAY_CEILING — efter r4

## I recalculated r4 and it holds
From the four arm values per preparation (C/A/B/D), I reproduced both your columns with a closed
form I derived myself: **the Bliss infarct is `A·B/C`** and **the union gap is `q = D − (A+B−C)`**.

| preparation | C/A/B/D | ceiling | `A·B/C` | you recorded | `q` | you recorded |
|---|---|---|---|---|---|---|
| healthy rat | 41/28/23/19 | 23 | 15,707317 | 15,7073 | +9 | +9 |
| diabetic rat | 39/35/32/21 | 32 | 28,717949 | 28,7179 | −7 | −7 |
| hyperglycemic Levo/CsA | 53/56/50/35 | 50 | 52,830189 | 52,8302 | −18 | −18 |

Six numbers, six hits. The portable ceiling falls: the combination beats the best single treatment by 15 pp in the
hyperglycemic preparation, and the new HG stratum gives Holm `.00693754`.

## The obstacle, exactly
Each individual p you report relies on approximate Welch inference, and **your own simultaneous
t-box spans zero in all three preparations**: healthy [−21,8610; 39,8610], diabetic
[−37,8610; 23,8610], Levo [−59,6571; 23,6…]. Thus: under the correction you yourself call
strict, none of the three is decisive. More preparations under the same contract do not move that.

## Changed operation
Calculate backward instead of forward. From the published spreads (SD 4/2/10/4 pp for Levo,
SEM and n = 6/arm for the mechanical ones), determine **which n per arm** is required for the simultaneous
box to exclude zero for each preparation. One number per preparation, in animals per arm.

## Strongest control
Additivity `q = 0`, thus exactly the hypothesis that `q` measures deviation from. It is the correct
null hypothesis here, not "no effect".

## Falsifier
If no achievable `n` — say below 200 animals per arm — makes the box exclude zero at the
published spreads, then the ceiling question is **not decidable from summary data**, and that is
the result. Write it with the three n values as evidence instead of looking for a fourth preparation.
