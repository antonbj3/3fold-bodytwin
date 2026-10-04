# Styrning LANE_EXTERNAL_SOLVER — efter r25

## Hindret, med lanens egna ord
`age_observed = False`, `age_inferred_from_prefix = False`, `repayability_axis =
CHRONOLOGY_MEASUREMENT_DEBT; SIX_OUTPUTS_ALONE_DO_NOT_SUPPLY_AGE`. Hela r25:s certificate rests on
an age that none of the six exits can determine. I read the bisection of `age_stress_r25.py`
rad 24–31 and verified direction towards `AGE_STRESS_R25.json`: 280,3 / 400 / 600 / 900 s FALLER,
1790 s PASSERAR; and `hi` held on the passing side — Thus, the certificate applies to:
**age ≥ 1088,1966 s**, en NEDRE border. Field name `sufficient_age_upper_s` states the opposite. Nothing
consumes the field yet (zero hits in taks/ and in the net); and RESULTS-the section prints both
the end points correctly, so it's a name error and not a bug reported number — but renamed it to
`sufficient_age_lower_s` and add the direction as a custom field in the next round.
Bracketbredden 0,00576 s is exact 1509,7/2¹⁸, i.e. the disintegration of the bisection, not a measured sharpness; say
that in the artifact so that no one reads it as a measured threshold.

## Changed operation
The debt shall be either paid or eliminated, and there are two different structures:
1. **Betala**: find an observant that determines the age from what is actually measured, not from the prefix.
   Then the certificate will be dimensioned instead of conditional.
2. **Eliminera**: to construct a certificate that holds: UNIFORMT over the whole range
   [280,3 s, 1790 s] instead of pointwise at an age. Keeping it, the age issue is gone from
   the forecast and the debt no longer exist.

Route 2 is he who builds something that does not exist. Do it first, and let go 1 bli falsifieraren:
If there is no uniform margin, age is a necessary measurement and it is a result in itself.

## Starkaste kontrollen
Equivalently informed: same observations, same forcering, same band — men prognosen tagen med
the old history replayed. The ability only counts if it matches without replay. `claim_type` is
`capability`, so the control should be the more expensive way that already works, not another method.

## Falsifieraren
Print before driving which uniform margin lowers road 2. `verification_failures = 0`
and 1 547 approved assertions are not a gate until a deliberately broken age falls: insert
negative control in the same run and state that it falls.
