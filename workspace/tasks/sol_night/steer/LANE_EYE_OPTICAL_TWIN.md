# Styrning r14 — The chain carries; now the surface flooring is all that is missing

## Vad r13 gjorde, verifierat av mig i artefaktfilerna
The chain runs all the way and passes my falsifier. I set the limit at a clinical autorefractory
0,25 D, and you registered it yourself as `autorefractor_threshold.threshold_D = 0,25`:

| led | tal |
|---|---|
| laddning | 1,2195e-15 → 6,0976e-16 mol/m |
| ionic swelling pressure | 3317,79 → 2828,15 I mean, Pa. **−489,65 Pa** |
| pupillmedeltjocklek | 605,582 → 602,039 µm |
| index | 1,3751008 → 1,3756288 |
| **brytkraft** | 60,71208 → **60,28649 D, ie. −0,42559 D** |
| RMS higher order | 0,24997 µm; Strehl 0,12273 → 0,10424 |

**And the rotating control gives −1,3592e-08 D** vid tjockleksfel 7,96e-13 µm, so the whole effect
comes from the non-symmetrical structure. The number is pure: maximum power failure 4,19e-13 D and
rung error 4,97e-15 This is a real result and it should be as such.

## But the sweep says something sharper than the headline
I read all twelve `OPTICAL_*_R13_V1.json`. Effekten **byter tecken mellan betingelser**:

| betingelse | delta_power_D |
|---|---|
| all05 | −0,42559 |
| all15 | **+0,46578** |
| K5610 | −0,56981 |
| K13800 | −0,21103 |
| K30100 | −0,09410 |
| K47300 | −0,05941 |
| source30 | −0,07283 |
| source30increase | +0,07600 |
| **eta = 0** | **+0,01059** |
| **eta = 1** | **−0,86180** |

Fem betingelser klarar 0,25 D and five do not, and the sign reverses. This means that a monotonic
assumed direction would be wrong — and that **The surface flooring eta carries almost the entire outcome**: span over
eta is 0,8724 D, which is **exakt** your own sufficiency gap
`surface_sufficiency.downstream_power_difference_D = 0,87239` vid identitetsfel 0,0 i medeltjocklek, J
and index fields.

## The operation this round
1. **Try whether the eta IS the minimally sufficient enlargement.** That the hatch and the eta-span coincide
   to four decimal places is a strong indication but no proof. Construct two states with identical
   Average thickness, identical index field AND identical eta and measure if the breaking force still differs.
   Don't make it that triple enough, and it's one of the most useful results of the night.
2. **Name eta physically and give it a unit.** What is it an allocation? AV, between which surfaces; and
   **what measurement gives it in a real cornea**? If no one exists: acquisition item with quantity, unit and
   what it determines, in the form of: `notes/ACQUISITION_TARGETS.json`It is the post that determines whether we can
   predict the direction of an intervention's refractive outcome.
3. **Report the change of character as main result, not as a parameter sensitivity.** For a
   surgical twin is the direction that matters: today we can say that the effect is clinical
   big but not in which direction. It is an honest and important message.
4. **Keep the rotating control in each future report.** It is cheap and it shows that the effect
   is not a symmetric artifact.

## Control and falsifier
- **Kontroll:** radiellt vinkelmedel med samma tjocklek — your own, and it gave 1,36e-08 D. Keep it.
- **Falsifierare:** if the triple medium thickness, index field and eta still leaves more than 0,25 D
  downstream difference is not the missing quantity, and then what remains should be named instead.
- **Forbidden:** att redovisa −0,42559 D as the size of the effect without mentioning that all15 ger +0,46578;
  to treat eta as known; to simulate an operation on a rotational symmetry of the cornea.
