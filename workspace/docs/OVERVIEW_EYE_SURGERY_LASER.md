# Intraocular lens power, electrosurgical incision, laser ablation

Three places where a tissue model has to choose a number a clinician would otherwise choose by
judgement: the power of an implanted lens, the setting on an electrosurgical generator, and the
fluence of an ablating laser. Each choice is scored against a measurement published by someone else.

All results are pending independent review. None is a claim of clinical validation.

## Intraocular lens power

The quantity a surgeon cannot measure before operating is where the implanted lens will sit. From
four preoperative measurements — anterior chamber depth, axial length, lens thickness and central
corneal thickness — the position is predicted by a leave-one-out fit, with no eye in its own fit:
**0.114 mm mean absolute error**, against a spread of **0.284 mm** in the quantity itself.

That position feeds the power calculation. Scored against the power that hindsight shows was correct,
across **89 eyes**:

| reading | mean miss |
|---|---|
| power chosen from the predicted position | **0.535 D** |
| power actually implanted | 1.130 D |

The predicted choice is closer in **36** eyes, further in **22**, and equal in **31**. It falls within
tolerance in **29 of 89**.

Two results below the headline are worth more than it.

The 0.5 D manufacturing grid is often given as the floor under any power decision. It is not large
enough to be that floor: rounding to a 0.5 D lattice contributes **0.125 D**, which is **23.4 %** of
the achieved 0.535 D. The remaining **0.520 D** is not explained by the grid, and for the grid to
explain it the label tolerance would have to be five times the largest published bench scatter [1].

Converting axial length between devices removes **45 %** of the axial-length slope and **0.037 D** of
the error overall — but the twenty shortest eyes get **1.93 times worse**. So the conversion helps on average and hurts the eyes
where the margin is smallest.

For astigmatism, the residual is **0.605 D** using the alternative axis convention against **0.892 D**
using a population estimate.

## Electrosurgical incision

An incision is usually discussed in terms of how hard the instrument is pressed. The energy balance
says that is the wrong variable. For a 1 mm × 1 mm kerf, vaporising the tissue in the cut costs
**2260 J per metre**, while mechanically breaking it costs **4 J per metre** — **0.18 %**. Warming the
kerf from body temperature to boiling adds **10.0 %**.

Over a hand force of 0.5–6 N and contact areas of 0.2–1.0 mm², the mechanical share stays at or below
**0.07**. Pressing harder cannot separate two settings that differ by millimetres in collateral
damage, and an instrument that teaches the hand to press harder is teaching a quantity that does not
carry the outcome.

What does set the damage is how far heat spreads beside the cut. Published thermal damage depth is
**0.56–1.70 mm at 1–3 s** in human tissue [2]. Read as a diffusion depth, that gives an effective
diffusivity of **3.1 × 10⁻⁷ to 9.6 × 10⁻⁷ m²/s**, which is **2.2 to 6.9 times** conduction alone. The
measured damage travels further than conduction permits, which places the evaporative channel in the
measurement rather than in the model.

Two further results constrain what can be calibrated at all.

The generator setting is not the delivered power. At a nominal 25 W, the delivered time-mean power
measures **15.9 to 21.9 W**; at 50 W, **33.3 to 45.5 W**. The lowest delivery is **64 %** and **67 %**
of nominal, and the relative spread is the same at both settings — a constant loss, not one that
worsens with setting.

Geometry dominates electrical configuration. Varying the geometry moves the peak surface temperature
by **27.6 K**; varying the number of internal measurement points moves it by **3.1 K**; moving the
return electrode from 8 to 16 mm moves it by **0.12 K**.

Above a generator reading of 60 the model returns no setting at all. The published damage
correlations are measurably weaker in that range, and a number produced there would be an
extrapolation presented as a choice.

## Laser ablation

The external anchor here is a law rather than a value: a published mass-loss regression with both an
intercept and a slope, held out of the fit [3]. Its two coefficients are measured independently of
each other in the source, and they agree with each other — the slope of 267 µg/J is the reciprocal of
the reported heat of ablation of 3740 J/g to within 0.2 %, which is the first thing to check before
using a law and is rarely reported.

| | threshold | slope |
|---|---|---|
| published, held out | 1.15 J/cm² | 267 µg/J |
| this chain | 1.455 J/cm² | 200.3 µg/J |

The chain is **26.5 % high** on the threshold and **25.0 % low** on the slope, and the two errors have
**opposite signs**. Reported as two relative errors that looks like agreement within a quarter. Inverted
into a decision — choose the fluence that removes a target mass — it does not:

| target removed, µg/cm² | fluence chosen | actually removed | error | error of a fixed fluence |
|---|---|---|---|---|
| 25 | 1.580 J/cm² | 114.7 | +359 % | +1550 % |
| 100 | 1.954 | 214.7 | +115 % | +313 % |
| 400 | 3.452 | 614.5 | +54 % | **+3 %** |
| 800 | 5.448 | 1147.6 | +44 % | −48 % |

The model beats a fixed setting by a wide margin at low targets and **loses to it at 400 µg/cm²**,
where a constant fluence happens to be nearly right. That crossing is the useful part. It says where
a model earns its place and where a device setting already suffices, which is a narrower and more
useful claim than a single accuracy figure.

One convention is in the way and is reported rather than resolved. The same threshold computes as
**1.455 J/cm²** for a Gaussian beam and **1.908** for a top-hat — a **31 %** spread that comes from the
beam profile alone.

## References

1. Hoffer KJ, Calogero D, Faaland RW, Ilev IK. Testing the dioptric power accuracy of
   exact-power-labeled intraocular lenses. *J Cataract Refract Surg* 2009;35(11):1995–9.
   PMID 19878834, doi:10.1016/j.jcrs.2009.06.021
2. Thermal damage depth in human tissue at 1–3 s exposure. PMID 12642261,
   doi:10.1177/03635465030310021601
3. Payne BP et al. Comparison of pulsed CO₂ laser ablation at 10.6 µm and 9.5 µm.
   *Lasers Surg Med* 1998;23(1):1. PMID 9694144,
   doi:10.1002/(SICI)1096-9101(1998)23:1<1::AID-LSM1>3.0.CO;2-T — Table 1 gives the threshold
   1.15 J/cm² (95 % interval 1.03–1.27) and the slope 267 µg/J (interval 242–292) after twelve
   pulses at 1 Hz across 28 specimens. The same paper reports a heat of ablation of 3740 J/g, and
   1/3740 g/J is 267.4 µg/J, so the slope and the heat of ablation are the same measurement read
   two ways — an internal check the chain passes.

Status: pending independent review. No claim of biological or clinical validation is made for any
number above.
