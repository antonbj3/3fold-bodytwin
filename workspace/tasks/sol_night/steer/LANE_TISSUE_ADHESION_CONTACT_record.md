# Styrning LANE_TISSUE_ADHESION_CONTACT — efter r1

## My falsifier released, and you were right to stop there.
I wrote: calculate the contact pressure first, and if we are outside the source area applies
the correction does not help us. You calculated it, and I recalculated it:

- declared patch, circular with diameter 1 mm → area **0,7853981633974483 mm²**
- 1 N over that area → **1273,2395447351628 kPa**
- mot gelens 35 kPa: **36,378273×**; mot 40 kPa: **31,830989×**
- against the maximum pressure any material in the source was measured at, 1200 kPa: **1,061033×**

`unguarded_adhesion_transfer_gate = FAIL`. The conclusion that adhesion dominates and Coulomb is
negligible is measured at near-zero pressure on the cutting flank. At our declared standard pressure we are
en faktor 36 over the gel and soon **over** the highest measuring point of the source, so the correction is not
transmissible to that geometry; and `source_comparison_is_not_a_transfer_threshold = True` is the
the correct accounting.

## Hindret, exakt
`pressure_span_kPa = [31,830988618379067; 40743,66543152521]` over 36 fall. Lower end **is**
in the gel area. So the question is not whether the correction applies, but **var** it applies: which of our
actual instrument geometries fall under the pressure range of the source.

## Changed operation
Turn the question. Enter the patch area and load that gives a contact pressure within the source's saturated area, and
determines if that geometry is an instrument geometry someone actually uses — en bredare backe, en
larger contact area, a lower planned load. If yes, the adhesion law is applicable there and shall be built
for that case. If no, the answer is that our printing area is outside all published
adhesion measurement, and then it is **det** to stand at the edge.

## Starkaste kontrollen
Coulomb med µ = 0,295 ± 0,056 vid 1273 kPa, i.e. at our actual pressure. It's the only one
The model is at all measured in our area, and the adhesion law should beat it there or not at all.

## Falsifierare
`second_dominance_sign_flip = True` med 1545 counter-example pairs: if the character change in the dominance order exists
left over the entire pressure range a real instrument geometry covers, so the dominance order is not
a property of the tissue but of where the measurement happens to be taken, and neither model can then support a direction.
