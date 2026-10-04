# Direction LANE_SPINE_LIGAMENT_STIFFNESS — after r5

## You answered the insensitivity question, and the answer is no
The document's own statement was that the decision is independent of the ligament stiffness within published
dispersion. That's not true. The work **changes signs**: `+9,991910446631359e−05 J` at common
offset −0,02 and `−0,00016259115993616925 J` at +0,02, i.e. a sign change with the amplitude ratio
1,6272. The critical offset bracket is `[0,003116457225125695; 0,003516457225125695]`, a width of
exactly **4,0e−04**, and translated into slack: **−0,10654162637971365 to −0,01912365363549752 mm**.

I recalculated it in micrometers: **19,12 to 106,54 µm**, a range of 87,42 µm. The decision
therefore turns within a tenth of a millimeter of ligament slack. `uniform_empirical_prediction = FAIL` is
the correct conclusion.

## The obstacle, exactly
A sign change at 19–107 µm is only useful if slack can be measured better than that. We know that
not. Nine linked hold-out cases are unavailable (`coupled_held_cases_unavailable = 9`), and
`human_ligament_admission`, `load_driven_continuous_spectrum_rotation_enclosure` and
`instrument_metrology` is all UNKNOWN — even though the archive at 25 733 007 changed with 97 members
has been downloaded since r4.

## Changed operation
Find the measurement uncertainty of ligament slack in the 97 archive members or in published literature, in
micrometer, and set it against 19–107 µm. Three outcomes and all three are answers:

1. the uncertainty is **under 19 µm** → the decision is decidable and should be run with slack as input;
2. it lies **within 19–107 µm** → the decision is decidable only in parts of the interval, specify which;
3. it is **over 107 µm** → the sign change is immeasurable and the decision must not carry a direction. Say it
   plainly, with the number.

## Strongest control
The decision with slack was set to zero, i.e. without the quantity the entire character change depends on. If not
differs from the slack-dependent decision outside the critical bracket is slack load-bearing only
in a 87 µm wide window, and it should be said with that delineation.

## Falsifier
If the character change disappears when you use the 97 archive members' own geometry instead
nominal, then it is an artifact of the nominal geometry and not a property of the ligament.
