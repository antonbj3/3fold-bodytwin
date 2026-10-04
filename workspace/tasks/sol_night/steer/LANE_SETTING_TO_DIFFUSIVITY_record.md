# Styrning LANE_SETTING_TO_DIFFUSIVITY — runda 1

Start in PMID 30294056. it is the only named source validated at multiple power levels, and
The whole question depends on whether power and residence time can be separated in the published posts.

Report first, before all modelling: stand power and time as **separata** columns in that source;
Or just as supplied energy? If only as a product the question is decided in the negative direction and then
you should say it directly instead of building a calibration on top of a product.

## Addition on the same day: what a combination must deliver, in numbers
I figured out what it takes, so you don't have to guess what's enough.

Dissemination factor at: `D_eff` is now **3,0718**, making the housing **2,0718×** by themselves
For the housing to be narrower than the actual value, a spreading factor is required. **under 2,0**,
thus two `δ²/t`-poster inom en faktor 2 av varandra. Konkret:

- on the upper `9,6333e−07 m²/s` standing still must the lower up **over 4,81665e−07**
- om den undre `3,1360e−07 m²/s` standing must be the upper bottom **under 6,272e−07**

And the inventory is done: `results/LANE_READ_CREATIVE/CONSUMABLE.json` berries 162 konsumerbara
observabler, 151 of them measurements with: DOI eller PMID — men **exakt noll termiska diffusiviteter**.
The three items that mention conductivity at all are hydraulic, so the thermal chain rests on
**two published depths–time points and nothing else**, which is why the factor is 3 and not 1,2.

The other two named sources measure thermal spread in the liver, muscle and skin (PMID 39161807,
doi `10.1016/j.helvyon.2024.e35266`; PMID 40442741, doi `10.1186/s12893-025-02969-8`). Att ta in
they do two things at once: narrows `D_eff` and test whether it is at all tissue independent.
`δ²/t` differ more between tissues than within a tissue is `D_eff` not a constant for the chain
and then the decision shall be conditional on tissue.
