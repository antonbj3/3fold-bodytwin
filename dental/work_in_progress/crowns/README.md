# Work in progress: full crowns

This folder is unfinished work. It is not part of the reviewed demos.

## Where it stands (4 October 2026)

- **Root cause found.** Earlier full-crown attempts failed because the rule-based virtual preparations reached outside the natural tooth at the cervical margin. The insertion cone cut through the tooth surface on 18 of 18 teeth.
- **Virtual preparations kept inside the tooth.** With this constraint, 7 of 18 teeth have a preparation and crown assembly that passes the geometry checks: containment, insertion, cement film and wall. No crown on these teeth has passed the full chain yet.
- **Scanned clinical preparations.** One crown passes the digital checks for cement gap, wall thickness, closed shell, self-intersection, export and straight insertion. This is material scenario KATANA ML on scanned preparation 6158-21 from Alsheghri et al. (2025, MIT licence).
  - An independent check reproduces these results under the stated scale assumption of 1 mm per source unit. The physical scale of the source scans is not verified.
  - Neighbouring teeth and the antagonist are not part of that dataset, so the crown has no contact or occlusion check.

## Files

- `real_preparation_6158-21/`: the crown certificate, the independent check and the source data licence.

## Status

Not a qualified crown. Not for manufacturing. Physical validation is pending.
