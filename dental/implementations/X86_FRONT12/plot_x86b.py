"""Standalone measured-runtime figure; no scientific status promotion."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
j = json.loads((HERE / 'raw/x86b/FRONT_REPLAY.json').read_text())
rows = j['records']
(fig, ax) = plt.subplots(figsize=(10, 6))
colors = ['#ac3434' if r['status'] != 'PASS' else '#d5922d' if r['demo_id'] in ['PROOF_LANE_FULL_CROWN_R3', 'GENCAD_V6', 'X53'] else '#267c93' for r in rows]
bars = ax.barh([r['demo_id'] for r in rows], [r['wall_seconds'] for r in rows], color=colors)
ax.invert_yaxis()
ax.set_xlabel("Measured wall time (s ), including profile resource queue")
ax.set_title("Twelve front profiles in new copies — actual driving cost")
maxval = max((r['wall_seconds'] for r in rows))
for (bar, r) in zip(bars, rows):
    ax.text(bar.get_width() + maxval * 0.015, bar.get_y() + bar.get_height() / 2, f"{r['status']} · {r['wall_seconds']:.2f} s", va='center', fontsize=9)
ax.set_xlim(0, maxval * 1.5)
ax.spines[['top', 'right']].set_visible(False)
fig.text(0.03, 0.01, "Blue: existing package profile. Orange: fast control + frozen data. Red: execution/gate error.\nPASS applies to the profile; physical validation and remaining research gates are displayed in FRONT.md.", fontsize=9)
fig.tight_layout(rect=[0, 0.07, 1, 1])
fig.savefig(HERE / 'PACKAGE_AND_LAB_X86B.png', dpi=150)
plt.close(fig)
