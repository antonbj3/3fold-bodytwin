"""Re-run the delivered demo while upstream geometry and fulltext paths are unavailable.

This tests reuse of pinned local input/transcription files, not a fresh OS install.
"""
from pathlib import Path
import prepare_demo, literature, run_demo, run_r2, run_r3, plot_results
from metrology import write_json
prepare_demo.SOURCE = Path('/NONEXISTENT_X55_UPSTREAM_CROWN')
prepare_demo.SINTER = Path('/NONEXISTENT_X55_UPSTREAM_SINTER')
literature.CORPUS = Path('/NONEXISTENT_X55_UPSTREAM_CORPUS')
run_demo.run()
run_r2.run()
run_r3.run()
plot_results.run()
write_json(run_demo.R / 'PORTABLE_VERIFICATION.json', dict(status='PASS', scope='Numerical pipeline re-run with upstream source/fulltext locators set to nonexistent paths', reused='Pinned local simulated meshes, labels and exact primary table transcriptions', limitation='Same installed Python/packages and OS; not fresh operating system or physical validation'))
