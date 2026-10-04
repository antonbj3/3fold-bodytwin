# Translation validation commands

Run from the dental directory with the versions declared in requirements.txt.

```bash
python3 -m venv --system-site-packages .venv
PYTHONNOUSERSITE=1 .venv/bin/python -B -m pip install --no-cache-dir -r requirements.txt
export DENTAL_PYTHON="$PWD/.venv/bin/python"
export TMPDIR="$PWD/.venv/runtime"
mkdir -p "$TMPDIR"
PYTHONDONTWRITEBYTECODE=1 bash run_all.sh
python3 -B tools/scrub.py --root .
git diff --check
```

The first scrub rejected a translated geometric term. Use “cavity” for that geometry. System SciPy 1.8.0 caused 6 setup errors in a 47-test run. With the declared dependencies, run_all.sh passed 53 tests and 23 computation profiles.

The text scan retains 43 files as intentional exceptions: machine values, executable text rules, stored lookup keys, the English glossary and accented citation names. It leaves 0 files with untranslated prose. TRANSLATION_AUDIT.json lists every exception and translated file.

Refresh release_sha256 in provenance/CODE_MANIFEST.json and file digests in FROZEN_PREDICTIONS.json and FROZEN_CAPABILITIES.json. Rebuild RELEASE_FILES.json from every release file except itself. Preserve source and review hashes and all numerical criteria.
