"""Resolve relocated source paths from explicit input configuration."""
import os
from pathlib import Path
import re
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]

class MissingInput(ValueError):
    pass

def expand(value):
    defaults = {'DENTAL_IMPLEMENTATIONS': str(ROOT/'implementations'),
                'DENTAL_PYTHON': sys.executable,
                'DENTAL_WORK_ROOT': str(Path(tempfile.gettempdir())/'dental-release-work'),
                'DENTAL_CASE_ID': 'synthetic'}
    def replacement(match):
        key = match[1]
        result = os.environ.get(key, defaults.get(key))
        if result is None:
            raise MissingInput('Set '+key+' to the externally supplied input location; see docs/DATA.md')
        return result
    return re.sub(r'@([A-Z0-9_]+)@', replacement, value)
