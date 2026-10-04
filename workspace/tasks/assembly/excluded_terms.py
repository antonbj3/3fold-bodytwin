"""The one place the admission filters get their block list from, loaded from outside the repo.

The list itself is not in the repo on purpose. A list of excluded terms sitting in a public tree
states publicly which subject was withheld, which gives away as much as the withheld material
would. So the terms live in a file only this machine has, the filters read them from there, and
the repo carries the mechanism without the vocabulary.

Fail closed, in the direction that costs least when it is wrong: if the list cannot be read, the
pattern matches everything and every record is rejected. The opposite default is how the material
a filter exists to stop gets through quietly -- a filter that silently passes everything reports
the same green as one with nothing to reject.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

LIST_PATH = Path(os.environ.get('BODYTWIN_EXCLUSION_TERMS',
                                Path.home() / '.bodytwin' / 'exclusion_terms.txt'))
NAMES_PATH = Path(os.environ.get('BODYTWIN_EXCLUDED_NAMES',
                                 Path.home() / '.bodytwin' / 'excluded_names.txt'))
MATCH_EVERYTHING = r'(?s).*'


def _load(path: Path) -> list[str]:
    try:
        lines = path.read_text(encoding='utf-8').splitlines()
    except Exception:
        return []
    return [l.strip() for l in lines if l.strip() and not l.lstrip().startswith('#')]


def _pattern(terms: list[str], word_bounded: bool) -> re.Pattern[str]:
    if not terms:
        return re.compile(MATCH_EVERYTHING)
    parts = [(r'\b' + re.escape(t) + r'\b') if word_bounded else re.escape(t).replace(r'\ ', r'\s+')
             for t in terms]
    return re.compile('|'.join(parts), re.I)


def block_pattern() -> re.Pattern[str]:
    """Subject matter that does not enter. Substring match: the terms are stems on purpose."""
    return _pattern(_load(LIST_PATH), word_bounded=False)


def person_pattern() -> re.Pattern[str]:
    """Named individuals. Word bounded, because a short name is a substring of ordinary words."""
    return _pattern(_load(NAMES_PATH), word_bounded=True)


def available() -> bool:
    return bool(_load(LIST_PATH))


def reason() -> str:
    if available():
        return f'block list loaded from {LIST_PATH}'
    return (f'block list missing at {LIST_PATH}; filtering fail-closed, every record rejected '
            f'(set BODYTWIN_EXCLUSION_TERMS to point at it)')
