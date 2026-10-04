"""Validate the existing native/prosthetic correction against its source identity."""
import re, zipfile
from common import ZIP, digest

def source_guard(row, polarity):
    if polarity != 1 or row['source_member'].split('/')[0] != row['case_id']:
        return False
    match = re.search('absence of tooth\\s+(\\d+).*prosthetic element is present', row['clause'], re.I)
    if not match or int(match.group(1)) != row['fdi']:
        return False
    with zipfile.ZipFile(ZIP) as z:
        if row['source_member'] not in z.namelist():
            return False
        blob = z.read(row['source_member'])
    return digest(blob) == row['source_sha256'] and row['clause'].lower() in blob.decode('utf8').lower()
