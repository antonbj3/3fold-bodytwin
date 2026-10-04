"""Only deidentified case locators and relevant clinical spans leave the ZIP.

No demographics, history, complaints, treatment text or names are exported.
"""
import csv
import gzip
import io
import re
import struct
import zipfile
from common import ARCHIVE, MEMBER, sha
FDI = re.compile('(?<![\\dA-Za-z])([1-4][1-8])(?![\\dA-Za-z])')
FDI_R2 = re.compile('(?<![\\dA-Za-z])([1-4][1-8])(?=$|[^\\dA-Za-z]|[DMBL](?=[^A-Za-z]|$)|CBCT|CT)')
PATTERNS = {'PD': re.compile('(?:periodontal\\s+pockets?|\\bPD\\b)[^.;]{0,75}?(\\d+(?:\\.\\d+)?)\\s*(?:-\\s*(\\d+(?:\\.\\d+)?))?\\s*mm', re.I), 'CAL': re.compile('attachment\\s+loss[^.;]{0,35}?(\\d+(?:\\.\\d+)?)\\s*(?:-\\s*(\\d+(?:\\.\\d+)?))?\\s*mm', re.I), 'RBL_LOCATION': re.compile('(?:alveolar\\s+(?:bone|ridge)|bone)[^.;]{0,80}?(?:resorp\\w*|absor\\w*)[^.;]{0,70}?(?:apical\\s+1/3|root\\s+(?:apex|tip)|middle\\s+1/2|middle\\s+third)', re.I), 'PULP_PENETRATION_REPORTED': re.compile('(?:deep\\s+)?caries[^.;]{0,45}?penetrat\\w*\\s+(?:the\\s+)?pulp', re.I), 'NEAR_PULP_REPORTED': re.compile('(?:caries|carious)[^.;]{0,75}?(?:near|close\\s+to)\\s+(?:the\\s+)?(?:pulp|marrow)', re.I)}
PATTERNS_R2 = {'FURCATION': re.compile('(?:furcation|bifurcation)[^.;]{0,35}?(III|II|I|[1-3])(?:\\s*(?:degree|class))?', re.I), 'BONE_LOSS_REPORTED': re.compile('(?:alveolar\\s+bone|bone)[^.;]{0,50}?(?:resorp\\w*|absor\\w*)', re.I), 'CARIES_DEPTH_REPORTED': re.compile('deep\\s+caries|caries\\s+(?:damaged|reaches|damage\\s+reaches)[^.;]{0,45}?dentin(?:e)?\\s+(?:layer|third)', re.I), 'PULP_NOT_REACHED_REPORTED': re.compile('caries[^.;]{0,35}?(?:does\\s+not|did\\s+not)\\s+reach\\s+(?:the\\s+)?(?:pulp|medullary)\\s+cavity', re.I)}

def source_rows():
    with zipfile.ZipFile(ARCHIVE) as z:
        data = z.read(MEMBER)
        infos = {x.filename: x for x in z.infolist() if x.filename.endswith('.nii.gz')}
        rows = list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'))))
    return (rows, infos, sha(data))

def case_id(filename):
    return 'M-' + sha(('X57-MMDental-v1:' + filename).encode())[:12]

def header(z, member):
    with z.open(member) as f, gzip.GzipFile(fileobj=f) as g:
        h = g.read(352)
    endian = '<' if struct.unpack('<i', h[:4])[0] == 348 else '>'
    if struct.unpack(endian + 'i', h[:4])[0] != 348:
        raise ValueError('Unsupported NIfTI header')
    dim = struct.unpack(endian + '8h', h[40:56])
    pixdim = struct.unpack(endian + '8f', h[76:108])
    if h[344:348] not in [b'n+1\x00', b'ni1\x00']:
        raise ValueError('Wrong NIfTI magic')
    return {'shape': list(dim[1:4]), 'spacing_header': list(pixdim[1:4]), 'xyzt_units_byte': h[123], 'header_sha256': sha(h), 'physical_frame': 'NIfTI sform/qform retained in source; no landmarks measured', 'resolution': 'PER_ARCH', 'data_quantity': 'header voxel grid only'}

def tooth_scope(text, start, mode='R1'):
    preceding = list((FDI_R2 if mode == 'R2' else FDI).finditer(text[:start]))
    if not preceding:
        return ([], 'PER_ARCH', 'NO_PRECEDING_FDI')
    last = preceding[-1]
    if start - last.end() > 260:
        return ([], 'PER_ARCH', 'FDI_TOO_DISTANT')
    group = [last]
    for m in reversed(preceding[:-1]):
        between = text[m.end():group[0].start()]
        if mode == 'R2':
            between = re.sub('^[DMBL](?=[,*/\\s])', '', between)
        if re.search('[A-Za-z]', between) or len(between) > 12:
            break
        group.insert(0, m)
    local = text[last.end():start]
    if mode == 'R2' and re.search('whole\\s+mouth|throughout.*mouth|full\\s+mouth|anterior\\s+teeth|posterior\\s+teeth|overall', local, re.I):
        return ([], 'PER_ARCH', 'EXPLICIT_ARCH_REGION')
    if mode == 'R2' and re.search('adjacent\\s+teeth|remaining\\s+teeth|rest\\s+of\\s+(?:the\\s+)?teeth', local, re.I):
        return ([], 'PER_ARCH', 'REFERENT_CHANGED_TO_OTHER_TEETH')
    teeth = sorted(set((int(m.group(1)) for m in group)))
    resolution = 'PER_SURFACE_REGION' if re.search('mesial|distal|buccal|lingual', local[-80:], re.I) else 'PER_TOOTH'
    return (teeth, resolution, 'PRECEDING_EXPLICIT_FDI_GROUP')

def extract_row(row, row_number, csv_hash, mode='R1'):
    text = row['Oral Check'].replace('\\n', ' ').replace('\n', ' ')
    patient = case_id(row['Filename'])
    found = []
    patterns = {**PATTERNS, **(PATTERNS_R2 if mode == 'R2' else {})}
    for (quantity, pattern) in patterns.items():
        for m in pattern.finditer(text):
            if quantity in ['BONE_LOSS_REPORTED', 'CARIES_DEPTH_REPORTED']:
                finer = ['RBL_LOCATION'] if quantity == 'BONE_LOSS_REPORTED' else ['PULP_PENETRATION_REPORTED', 'NEAR_PULP_REPORTED']
                if any((a.start() <= m.start() < a.end() for key in finer for a in PATTERNS[key].finditer(text))):
                    continue
            (teeth, resolution, scope_rule) = tooth_scope(text, m.start(), mode)
            surface = None
            if mode == 'R2' and teeth:
                region_text = m.group()
                if quantity == 'CARIES_DEPTH_REPORTED':
                    region_text = text[m.start():m.end() + 60]
                sm = re.search('mesial|distal|buccal|lingual', region_text, re.I)
                if sm:
                    resolution = 'PER_SURFACE_REGION'
                    surface = sm.group().lower()
                if quantity == 'FURCATION':
                    resolution = 'PER_SURFACE_REGION'
                    surface = 'furcation'
                if surface is None:
                    tokens = list(FDI_R2.finditer(text[:m.start()]))
                    if tokens:
                        suffix = text[tokens[-1].end():tokens[-1].end() + 1]
                        if suffix in ['D', 'M', 'B', 'L']:
                            resolution = 'PER_SURFACE_REGION'
                            surface = {'D': 'distal', 'M': 'mesial', 'B': 'buccal', 'L': 'lingual'}[suffix]
            reason = None
            interval = None
            outquantity = quantity
            if quantity in ['PD', 'CAL']:
                interval = [float(m.group(1)), float(m.group(2) or m.group(1))]
                if interval[0] > interval[1] or not 0 <= interval[0] <= interval[1] <= 20:
                    reason = 'INVALID_MM_INTERVAL'
                if quantity == 'CAL' and re.search('averag', m.group(), re.I):
                    reason = 'MEAN_CAL_NOT_WORST_SITE'
                if mode == 'R2' and re.search('not\\s+(?:reached|detected)|no\\s+periodontal', m.group(), re.I):
                    reason = 'NEGATED_OBSERVATION'
            elif quantity == 'RBL_LOCATION':
                outquantity = 'RBL'
                interval = [2 / 3, 1] if re.search('apical\\s+1/3', m.group(), re.I) else [0.5, 1] if re.search('middle', m.group(), re.I) else [1, 1]
            elif quantity == 'FURCATION':
                token = m.group(1).upper()
                grade = {'I': 1, 'II': 2, 'III': 3}.get(token, int(token) if token.isdigit() else 0)
                interval = [grade, grade]
            before = text[max(0, m.start() - 65):m.start()]
            if mode == 'R2' and re.search('no\\s+(?:obvious\\s+)?(?:alveolar|bone|caries)|no\\s+(?:caries|pulp)', before[-30:], re.I):
                reason = 'NEGATION_SCOPE_AMBIGUOUS'
            found.append({'observation_id': 'O-' + sha(f'{csv_hash}:{row_number}:{quantity}:{m.start()}'.encode())[:16], 'patient_id': patient, 'snapshot_row': row_number, 'chain': 'BD10' if quantity in ['PD', 'CAL', 'RBL_LOCATION', 'FURCATION', 'BONE_LOSS_REPORTED'] else 'BD08', 'quantity': outquantity, 'interval': interval, 'unit': 'mm' if quantity in ['PD', 'CAL'] else 'root-length fraction semantic scenario' if quantity == 'RBL_LOCATION' else 'furcation grade' if quantity == 'FURCATION' else 'category', 'knowledge_status': 'CONSTITUTIVE_CLOSURE' if quantity == 'RBL_LOCATION' else 'EXTERNALLY_MEASURED_REPORT_ONLY', 'match_text': m.group(), 'teeth': teeth, 'resolution': resolution, 'surface': surface, 'scope_rule': scope_rule, 'timescale': 'SIMULTANEOUS', 'accepted': reason is None, 'rejection': reason, 'locator': {'archive': str(ARCHIVE), 'member': MEMBER, 'csv_sha256': csv_hash, 'row_1based_including_header': row_number, 'column': 'Oral Check', 'original_column_sha256': sha(row['Oral Check'].encode()), 'canonicalization': 'literal backslash-n and LF replaced by single space', 'span_start': m.start(), 'span_end': m.end()}, 'uncertainty': 'Text semantics/clinical measurement error UNKNOWN; intervals are reported range or semantic scenarios', 'measurement_debt': 'Independent CEJ/crest/apex landmarks needed to replace root-third scenario' if quantity == 'RBL_LOCATION' else None, 'local_context_for_review': text[max(0, m.start() - 75):min(len(text), m.end() + 50)]})
    return found
