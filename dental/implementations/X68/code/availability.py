"""Read archive schemas/headers and aggregate text only; no patient text exported."""
from dental_release.paths import expand as _release_expand
import pathlib, json, zipfile, gzip, struct, csv, io, re, hashlib, time
ROOT = pathlib.Path(__file__).resolve().parents[1]
TF = pathlib.Path(_release_expand('@DENTAL_CORPUS_ROOT@/geometry/ToothFairy2/ToothFairy2_Dataset.zip'))
MM = pathlib.Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/datasets/dental_3fold_extra/MMDental/MMDental.zip'))

def write(n, o):
    (ROOT / n).write_text(json.dumps(o, indent=2, ensure_ascii=False) + '\n')

def main():
    start = time.perf_counter()
    out = {}
    with zipfile.ZipFile(TF) as z:
        names = z.namelist()
        ds = next((n for n in names if n.endswith('dataset.json')))
        meta = json.loads(z.read(ds))
        labels = meta.get('labels', {})
        im = next((n for n in sorted(names) if '/imagesTr/' in n and n.endswith('.mha')))
        with z.open(im) as f:
            head = f.read(4096)
        h = head[:head.index(b'ElementDataFile')].decode()
        spacing = re.search('ElementSpacing = ([^\\n]+)', h).group(1)
        out['ToothFairy2'] = {'archive_path': str(TF), 'archive_size_bytes': TF.stat().st_size, 'training_image_count': sum(('/imagesTr/' in n and n.endswith('.mha') for n in names)), 'test_image_count': sum(('/imagesTs/' in n and n.endswith('.mha') for n in names)), 'count_scope': 'Actual ZIP members, not published challenge totals; local archive has imagesTr/labelsTr only', 'metadata_numTraining': meta.get('numTraining'), 'metadata_numTest': meta.get('numTest'), 'training_label_count': sum(('/labelsTr/' in n and n.endswith('.mha') for n in names)), 'labels': labels, 'header_spacing_mm': spacing, 'metadata_member': ds, 'metadata_sha256': hashlib.sha256(z.read(ds)).hexdigest(), 'cortex_inner_boundary_label': False, 'regional_cortex_thickness': 'Not supplied as measurement; images could be re-segmented, requires validated inner boundary at limited voxel resolution', 'absolute_density': 'UNKNOWN_UNCALIBRATED_GRAY; prior self-calibration failed', 'final_bore_profile': 'ABSENT', 'thread_core_geometry': 'ABSENT', 'paired_insertion_torque': 'ABSENT_FROM_IMAGE_ARCHIVE', 'license': 'CC BY-SA 4.0; verify attribution/ShareAlike for derived redistributed dataset'}
    with zipfile.ZipFile(MM) as z:
        names = z.namelist()
        nm = next((n for n in names if n.endswith('medical_records.csv')))
        b = z.read(nm)
        text = b.decode('utf-8-sig')
        reader = csv.DictReader(io.StringIO(text))
        fields = reader.fieldnames
        records = list(reader)
        patterns = {'insertion_torque_English': 'insertion\\s+torque', 'insertion_torque_Chinese': '植入扭矩|植入转矩|植入轉矩|植入扭力', 'drill_diameter_English': '(?:final\\s+drill|drill\\s+diameter|under.?drill)', 'torque_units_any_context': '(?:N\\s*[·.\\-]?\\s*cm|Ncm|牛[顿頓].?厘米)', 'cortex_thickness_English': 'cortical\\s+(?:bone\\s+)?thickness'}
        counts = {key: sum((bool(re.search(pat, ' '.join((str(v) for v in r.values())), re.I)) for r in records)) for (key, pat) in patterns.items()}
        nii = [n for n in names if n.endswith('.nii.gz')]
        mask = [n for n in names if re.search('mask|label|segment', n, re.I)]
        with z.open(nii[0]) as f:
            with gzip.GzipFile(fileobj=f) as g:
                hd = g.read(352)
        pix = list(struct.unpack('<8f', hd[76:108])[1:4])
        scl = list(struct.unpack('<2f', hd[112:120]))
        out['MMDental'] = {'archive_path': str(MM), 'archive_size_bytes': MM.stat().st_size, 'volume_count': len(nii), 'mask_or_segmentation_named_members': len(mask), 'journal_schema_columns': fields, 'journal_rows': len(records), 'journal_sha256': hashlib.sha256(b).hexdigest(), 'phrase_counts_aggregate_only': counts, 'one_volume_spacing_mm': pix, 'one_volume_slope_intercept': scl, 'cortex_inner_boundary_label': False, 'cortex_thickness': 'Images permit new geometric measurement, no local published cortex layer mask; no accuracy guarantee', 'absolute_density': 'UNKNOWN_UNCALIBRATED_GRAY', 'final_bore_profile': 'No structured final-drill field; zero tested English phrase matches is not exhaustive proof of absence', 'thread_core_geometry': 'No structured pitch, depth or core profile field', 'paired_insertion_torque': 'No structured measured insertion-torque field; journal torque mentions cannot be relabeled from prescribed abutment tightening. Prior independent lane audit found no paired insertion outcomes.', 'license': 'CC BY 4.0 per dataset publication, DOI 10.1038/s41597-025-05398-7; no journal text/person-level records exported'}
    out['shared'] = {'resolution_of_new_needed_edge': 'PER_SURFACE_REGION', 'time_scale_of_bore_to_insertion': 'HANDOVER', 'time_scale_of_density_region_to_contact': 'SIMULTANEOUS', 'replacement_measurements_for_PHENOMENOLOGICAL_debt': ['Calibrated local bone density/mechanics', 'Inner/outer cortical surfaces along planned implant contact', 'Measured final osteotomy depth profile', 'Thread/core/flute geometry', 'Paired motor peak and torque-depth record'], 'elapsed_seconds': time.perf_counter() - start, 'volume_arrays_extracted': 0, 'large_arrays_written': 0}
    write('CBCT_AVAILABILITY.json', out)
    print(json.dumps({k: {f: v for (f, v) in x.items() if f in ['training_image_count', 'volume_count', 'journal_rows', 'mask_or_segmentation_named_members', 'phrase_counts_aggregate_only', 'header_spacing_mm']} for (k, x) in out.items()}, indent=2))
if __name__ == '__main__':
    main()
