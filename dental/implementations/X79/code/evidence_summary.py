from common import *
from data import reports, load_geometry, treatment_rows
import collections, zipfile, re

def main():
    st = time.perf_counter()
    cpu = time.process_time()
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    g = load_geometry()
    ids = m['train'] + m['calibration'] + m['test']
    summaries = {}
    rr_by_modality = {}
    for modality in ['ios', 'intraoral-photo']:
        rr = reports(ids, modality)
        rr_by_modality[modality] = rr
        ars = [a for reps in rr.values() for r in reps for a in r['parsed']['assertions']]
        ex = [a for reps in rr.values() for r in reps for a in r['parsed']['excluded']]
        summaries[modality] = dict(reports=sum((len(v) for v in rr.values())), report_patients=sum((bool(v) for v in rr.values())), accepted_assertion_clauses=len(ars), excluded_clause_events=len(ex), candidate_clause_event_dropout_fraction=len(ex) / (len(ex) + len(ars)) if ex or ars else None, exclusion_reasons=dict(collections.Counter((r['reason'] for r in ex))), positive_named_tooth_occurrences=sum((len(a['teeth']) for a in ars if a['polarity'])), positive_named_patients=sum((any((a['polarity'] and a['teeth'] for r in reps for a in r['parsed']['assertions'])) for reps in rr.values())), crown_named_tooth_occurrences=sum((len(a['teeth']) for a in ars if a['polarity'] and a['kind'] == 'crown')), global_explicit_negative_reports=sum((any((a['scope'] == 'ALL_REPORTED_TEETH' for a in r['parsed']['assertions'])) for reps in rr.values() for r in reps)), hedged_or_unlocalized_are_unknown=True)
    photo = rr_by_modality['intraoral-photo']
    ios = rr_by_modality['ios']
    paired = []
    new = []
    for c in ids:
        if not photo[c]:
            continue
        (prows, _) = treatment_rows({c: [photo[c][0]]}, g)
        ipos = {f for a in ios[c][0]['parsed']['assertions'] if a['polarity'] for f in a['teeth']} if ios[c] else set()
        for r in prows:
            if r['y'] and r['fdi'] not in ipos:
                new.append(r)
        if len(photo[c]) >= 2:
            (r1, _) = treatment_rows({c: [photo[c][0]]}, g)
            (r2, _) = treatment_rows({c: [photo[c][1]]}, g)
            a = {r['fdi']: r['y'] for r in r1}
            b = {r['fdi']: r['y'] for r in r2}
            for f in sorted(set(a) & set(b)):
                paired.append(dict(case_id=c, fdi=f, report1=a[f], report2=b[f], member1=photo[c][0]['member'], member2=photo[c][1]['member'], sha256_1=photo[c][0]['sha256'], sha256_2=photo[c][1]['sha256']))
    repeat = dict(paired_teeth=len(paired), paired_patients=len({r['case_id'] for r in paired}), discordant_teeth=sum((r['report1'] != r['report2'] for r in paired)), rater_noise_floor='UNKNOWN; author, repeats, translation and timing provenance not supplied; report consistency only')
    wear = []
    with zipfile.ZipFile(ZIP) as z:
        for n in sorted(z.namelist()):
            if not n.endswith('.txt') or not any(('/reports_' + mod + '_en/' in n for mod in ['ios', 'intraoral-photo'])):
                continue
            b = z.read(n)
            txt = b.decode()
            for s in re.split('(?<=[.!?])\\s+', txt):
                if re.search('\\bwear\\b|\\bworn\\b|\\battrition\\b|\\babrasion\\b|\\babfraction\\b', s, re.I):
                    teeth = sorted(set((int(f) for f in re.findall('\\b([1-8][1-8])\\b', s))))
                    subtype = 'cervical_abfraction_mention' if re.search('abfraction', s, re.I) else 'abrasion_mention' if re.search('abrasion', s, re.I) else 'primary_incisal_wear_mention' if re.search('primary', s, re.I) else 'incisal_wear_mention' if re.search('incis', s, re.I) else 'unlocalized_wear_mention'
                    wear.append(dict(case_id=n.split('/')[0], member=n, sha256=digest(b), clause=s, named_teeth=teeth, reported_subtype=subtype, uncertain=bool(re.search('seem|appear|may|cannot', s, re.I)), resolution='PER_TOOTH' if teeth else 'PER_ARCH', timescale='SIMULTANEOUS', mechanism='UNKNOWN; reported wear/abrasion/abfraction is not a measured wear rate or inferred force'))
    write('raw/RESTORATION_REPEAT_PAIRS.json', paired)
    write('raw/NEW_PHOTO_RESTORATION_LINKS.json', new)
    write('raw/WEAR_REPORT_SOURCE_CLAUSES.json', wear)
    summary = dict(claim_type='information_link', external_referent=read('PREREG_R1_RESTORATION_GEOMETRY.json')['external_referent'], report_observability=summaries, new_first_photo_teeth_not_in_first_IOS=len(new), new_patients=len({r['case_id'] for r in new}), repeated_photo_observations=repeat, wear=dict(clauses=len(wear), patients=len({r['case_id'] for r in wear}), test_patients=len({r['case_id'] for r in wear if r['case_id'] in m['test']}), outcome='UNKNOWN_INSUFFICIENT_SOURCE_AND_NO_UNWORN_REFERENCE', specificity='UNKNOWN; silence is not wear absence', measurement_needed='Complete scored wear exam plus baseline/unworn shape or longitudinal IOS. Color/clinical history may distinguish restoration/erosion.'), cost=cost(st, cpu))
    write('raw/REPORT_OBSERVABILITY.json', summary)
    inherited = {}
    for (label, path) in [('crowding', X76 / 'raw/R2_RESULTS.json'), ('crossbite_scissor', DENT / 'results/LANE_X50_BITE_SCREENING/results.json')]:
        r = read(path)
        inherited[label] = dict(path=str(path), sha256=sha(path), result=r, own_new_fit=False)
    write('raw/INHERITED_FINDINGS.json', inherited)
    print('New photo-only first-report teeth', len(new), 'patients', len({r['case_id'] for r in new}), 'repeat', repeat, flush=True)
if __name__ == '__main__':
    main()
