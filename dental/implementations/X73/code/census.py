from common import *

def main():
    check_frozen()
    start = time.perf_counter()
    records = []
    metadata = []
    with zipfile.ZipFile(TF1) as a, zipfile.ZipFile(MAX) as b, zipfile.ZipFile(TF2) as c:
        indexes = {}
        for (key, z) in [('tf1', a), ('maxillo', b), ('tf2', c)]:
            idx = [dict(member=i.filename, crc=i.CRC, bytes=i.file_size, compressed_bytes=i.compress_size, method=i.compress_type) for i in z.infolist() if not i.is_dir()]
            indexes[key] = idx
            for i in idx:
                n = i['member']
                if n.endswith(('README.md', 'dataset.json', 'masks.json')):
                    data = member(z, n)
                    path = ROOT / 'sources' / key / n
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(data)
                    obj = json.loads(data) if n.endswith('.json') else None
                    metadata.append(dict(dataset=key, member=n, sha256=sha(path), bytes=len(data), json_root_keys=list(obj) if isinstance(obj, dict) else None, independent_reader_attribution=False))

        def ix(z, s):
            return {i.filename.split('/')[-2]: i.filename for i in z.infolist() if i.filename.endswith(s)}
        ai = ix(a, 'data.npy')
        ad = ix(a, 'gt_alpha.npy')
        asp = ix(a, 'gt_sparse.npy')
        bi = ix(b, 'data.npy')
        bd = ix(b, 'gt_alpha.npy')
        bsp = ix(b, 'gt_sparse.npy')
        labels = {n.split('/')[-1][:-4]: n for n in c.namelist() if '/labelsTr/' in n and n.endswith('.mha')}
        for p in sorted(set(ai) | set(bi), key=lambda x: int(x[1:])):
            case = f'ToothFairy2P_{int(p[1:]):03d}'
            walls = []
            sparse = []
            for (k, m) in [('tf1', ad), ('maxillo', bd)]:
                if p in m:
                    walls.append(dict(dataset=k, member=m[p], status='PROTOCOL_REVISION_OR_REUSED_CONSENSUS', independent=False))
            if case in labels:
                walls.append(dict(dataset='tf2', member=labels[case], status='ONE_CONSENSUS_VOLUME; PER_READER_LINEAGE_UNKNOWN', independent=False))
            for (k, m) in [('tf1', asp), ('maxillo', bsp)]:
                if p in m:
                    sparse.append(dict(dataset=k, member=m[p], status='DIFFERENT_OBSERVABLE_CENTERLINE_NOT_WALL', independent=False))
            records.append(dict(patient=p, case=case, walls=walls, sparse=sparse, dense_pair_candidate=len(walls) >= 2, independent_pair_admitted=False, rejection_reason='no per-reader independent wall-mask attribution', original_images_crc_equal=a.getinfo(ai[p]).CRC == b.getinfo(bi[p]).CRC and a.getinfo(ai[p]).file_size == b.getinfo(bi[p]).file_size if p in ai and p in bi else None))
        dump(ROOT / 'raw/ARCHIVE_INDEX.json', indexes)
        dump(ROOT / 'raw/METADATA_AUDIT.json', metadata)
        dump(ROOT / 'raw/PATIENT_CENSUS.json', records)
        result = dict(claim_type='information_link', wall_patients_tf1=len(ad), wall_patients_maxillo=len(bd), wall_volumes_tf2=len(labels), tf1_image_patients=len(ai), maxillo_image_patients=len(bi), tf1_sparse_patients=len(asp), maxillo_sparse_patients=len(bsp), dense_pair_candidates=sum((r['dense_pair_candidate'] for r in records)), all_multi_annotation_candidates=sum((len(r['walls']) + len(r['sparse']) >= 2 for r in records)), metadata_files_checked=len(metadata), independent_pairs_admitted=0, independent_pair_rejection_fraction=1.0, published_independent_substudy=dict(locator='https://www.federicobolelli.it/media/publications/pdfs/2023iciap_iacat2.pdf section4.1', patients=6, available_per_reader_masks=False), gate='UNKNOWN_INDEPENDENT_SPREAD_NOT_IDENTIFIABLE', cost=dict(wall_seconds=time.perf_counter() - start, threads=1, gpu=False))
        dump(ROOT / 'raw/R1_RESULT.json', result)
        (ROOT / 'HANDOFF_R1.md').write_text("# R1 : independent wall variation missing identifiable pairs\n\nThe entire local membership and annotation metadata has been read. Clear labels are midlines, frequent labels audits/consensus. The published double annotation of six patients has no identified separate reader masks here. See raw/R1_RESULT.json for denominator .\n\nThe next design R2 retains support and direction of all frequent audit pairs; it does not call them independent annotator variation.\n")
        state('R1_COMPLETE', result['gate'], 'R2 stream all dense revisions and retain coverage loss', dense_candidates=result['dense_pair_candidates'])
        print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
