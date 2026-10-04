from score_ray import *

def run():
    p = R / 'FROZEN_PREDICTIONS_D_STREAM.jsonl'
    assert not p.exists()
    rows = []
    seen = set()
    previous = '0' * 64
    st = time.perf_counter()
    while True:
        records = read(R / 'raw/D_GENERATION.json') if (R / 'raw/D_GENERATION.json').exists() else []
        for rec in records:
            if rec['key'] in seen:
                continue
            record = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), previous_record_sha256=previous, prereg_sha256=sha(R / 'PREREG_D.json'), generator_sha256=sha(R / 'code/oriented_crown.py'), row=rec)
            raw = json.dumps(record, sort_keys=True)
            previous = hashlib.sha256(raw.encode()).hexdigest()
            with p.open('a') as h:
                h.write(raw + '\n')
                h.flush()
                os.fsync(h.fileno())
            if rec['status'] == 'GENERATED':
                try:
                    out = evaluate_c(rec)
                except Exception as e:
                    out = dict(**rec, score_error=repr(e), geometric_conjunction=False)
            else:
                out = dict(**rec, geometric_conjunction=False)
            out['individual_prediction_sha256'] = previous
            rows.append(out)
            seen.add(rec['key'])
            dump(R / 'raw/D_SCORE.json', rows)
            dump(R / 'RESULTS_D_PARTIAL.json', dict(claim_type='capability', rows=rows, external_referent=read(R / 'PREREG_D.json')['external_referent']))
            print(out['key'], out.get('gates'), out.get('score_error'), flush=True)
        if (R / 'FROZEN_PREDICTIONS_D.json').exists() and len(rows) == len(read(R / 'FROZEN_PREDICTIONS_D.json')['rows']):
            break
        time.sleep(3)
    p.with_suffix('.sha256').write_text(sha(p) + '  ' + p.name + '\n')
    dump(R / 'RESULTS_D.json', dict(claim_type='capability', rows=rows, seconds_including_wait=time.perf_counter() - st, external_referent=read(R / 'PREREG_D.json')['external_referent']))
if __name__ == '__main__':
    run()
