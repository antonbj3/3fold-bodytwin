from common import *
import argparse, resource, subprocess

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ratings')
    ap.add_argument('--verify-only', action='store_true')
    args = ap.parse_args()
    from integrity import verify
    tic = time.perf_counter()
    receipt = verify()
    dump(ROOT / 'raw/INPUT_VERIFICATION.json', receipt)
    if args.verify_only:
        print(receipt)
        return
    if args.ratings:
        from lab import analyze
        result = analyze(args.ratings)
        dump(ROOT / 'LAB_RATINGS_ANALYSIS.json', result)
        print(result)
        return
    import item_analysis, boundary, functional, external_support, physics, controls, norms, contrast, lab, x60_table, report
    costs = {}
    for (name, fn) in [('items', item_analysis.run), ('boundary', boundary.run), ('functional', functional.run), ('external_support', external_support.run), ('physical', physics.run), ('controls', controls.run), ('norms', norms.run), ('contrast', contrast.run), ('X60', x60_table.run), ('images', lab.make_pack)]:
        t = time.perf_counter()
        result = fn()
        costs[name] = time.perf_counter() - t
        if name == 'images':
            dump(ROOT / 'raw/LAB_PACKAGE.json', result)
        print(name, 'completed', round(costs[name], 2), 's', flush=True)
    dump(ROOT / 'raw/RUN_COST.json', dict(wall_seconds=time.perf_counter() - tic, stage_seconds=costs, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, threads=4, new_model_fit=False, external_generation='Previously actually executed, frozen output hash-checked; separate source-model replay via replay_external.sh', setup_and_inherited_cost='See source/model receipts; prior agent discovery cost UNKNOWN'))
    result = report.create_report()
    dump(ROOT / 'raw/FINAL_INPUT_VERIFICATION.json', verify())
    state('DEMO_COMPLETE_PARTIAL_PHYSICS', result['headline_results'], 'Independent review; matched specimen and seated-film measurements; expand external pilot after input-alignment audit')
    print(result['status'])
if __name__ == '__main__':
    main()
