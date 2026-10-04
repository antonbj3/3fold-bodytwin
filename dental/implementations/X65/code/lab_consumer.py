"""Read a region-level measured stress port; reject missing empirical calibration."""
import argparse, json
from review_provenance import measured_query

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--R', type=float, required=True)
    a = p.parse_args()
    print(json.dumps(measured_query(json.load(open(a.input)), a.R), indent=2))
if __name__ == '__main__':
    main()
