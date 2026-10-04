"""Query exported fields or a future independently registered same-subject pair."""
import argparse, json
from region_field import RegionField, registered_query
from common import load, clean

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field', required=True, help='Source field descriptor JSON')
    parser.add_argument('--points', required=True, help='JSON file containing an N x 3 coordinate list, mm')
    parser.add_argument('--patient', required=True)
    parser.add_argument('--frame', required=True)
    parser.add_argument('--pose-mm', type=float, default=0)
    parser.add_argument('--bridge')
    parser.add_argument('--target-field')
    args = parser.parse_args()
    field = RegionField.open(args.field)
    points = load(args.points)
    field.query(points, patient_id=args.patient, frame_id=args.frame, theta=dict(pose_delta_mm=args.pose_mm))
    if args.bridge:
        if not args.target_field:
            parser.error('--target-field required with --bridge')
        result = registered_query(load(args.bridge), field, RegionField.open(args.target_field), points)
    else:
        result = field.query(points, patient_id=args.patient, frame_id=args.frame, theta=dict(pose_delta_mm=args.pose_mm))
    print(json.dumps(clean(result), indent=2, ensure_ascii=False, allow_nan=False))
if __name__ == '__main__':
    main()
