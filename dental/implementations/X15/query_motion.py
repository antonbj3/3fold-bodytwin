import argparse, json
from pathlib import Path
from motion_envelope import query
p = argparse.ArgumentParser(description='Research section-plane geometry; physical status always UNDETERMINED')
p.add_argument('--case', default='ToothFairy2F_001')
p.add_argument('--tooth', type=int, default=44)
p.add_argument('--source', choices=['label', 'image'], default='label')
p.add_argument('--translation-mm', type=float, default=0.5)
p.add_argument('--angle-deg', type=float, default=5)
p.add_argument('--epsilon-mm', type=float, default=0.608)
a = p.parse_args()
models = json.load(open(Path(__file__).with_name('MOTION_MODELS.json')))
model = next((m for m in models if (m['case'], m['tooth'], m['source']) == (a.case, a.tooth, a.source)), None)
print(json.dumps(query(model, a.translation_mm, a.angle_deg, a.epsilon_mm) if model else {'anatomical_status': 'UNDETERMINED', 'reason': 'PROFILE_NOT_COMPLETE'}, indent=2))
