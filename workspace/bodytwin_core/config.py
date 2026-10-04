from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]

def paths(config=None):
    cfg = json.loads(Path(config or Path(__file__).with_name('config.json')).read_text())
    return {k: (Path(v) if Path(v).is_absolute() else ROOT / v) for k, v in cfg.items()}
