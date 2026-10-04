"""Manual primary table cell transcription is independent of extraction paths."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from literature import extract
EXPECTED = {('PMC10314363', '3D printed (SLM)', 'Mean'): 27.2, ('PMC10314363', 'CAD/CAM milling', 'Mean'): 44.4, ('PMC10333096', 'EZIS', 'Aegis HM'): 39.68, ('PMC10333096', 'EZIS', 'Trione Z'): 39.2, ('PMC10333096', 'EZIS', 'Motion 2'): 36.19, ('PMC10333096', '3Shape', 'Aegis HM'): 43.78, ('PMC10333096', '3Shape', 'Trione Z'): 42.0, ('PMC10333096', '3Shape', 'Motion 2'): 37.15, ('PMC10333096', 'Exocad', 'Aegis HM'): 41.73, ('PMC10333096', 'Exocad', 'Trione Z'): 37.67, ('PMC10333096', 'Exocad', 'Motion 2'): 36.79}

def primary_gate(cells):
    actual = {(c['study'], c['row'], c['column']): c['mean_um'] for c in cells if c['object'] in ['crown', 'crown_coping']}
    return actual.keys() == EXPECTED.keys() and max((abs(actual[k] - v) for (k, v) in EXPECTED.items())) <= 0.01

class SourceTests(unittest.TestCase):

    def test_all_cells_and_corruption(self):
        (cells, _) = extract()
        self.assertTrue(primary_gate(cells))
        broken = [dict(c) for c in cells]
        broken[0]['mean_um'] += 1
        self.assertFalse(primary_gate(broken))
if __name__ == '__main__':
    unittest.main()
