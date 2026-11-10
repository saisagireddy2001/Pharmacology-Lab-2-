import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from analyze import response_percent, inhibition_percent, log_molar, analyze_csv, make_suxa_table, ROOT

class TestAnalysis(unittest.TestCase):
    def test_percentages(self):
        self.assertAlmostEqual(response_percent(49, 100),49)
        self.assertAlmostEqual(inhibition_percent(49,100),51)
    def test_invalid_baseline(self):
        with self.assertRaises(ValueError):response_percent(1,0)
    def test_log_molar(self):
        self.assertAlmostEqual(log_molar(1),-6)
    def test_eight_dose_response(self):
        _,summary=analyze_csv(ROOT/'data/simulated_csv/Part_III_Suxamethonium_Eight_Dose_Curve.csv')
        table=make_suxa_table(summary)
        self.assertEqual(len(table),8)
        self.assertEqual(table.concentration_uM.iloc[-1],100)

if __name__ == '__main__':unittest.main()
