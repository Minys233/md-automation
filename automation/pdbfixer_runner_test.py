import unittest
from pdbfixer_runner import PDBFixerRunner
from pathlib import Path
import logging

logging.getLogger().setLevel(logging.ERROR)

curr_path = Path(__file__).parent


class TestPDBFixerRunner(unittest.TestCase):
    def test_parse(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '10GS.pdb', name='10GS')
        self.assertEqual(len(fixer.title_lines), 313)
        self.assertEqual(len(fixer.mainbody_lines), 3670)
        self.assertEqual(len(fixer.other_lines), 14)
        with open(curr_path / 'test_data' / '10GS.pdb') as f:
            self.assertEqual(len(f.readlines()), len(fixer.title_lines) + len(fixer.mainbody_lines) + len(fixer.other_lines) + len(fixer.notused_lines))
        fixer = PDBFixerRunner(curr_path / 'test_data' / '1PLC.pdb', name='1PLC')
        self.assertEqual(len(fixer.title_lines), 349)
        self.assertEqual(len(fixer.mainbody_lines), 1683)
        self.assertEqual(len(fixer.other_lines), 6)
        with open(curr_path / 'test_data' / '1PLC.pdb') as f:
            self.assertEqual(len(f.readlines()), len(fixer.title_lines) + len(fixer.mainbody_lines) + len(fixer.other_lines) + len(fixer.notused_lines))

    def test_run(self):
        pass

if __name__ == '__main__':
    unittest.main()