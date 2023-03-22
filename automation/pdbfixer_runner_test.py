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
    # no example found for this, ie, chiain id start from other than A
    # def test_keep_chainid(self):
    #     pass
    def test_keep_resid(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '1B07.pdb', name='1B07')
        # keep everything
        fixer.run(remove_het=False, keep_water=False)
        self.assertTrue(fixer.success)
        original_residx2name = dict()
        after_residx2name = dict()
        fixer.pdbin_handle.seek(0)
        for i in fixer.pdbin_handle.readlines():
            if i.startswith('ATOM') or i.startswith('HETATM'):
                original_residx2name[i[21:26].strip()] = i[17:20].strip()
        for i in fixer.pdbout_str.split('\n'):
            if i.startswith('ATOM') or i.startswith('HETATM'):
                after_residx2name[i[21:26].strip()] = i[17:20].strip()
        self.assertTrue(len(original_residx2name) <= len(after_residx2name))
        for idx in original_residx2name:
            self.assertEqual(original_residx2name[idx], after_residx2name[idx])
        
    def test_disulfide_bond(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '1IL8.pdb', name='1IL8')
        fixer.run(remove_het=True, keep_water=True)
        self.assertTrue(fixer.success)
        prefixset = set()
        for line in fixer.pdbout_str.split('\n'):
            prefixset.add(line[:6].strip())
        self.assertIn('CONECT', prefixset)
    
    def test_remove_het(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '10GS.pdb', name='10GS')
        fixer.run(remove_het=True, keep_water=True)
        self.assertTrue(fixer.success)
        resname = set()
        for line in fixer.pdbout_str.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resname.add(line[17:20].strip())
        self.assertNotIn('VWW', resname)
        self.assertNotIn('MES', resname)
        self.assertIn('HOH', resname)
    
    def test_keep_het(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '10GS.pdb', name='10GS')
        fixer.run(remove_het=False, keep_water=True)
        self.assertTrue(fixer.success)

        inputresname = dict()
        fixer.pdbin_handle.seek(0)
        for line in fixer.pdbin_handle.readlines():
            if line.startswith('ATOM') or line.startswith('HETATM'):
                if line[17:20].strip() == 'HOH' and 'O   HOH' not in line:
                    continue
                inputresname[line[17:20].strip()] = inputresname.get(line[17:20].strip(), 0) + 1
        resname = dict()
        for line in fixer.pdbout_str.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                if line[17:20].strip() == 'HOH' and 'O   HOH' not in line:
                    continue
                resname[line[17:20].strip()] = resname.get(line[17:20].strip(), 0) + 1
        self.assertIn('VWW', resname)
        self.assertIn('MES', resname)
        self.assertIn('HOH', resname)
        self.assertEqual(inputresname['VWW'], resname['VWW'])
        self.assertEqual(inputresname['MES'], resname['MES'])
        self.assertEqual(inputresname['HOH'], resname['HOH'])

    def test_dummy_keep_water(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '10GS.pdb', name='10GS')
        fixer.run(remove_het=False, keep_water=False)
        self.assertTrue(fixer.success)
        resname = set()
        for line in fixer.pdbout_str.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resname.add(line[17:20].strip())
        self.assertIn('VWW', resname)
        self.assertIn('MES', resname)
        self.assertIn('HOH', resname)

        fixer = PDBFixerRunner(curr_path / 'test_data' / '10GS.pdb', name='10GS')
        fixer.run(remove_het=False, keep_water=True)
        self.assertTrue(fixer.success)
        resname2 = set()
        for line in fixer.pdbout_str.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resname2.add(line[17:20].strip())
        self.assertIn('VWW', resname2)
        self.assertIn('MES', resname2)
        self.assertIn('HOH', resname2)
        self.assertEqual(resname, resname2)



if __name__ == '__main__':
    unittest.main()