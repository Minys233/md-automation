import unittest
from pathlib import Path
import logging
from pdb4amber_runner import PDB4AmberRunner
from pdbfixer_runner import PDBFixerRunner

logging.getLogger().setLevel(logging.ERROR)
curr_path = Path(__file__).parent

class TestPDB4AmberRunner(unittest.TestCase):
    def test_pdb4amber_simple_run(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '10GS.pdb', '10GS.pdb', remove_het=False, keep_water=True)
        fixer.run()
        pdbin = fixer.pdbout_str
        pdb4amber = PDB4AmberRunner(pdbin, '10GS.pdb')
        pdb4amber.run()
        self.assertTrue(pdb4amber.success)
    
    def test_pdb4amber_renum(self):
        with open(curr_path / 'test_data' / '10GS.pdb', 'r') as f:
            pdbin = f.read()
        pdb4amber = PDB4AmberRunner(pdbin, '10GS.pdb',)
        pdb4amber.run()
        self.assertTrue(pdb4amber.success)
        self.assertTrue(pdb4amber.renum)

        # missing residue 1 for A
        for name, chain, number, namere, numberre in pdb4amber.renum:
            if chain == 'A' and name not in  ['VWW', 'MES', 'HOH']:
                self.assertEqual(name, namere)
                self.assertEqual(number - numberre, 1)
    
    def test_pdb4amber_rm_water(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '10GS.pdb', '10GS.pdb', remove_het=False, keep_water=True)
        fixer.run()
        pdbin = fixer.pdbout_str
        pdb4amber = PDB4AmberRunner(pdbin, '10GS.pdb', rm_water=True)
        pdb4amber.run()
        self.assertTrue(pdb4amber.success)
        self.assertTrue(pdb4amber.water)
        for line in pdb4amber.pdbout.split('\n'):
            if line.startswith('HETATM'):
                self.assertTrue(line[17:20] != 'HOH')
    
    def test_pdb4amber_prot_only(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '10GS.pdb', '10GS.pdb', remove_het=False, keep_water=True)
        fixer.run()
        pdbin = fixer.pdbout_str
        pdb4amber = PDB4AmberRunner(pdbin, '10GS.pdb', prot_only=True)
        pdb4amber.run()
        self.assertTrue(pdb4amber.success)
        self.assertTrue(not pdb4amber.water)
        for line in pdb4amber.pdbout.split('\n'):
            self.assertTrue(not line.startswith('HETATM'))
        for line in pdb4amber.nonprot.split('\n'):
            self.assertTrue(not line.startswith('ATOM'))

    def test_pdb4amber_fix_and_renum(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '1QG8.pdb', '1QG8.pdb', remove_het=False, keep_water=True)
        fixer.run()
        self.assertTrue(fixer.success)
        fixed_pdbin = fixer.pdbout_str
        with open(curr_path / 'test_data' / '1QG8.pdb', 'r') as f:
            unfixed_pdbin = f.read()
        
        pdb4amber_fixed = PDB4AmberRunner(fixed_pdbin, '1QG8_fixed.pdb', rm_water=True, prot_only=True)
        pdb4amber_fixed.run()
        self.assertTrue(pdb4amber_fixed.success)
        pdb4amber_unfixed = PDB4AmberRunner(unfixed_pdbin, '1QG8_unfixed.pdb', rm_water=True, prot_only=True)
        pdb4amber_unfixed.run()
        self.assertTrue(pdb4amber_unfixed.success)
        renum_fixed = pdb4amber_fixed.renum
        renum_unfixed = pdb4amber_unfixed.renum
        missing_residue_num = [134, 135, 136, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231]
        for missing_residue in missing_residue_num:
            self.assertTrue(missing_residue not in [i[2] for i in renum_unfixed])
            self.assertTrue(missing_residue in [i[2] for i in renum_fixed])
        # for 1QG8, residue 2-256, the missing residues are marked in REMARK 465
        # REMARK 465   M RES C SSSEQI                                                     
        # REMARK 465     GLU A   134                                                      
        # REMARK 465     ASN A   135                                                      
        # REMARK 465     ARG A   136                                                      
        # REMARK 465     ASP A   218                                                      
        # REMARK 465     GLN A   219                                                      
        # REMARK 465     SER A   220                                                      
        # REMARK 465     ILE A   221                                                      
        # REMARK 465     HIS A   222                                                      
        # REMARK 465     PHE A   223                                                      
        # REMARK 465     GLN A   224                                                      
        # REMARK 465     LEU A   225                                                      
        # REMARK 465     PHE A   226                                                      
        # REMARK 465     GLU A   227                                                      
        # REMARK 465     LEU A   228                                                      
        # REMARK 465     GLU A   229                                                      
        # REMARK 465     LYS A   230                                                      
        # REMARK 465     ASN A   231                                                      
    def test_pdb4amber_sslink(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '1IL8.pdb', '1IL8.pdb', remove_het=False, keep_water=True)
        fixer.run()
        self.assertTrue(fixer.success)
        fixed_pdbin = fixer.pdbout_str
        with open(curr_path / 'test_data' / '1IL8.pdb', 'r') as f:
            unfixed_pdbin = f.read()
        pdb4amber_fixed = PDB4AmberRunner(fixed_pdbin, '1IL8_fixed.pdb', rm_water=True, prot_only=True)
        pdb4amber_fixed.run()
        self.assertTrue(pdb4amber_fixed.success)
        pdb4amber_unfixed = PDB4AmberRunner(unfixed_pdbin, '1IL8_unfixed.pdb', rm_water=True, prot_only=True)
        pdb4amber_unfixed.run()
        self.assertTrue(pdb4amber_unfixed.success)

        sslink_fixed = pdb4amber_fixed.sslink
        sslink_unfixed = pdb4amber_unfixed.sslink
        # Missing residue 1 for both chains
        # DBREF  1IL8 A    1    72  UNP    P10145   IL8_HUMAN       28     99             
        # DBREF  1IL8 B    1    72  UNP    P10145   IL8_HUMAN       28     99             
        # SSBOND   1 CYS A    7    CYS A   34                          1555   1555  2.02  
        # SSBOND   2 CYS A    9    CYS A   50                          1555   1555  2.02  
        # SSBOND   3 CYS B    7    CYS B   34                          1555   1555  2.02  
        # SSBOND   4 CYS B    9    CYS B   50                          1555   1555  2.02  
        self.assertTrue((7, 34) in sslink_fixed)
        self.assertTrue((9, 50) in sslink_fixed)
        self.assertTrue((72+7, 72+34) in sslink_fixed)
        self.assertTrue((72+9, 72+50) in sslink_fixed)

        self.assertTrue((6, 33) in sslink_unfixed)
        self.assertTrue((8, 49) in sslink_unfixed)
        self.assertTrue((71+6, 71+33) in sslink_unfixed)
        self.assertTrue((71+8, 71+49) in sslink_unfixed)
    
    def test_pdb4amber_gaplist1(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '1QG8.pdb', '1QG8.pdb', remove_het=False, keep_water=True)
        fixer.run()
        self.assertTrue(fixer.success)
        fixed_pdbin = fixer.pdbout_str
        with open(curr_path / 'test_data' / '1QG8.pdb', 'r') as f:
            unfixed_pdbin = f.read()
        
        pdb4amber_fixed = PDB4AmberRunner(fixed_pdbin, '1QG8_fixed.pdb', rm_water=True, prot_only=True)
        pdb4amber_fixed.run()
        self.assertTrue(pdb4amber_fixed.success)
        pdb4amber_unfixed = PDB4AmberRunner(unfixed_pdbin, '1QG8_unfixed.pdb', rm_water=True, prot_only=True)
        pdb4amber_unfixed.run()
        self.assertTrue(pdb4amber_unfixed.success)
        self.assertTrue(not pdb4amber_fixed.gaplist)
        self.assertTrue(pdb4amber_unfixed.gaplist)
        self.assertTrue(len(pdb4amber_unfixed.gaplist) == 2)
        
    def test_pdb4amber_gaplist2(self):
        fixer = PDBFixerRunner(curr_path / 'test_data' / '10GS.pdb', '10GS.pdb', remove_het=False, keep_water=True)
        fixer.run()
        self.assertTrue(fixer.success)
        fixed_pdbin = fixer.pdbout_str
        with open(curr_path / 'test_data' / '10GS.pdb', 'r') as f:
            unfixed_pdbin = f.read()
        
        pdb4amber_fixed = PDB4AmberRunner(fixed_pdbin, '10GS_fixed.pdb', rm_water=True, prot_only=True)
        pdb4amber_fixed.run()
        self.assertTrue(pdb4amber_fixed.success)
        pdb4amber_unfixed = PDB4AmberRunner(unfixed_pdbin, '10GS_unfixed.pdb', rm_water=True, prot_only=True)
        pdb4amber_unfixed.run()
        self.assertTrue(pdb4amber_unfixed.success)
        # missing residues at starting / ending of chain, no gap
        self.assertTrue(not pdb4amber_fixed.gaplist)
        self.assertTrue(not pdb4amber_unfixed.gaplist)

    def test_pdb4amber_missing_atom_residues(self):
        import random
        random.seed(2022)
        fixer = PDBFixerRunner(curr_path / 'test_data' / '10GS.pdb', '10GS.pdb', remove_het=False, keep_water=True)
        fixer.run()
        self.assertTrue(fixer.success)
        fixed_pdbin = fixer.pdbout_str
        with open(curr_path / 'test_data' / '10GS.pdb', 'r') as f:
            lines = f.readlines()
            f.seek(0)
            atomlineidx = [i for i, _ in enumerate(f) if _.startswith('ATOM')]
            choseidx = random.choices(atomlineidx, k=10)
            unfixed_pdbin = []
            choselines = []
            for idx, line in enumerate(lines):
                if idx in choseidx:
                    choselines.append(line)
                    continue
                unfixed_pdbin.append(line)
            unfixed_pdbin = ''.join(unfixed_pdbin)
            # choselines = ''.join(choselines)
            missing_residue = [l[21:26] for l in choselines]
        # ATOM     24  CB  THR A   4      24.097   3.357  38.132  1.00 18.80           C  
        # ATOM    182  CD  GLN A  24      33.820  -3.268  23.968  1.00 20.22           C  
        # ATOM    238  C   GLU A  31      16.875   0.323  38.411  1.00 23.26           C  
        # ATOM    290  CA  TRP A  38       7.305   6.459  38.203  1.00 39.46           C  
        # ATOM    589  C   THR A  75      36.019   5.137  26.492  1.00 25.60           C  
        # ATOM   1152  CD  GLN A 147      38.307   0.998  15.731  1.00 55.59           C  
        # ATOM   1431  CD1 LEU A 183      27.252  -6.391  15.057  1.00 24.76           C  
        # ATOM   1875  OE1 GLU B  31      25.642  23.946  -6.515  1.00 63.47           O  
        # ATOM   2782  CB  GLN B 147      34.542  27.154  24.896  1.00 35.52           C  
        # ATOM   2841  CB  ASN B 154      21.522  22.838  17.535  1.00  9.31           C
        pdb4amber_fixed = PDB4AmberRunner(fixed_pdbin, '10GS_fixed.pdb', rm_water=True, prot_only=True)
        pdb4amber_fixed.run()
        pdb4amber_unfixed = PDB4AmberRunner(unfixed_pdbin, '10GS_unfixed.pdb', rm_water=True, prot_only=True)
        pdb4amber_unfixed.run()
        self.assertTrue(pdb4amber_fixed.success)
        self.assertTrue(not pdb4amber_fixed.missing_atom_residues)
        self.assertTrue(pdb4amber_unfixed.success)
        self.assertTrue(pdb4amber_unfixed.missing_atom_residues)
        # [('THR', 3, 1), ('GLN', 23, 1), ('GLU', 30, 1), ('TRP', 37, 1), ('THR', 74, 1), ('GLN', 146, 1), ('LEU', 182, 1), ('GLU', 238, 1), ('GLN', 354, 1), ('ASN', 361, 1)]
        self.assertEqual(sum([i[2] for i in pdb4amber_unfixed.missing_atom_residues]), len(missing_residue))

if __name__ == '__main__':
    unittest.main()