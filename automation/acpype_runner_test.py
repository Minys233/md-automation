import unittest
from pdbfixer_runner import PDBFixerRunner
from pdb4amber_runner import PDB4AmberRunner
from acpype_runner import AcpypeRunner, add_H4pH_ligand_obabel
from util import ligand_from_rcsb_model, ligand_from_ligand_expo, find_ligand
from pathlib import Path
from typing import Union
from rdkit import Chem

import logging


logging.basicConfig()
logger = logging.getLogger(__name__)
curr_path = Path(__file__).parent


class TestAcepypeRunner(unittest.TestCase):
    def test_fetch_check_small_molecule_ligand(self):
        pdb = Path(curr_path / 'test_data' / '10GS.pdb')
        fixer = PDBFixerRunner(pdb, '10GS', remove_het=False, keep_water=True,)
        fixer.run()

        pdb4amber = PDB4AmberRunner(fixer.pdbout, name='10GS_fix', rm_water=False, prot_only=False, amber_only=False)
        pdb4amber.run()

        ligands = find_ligand(pdb4amber.pdbout, 10)
        for idx, lig in enumerate(ligands):
            if len(lig.residues) == 1:
                ligand_str_rcsb_model = ligand_from_rcsb_model('10GS', lig.residues[0].name, lig.residues[0].chain, 'mol2')
                ligand_str_ligand_expo = ligand_from_ligand_expo('10GS', lig.residues[0].name, lig.residues[0].chain)
                mol2 = Chem.MolFromMol2Block(ligand_str_rcsb_model)
                mol = Chem.MolFromMolBlock(ligand_str_ligand_expo)
                self.assertTrue(mol is not None)
                self.assertTrue(mol2 is not None)
                mol2_adj = add_H4pH_ligand_obabel(ligand_str_rcsb_model, 7.0, 'mol2')
                mol_adj = add_H4pH_ligand_obabel(ligand_str_ligand_expo, 7.0, 'sdf')
                self.assertTrue(mol_adj)
                self.assertTrue(mol2_adj)

    def test_acpype_runner(self):
        molin = ligand_from_rcsb_model('10GS', 'VWW', 'A', format='mol')
        acpype_runner = AcpypeRunner(molin, name='VWW', pH=7.0, charge_type='gas',
                                    atom_type='gaff2', result_dir='/home/yaosen/10GS_ligand_A_VWW')
        acpype_runner.run()
        self.assertTrue(acpype_runner.success)


if __name__ == '__main__':
    unittest.main()