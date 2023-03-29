from pdbfixer import PDBFixer
from pdbfixer.pdbfixer import proteinResidues, dnaResidues, rnaResidues
from openmm.app import PDBFile, ForceField
import logging
import io
from typing import Union, List
from pathlib import Path
from subprocess_runner import Runner
from util import StringStream, PDB_TITLE_PREFIX, PDB_MAINBODY_PREFIX, PDB_OTHER_PREFIX

logging.basicConfig()
logger = logging.getLogger(__name__)


# def get_temperature_pH(remarks_lines, default_temperature=300, default_pH=7.0):
#     temp, pH = None, None
#     for i in remarks_lines:
#         if i.startswith('REMARK 200  PH'):
#             try:
#                 pH = float(i.split(':')[-1])
#             except Exception as e:
#                 logger.warning(f"Could not parse pH from {i.split(':')[-1]}, setting to default {default_pH}")
#         if i.startswith('REMARK 200  TEMPERATURE'):
#             try:
#                 temp = float(i.split(':')[-1])
#             except Exception as e:
#                 logger.warning(f"Could not parse temperature from {i.split(':')[-1]}, setting to default {default_temperature}")
#     if not temp:
#         temp = default_temperature
#     if not pH:
#         pH = default_pH
#     return temp, pH

class PDBFixerRunner(Runner):
    """A subclass of PDBFixer that can be used to fix PDB files.

    This subclass adds methods to fix PDB files and write them to a new file.
    """

    def __init__(self, input: Union[str, Path, io.IOBase, StringStream], remove_het: bool=False, keep_water: bool=True, replace_nonstandard: bool=False, addH: bool=False, pH: float=7.0, mutations: List[str] = []):
        self.pdbin = StringStream(input, name='input').read()
        # parameters for pdbfixer
        self.remove_het = remove_het
        self.keep_water = keep_water
        self.replace_nonstandard = replace_nonstandard
        self.addH = addH
        self.pH = pH
        self.mutations = mutations
        self.pdbout = ""

        cmd = f"pdbfixer" + " --keep-heterogens" if not remove_het else "" + \
              " water" if keep_water else " all" + \
              f" --add-residues --add-atoms --replace-nonstandard | fixer.applyMutations({mutations})"
        super().__init__(f'Equivalent to run: {cmd}', False)

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name})'

    def run(self, ignore_opreated: bool=False) -> None:
        """Run the PDBFixer pipeline, results are stored in self.pdbout as str.
        """
        if self._operated:
            if ignore_opreated:
                logging.warning(f'{self} has already been operated, running again.')
            else:
                logging.error(f'{self} has already been operated, skipping.')
                return
        self._operated = True

        fixer = PDBFixer(pdbfile=io.StringIO(self.pdbin))
        for mutstr in self.mutations:
            logger.warning(f"Applying mutation: {mutstr}")
            chainid, mut = mutstr.split(':')
            fixer.applyMutations(mutations=[mut], chain_id=chainid)
        
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()

        if self.replace_nonstandard:
            fixer.replaceNonstandardResidues()
            logger.warning(f"Replacing nonstandard residues.")

        logger.warning(f"Adding missing residues.")
        if self.remove_het:
            logger.warning(f"Removing heterogens.")
            fixer.removeHeterogens(keepWater=self.keep_water)
        else:
            logger.warning(f"Not removing any heterogens.")
        fixer.findMissingAtoms()
        logger.warning(f"Adding missing atoms.")
        logger.warning(f"Adding missing terminal atoms.")
        fixer.addMissingAtoms(seed=None)
        if self.addH:
            fixer.addMissingHydrogens(pH=self.pH, seed=None)
            logger.warning(f"Adding missing hydrogens.")
        pdbout = io.StringIO()
        PDBFile.writeFile(fixer.topology, fixer.positions, pdbout, keepIds=True,)
        self.pdbout = pdbout.getvalue()

    @property
    def success(self) -> bool:
        return self.pdbout != ""

if __name__ == '__main__':
    fixer = PDBFixerRunner('./test_data/10GS.pdb', '10GS.pdb', remove_het=False, keep_water=True)
    fixer.run()
    fixer.write('./test_data/10GS.fixed.pdb')

