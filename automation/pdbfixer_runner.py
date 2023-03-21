from pdbfixer import PDBFixer
from pdbfixer.pdbfixer import proteinResidues, dnaResidues, rnaResidues
from openmm.app import PDBFile, ForceField
import logging
import io
from typing import Union, List
from pathlib import Path
from subprocess_runner import Runner
from util import StringStream

logging.basicConfig()
logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)


def get_temperature_pH(remarks_lines, default_temperature=300, default_pH=7.0):
    temp, pH = None, None
    for i in remarks_lines:
        if i.startswith('REMARK 200  PH'):
            try:
                pH = float(i.split(':')[-1])
            except Exception as e:
                logger.warning(f"Could not parse pH from {i.split(':')[-1]}, setting to default {default_pH}")
        if i.startswith('REMARK 200  TEMPERATURE'):
            try:
                temp = float(i.split(':')[-1])
            except Exception as e:
                logger.warning(f"Could not parse temperature from {i.split(':')[-1]}, setting to default {default_temperature}")
    if not temp:
        temp = default_temperature
    if not pH:
        pH = default_pH
    return temp, pH


class PDBFixerRunner(Runner):
    """A subclass of PDBFixer that can be used to fix PDB files.

    This subclass adds methods to fix PDB files and write them to a new file.
    """
    # https://www.wwpdb.org/documentation/file-format-content/format33
    _TITLE_PREFIX = (
        # sect2.html Title section
        "HEADER", "OBSLTE", "TITLE", "SPLIT", "CAVEAT", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "REVDAT", "SPRSDE", "JRNL", "REMARK"
    )
    _MAINBODY_PREFIX = (
        # sect3.html Primary structure section
        "SEQRES", "DBREF", # "DBREF1", "DBREF2", "SEQADV", "MODRES", #  no longer used
        # sect4.html Heterogen section
        "HET", "HETNAM", "HETSYN", "FORMUL", # no longer used
        # sect5.html Secondary structure section
        # "HELIX", "SHEET", "TURN", # no longer used
        # sect8.html Crystallographic and Coordinate Transformation Section
        "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3", "MTRIX1", "MTRIX2", "MTRIX3",
        # sect9.html Coordinate section
        "MODEL", "ATOM", "TER", "HETATM", "ENDMDL", # "ANISOU", # no longer used
        # sect10.html Connectivity section
        "CONECT",
        # sect11.html Bookkeeping section
        "END", # "MASTER", # no longer used
    )
    _OTHER_PREFIX = (
        # sect6.html The connectivity annotation section
        "SSBOND", "LINK", "CISPEP",
        # sect7.html Miscellaneous Features Section
        "SITE",
    )

    def __init__(self, input: Union[str, Path, io.IOBase], name: str = '', default_temperature:float=300.15, default_pH: float=7.0):
        self.pdbin_handle = StringStream(input, name='input')
        self.name = name
        self.title_lines = []
        self.mainbody_lines = []
        self.other_lines = []
        self.notused_lines = []
        self.temp, self.pH = get_temperature_pH(self.title_lines, default_temperature=default_temperature, default_pH=default_pH)
        logger.info(f'Found in {name} - Temperature: {self.temp}, pH: {self.pH}')
        self._forcefield = None # ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        self._parse()
        self._operated = False
        self._cmd = "Not yet run"
        self.pdbout_str = ""

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name})'
    
    def _parse(self):
        """Parse the PDB file."""
        for line in self.pdbin_handle:
            if line.startswith(self._TITLE_PREFIX):
                self.title_lines.append(line)
            elif line.startswith(self._MAINBODY_PREFIX):
                self.mainbody_lines.append(line)
            elif line.startswith(self._OTHER_PREFIX):
                self.other_lines.append(line)
            else:
                logger.info(f'Skipping line with unused prefix: {line.strip()}')
                self.notused_lines.append(line)
        self.pdbin_handle.seek(0)

    def _het2delete(self, fixer, keep_water: bool) -> list:
        """calculate the residues to delete, for logging."""
        keep = set(proteinResidues).union(dnaResidues).union(rnaResidues)
        keep.add('N')
        keep.add('UNK')
        if keep_water:
            keep.add('HOH')
        removed_heterogens = []
        for residue in fixer.topology.residues():
            if residue.name not in keep:
                removed_heterogens.append(residue)
        return removed_heterogens

    def run(self, ignore_opreated: bool=False, remove_het: bool=False, keep_water: bool = True, mutations: List[str] = []) -> None:
        """Run the PDBFixer pipeline, results are stored in self.pdbout_str as str.

        :param remove_het: remove heterogens, don't do this for leap will recognize some of them, defaults to False
        :param keep_water: keep water in the structure, defaults to True
        :type keep_water: bool, optional
        :param mutations: _description_, defaults to []
        :type mutations: List[str], optional
        """
        if self._operated:
            if ignore_opreated:
                logging.warning(f'{self} has already been operated, running again.')
            else:
                logging.error(f'{self} has already been operated, skipping.')
                return
        self._operated = True
        self._cmd = f"remove_het={remove_het}, keep_water={keep_water}, mutations={mutations}"

        try:
            pdbstrin = ''.join(self.title_lines + self.mainbody_lines + self.other_lines)
            fixer = PDBFixer(pdbfile=io.StringIO(pdbstrin))
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            if remove_het:
                removed_heterogens = self._het2delete(fixer, keep_water=keep_water)
                fixer.removeHeterogens(keepWater=keep_water)

            logger.warning(f"Missing residues: {fixer.missingResidues}")
            logger.warning(f"Nonstandard residues: {fixer.nonstandardResidues}")
            logger.warning(f"Removed heterogens? {remove_het}: {removed_heterogens}")
            fixer.findMissingAtoms()
            fixer.addMissingAtoms(seed=None)
            fixer.addMissingHydrogens(self.pH, forcefield=self._forcefield)

            for mutstr in mutations:
                chainid, mut = mutstr.split(':')
                fixer.applyMutations(mutations=[mut], chain_id=chainid)

            pdbout = io.StringIO()
            PDBFile.writeFile(fixer.topology, fixer.positions, pdbout, keepIds=True)

            pdbout.seek(0)
            lines = pdbout.readlines()
            assert lines[1].startswith('CRYST1'), f"Popped a wrong line: {lines[1]}"
            lines[1] = f"REMARK   1 remove_het={remove_het} keep_water={keep_water} mutations={mutations}\n"
            lines = lines[:2] + self.title_lines + [l for l in self.mainbody_lines if l.startswith(self._MAINBODY_PREFIX[:12])] + self.other_lines + lines[2:]
            self.pdbout_str = ''.join(lines)
        except Exception as e:
            logger.critical(f'Error in PDBFixer: {e}')
            self.pdbout_str = ""

    def write(self, outfile: Union[str, Path]):
        """Write the fixed PDB file to a file.

        :param outfile: the output file
        :type outfile: Union[str, Path]
        """
        if self._operated:
            with open(outfile, 'w') as f:
                f.write(self.pdbout_str)
        else:
            logger.error(f'{self} has not been operated, skipping.')

    @property
    def operated(self) -> bool:
        return self._operated

    @property
    def cmd(self) -> str:
        return self._cmd
    
    @property
    def success(self) -> bool:
        return self.pdbout_str != ""

if __name__ == '__main__':
    fixer = PDBFixerRunner('test_data/1PLC.pdb', name='1PLC')
    print("====== TITLE lines ======")
    print(len(fixer.title_lines))
    print("====== MAINBODY lines ======")
    print(len(fixer.mainbody_lines))
    print("====== OTHER lines ======")
    print(len(fixer.other_lines))
    fixer.run(remove_het=True, keep_water=True) #mutations=['A:PRO-2-GLY', 'B:PRO-2-GLY'])
    # fixer.write('test_data/10GS_fixed.pdb')