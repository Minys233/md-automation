from runner import Runner
from subprocess_runner import SubprocessRunner
from typing import Union
from pathlib import Path
import io
from util import StringStream
import pdb4amber as _pdb4amber
import logging
import parmed
logging.basicConfig()
logger = logging.getLogger(__name__)

class PDB4AmberRunner(Runner):
    def __init__(self, pdbin: Union[str, Path, io.IOBase, StringStream], name: str = '', **kwargs):
        self.pdbin = StringStream(pdbin, name)
        self.kwargs = kwargs

        super().__init__(f'pdb4amber', False)
        # result stored here
        self.renum = None
        self.sslink = None
        self.gaplist = None
        self.missing_atom_residues = None
        self.pdbout = None

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.pdbin.name}), rm_water={self.rm_water}, prot_only={self.prot_only}, amber_only={self.amber_only}, custom_mask={self.custom_mask})'

    def run(self, ignore_opreated: bool=False):
        """This function is modified from pdb4amber.AmberPDBFixer.run() funtion.
        Most of the parameters are fixed as we has an already fixed pdb file.
        """
        if self._operated:
            if ignore_opreated:
                logging.warning(f'{self} has already been operated, running again.')
            else:
                logging.error(f'{self} has already been operated, skipping.')
                return
        self._operated = True

        logger.info(f"Using pdb4amber package with version: {_pdb4amber.__version__}")
        logger.info(f"Running pdb4amber for: {self.pdbin.name}")
        pdbin = self.pdbin._handle
        parm = parmed.read_PDB(pdbin)
        pdbfixer = _pdb4amber.AmberPDBFixer(parm)
        self.renum = [(residue.name, residue.chain, residue.number, residue.name, residue.idx + 1) for residue in pdbfixer.parm.residues]
        # running reduce program:=============================================
        # pdbfixer.add_hydrogen writes a logging file, so we need to reimplement it
        # reduce return -1 default, dont know why
        reduce_input = io.StringIO()
        pdbfixer.write_pdb(reduce_input)
        logger.info(f"Running reduce too adjust and add H")
        logger.info(f"Command: reduce -BUILD -NUC -NOFLIP -")
        reduce_runner = SubprocessRunner("reduce -BUILD -NUC -NOFLIP -", reduce_input.getvalue(), check=False)
        reduce_runner.run()
        reduce_log = reduce_runner.stderr
        logger.debug(f"Reduce log: \n{reduce_log}\n")
        pdbfixer.parm = parmed.read_PDB(io.StringIO(reduce_runner.stdout))
        sumdict = pdbfixer._summary()
        # find histidines that might have to be changed:=====================
        logger.info(f"Assigning histidines")
        pdbfixer.assign_histidine()
        # find possible S-S in the final protein:=============================
        logger.info(f"Finding disulfide bridges")
        sslist, cys_cys_atomidx_set = pdbfixer.find_disulfide()
        pdbfixer.rename_cys_to_cyx(sslist)
        self.sslink = [(idx0 + 1, idx1 + 1) for (idx0, idx1) in sslist]
        # find possible gaps:==================================================
        self.gaplist = pdbfixer.find_gaps()
        if self.gaplist: 
            logger.error(f"Gap(s) found: {self.gaplist}, there should not be any gap in the protein")
        # count heavy atoms:==================================================
        self.missing_atom_residues = [(residue.name, residue.idx + 1, n_missing) for (residue, n_missing) in pdbfixer.find_missing_heavy_atoms()]
        if self.missing_atom_residues:
            logger.error(f"Missing heavy atom(s): {self.missing_atom_residues}, there should not be any missing heavy atom in the protein")
        # =====================================================================
        # make final output to new PDB file
        # =====================================================================
        final_coordinates = pdbfixer.parm.get_coordinates()[0]  # only the first model
        write_kwargs = dict(coordinates=final_coordinates)

        write_kwargs['increase_tercount'] = False # so CONECT record can work properly
        if sumdict['has_altlocs']:
            logger.error('There should not be any alternate location in the protein, using the first one')
            write_kwargs = dict(altlocs='first')
        # remove altlocs label
        for atom in pdbfixer.parm.atoms:
            atom.altloc = ''
            for oatom in atom.other_locations.values():
                oatom.altloc = ''
        pdbout = pdbfixer._write_pdb_to_stringio(cys_cys_atomidx_set, disulfide_conect=True, noter=False, **write_kwargs)
        self.pdbout = pdbout.getvalue()

    @property
    def success(self) -> bool:
        return self.pdbout is not None
