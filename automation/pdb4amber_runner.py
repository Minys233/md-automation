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

def add_H4pH_protein_pdb2pqr(inp: Union[str, Path, io.IOBase, StringStream], pH: float, name: str = ''):
    """Adjust the H given input file, using pdb2pqr.
    Lose all information other than ATOM and HETATM.
    """
    inp = StringStream(inp, name)
    obabel_runner = SubprocessRunner(f'pdb2pqr --ff PARSE --keep-chain --pdb-output 10GS_fix.pqr.pdb --include-header --with-ph {pH} --pH {pH} --titration-state-method propka {pdbin} {pqrout}', stdin=inp.read(), check=True)
    obabel_runner.run()
    out = obabel_runner.stdout
    return out.getvalue()

# pdb2pqr --ff PARSE --keep-chain --pdb-output 10GS_fix.pqr.pdb --include-header --with-ph 7 --pH 7 --titration-state-method propka 10GS_fix.pdb 10GS_fix.pqr


class PDB4AmberRunner(Runner):
    def __init__(self, pdbin: Union[str, Path, io.IOBase, StringStream], name: str = '',
                 rm_water: bool=False, prot_only: bool=False, amber_only: bool=False,
                 custom_mask: str="", **kwargs):
        """PDB4AmberRunner is a wrapper for the pdb4amber command line tool. 
        It is used to convert a PDB file to another PDB file that can be used as input for tleap.

        :param input: The input pdb file
        :type input: Union[str, Path, io.IOBase]
        :param name: The name of the pdbin, defaults to ''
        :type name: str, optional
        """
        self.pdbin = StringStream(pdbin, name)
        self.rm_water = rm_water
        self.prot_only = prot_only
        self.amber_only = amber_only
        self.custom_mask = custom_mask
        self.kwargs = kwargs
        cmd = "pdb4amber --reduce --model 0" + " --dry" if rm_water else "" + " --prot" if prot_only else "" + \
              " --amber-compatible-residues" if amber_only else "" + " --strip " + custom_mask if custom_mask else ""
        super().__init__(f'Equivalent to run: {cmd}', False)
        # result stored here
        self.renum = None
        self.prot = None
        self.amber = None
        self.custom = None
        self.nonprot = None
        self.water = None
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
        # find and save non-standard Amber residues:==========================
        ns_names = pdbfixer.find_non_standard_resnames()
        logger.info(f"Non-standard residue names: {', '.join(ns_names)}")
        ns_mask = ':' + ','.join(ns_names)
        nonprot = io.StringIO()
        if ns_mask != ':':
            pdbfixer.parm[ns_mask].save(nonprot, format='pdb')
        self.nonprot = nonprot.getvalue()
        # remove water if -d option used:=====================================
        # remove water should be done before only protein residues, since waters are not protein residues
        water_mask = ':' + ','.join(parmed.residue.WATER_NAMES)
        water = io.StringIO()
        pdbfixer.parm[water_mask].save(water, format='pdb')
        self.water = water.getvalue()
        if self.rm_water:
            logger.info(f"Removing water and saving the waters")
            water = io.StringIO()
            pdbfixer.parm[water_mask].save(water, format='pdb')
            self.water = water.getvalue()
            pdbfixer.remove_water()
        # keep only protein residues:======================
        prot_mask = ':' + ','.join(_pdb4amber.pdb4amber.RESPROT)
        prot = io.StringIO()
        pdbfixer.parm[prot_mask].save(prot, format='pdb')
        self.prot = prot.getvalue()
        if self.prot_only:
            logger.info(f"Keeping only protein residues")
            pdbfixer.parm.strip('!' + prot_mask)
        # keep only Amber compatible residues:===============================
        amber_mask = ':' + ','.join(_pdb4amber.pdb4amber.AMBER_SUPPORTED_RESNAMES)
        amber = io.StringIO()
        pdbfixer.parm[amber_mask].save(amber, format='pdb')
        self.amber = amber.getvalue()
        if self.amber_only:
            logger.info(f"Keeping only Amber compatible residues")
            pdbfixer.parm.strip('!' + amber_mask)
        # strip atoms with given mask =====================================
        if self.custom_mask:
            logger.info(f"Stripping atoms with mask: {self.custom_mask}")
            custom = io.StringIO()
            pdbfixer.parm[self.custom_mask].save(custom, format='pdb')
            self.custom = custom.getvalue()
            pdbfixer.parm.strip(self.custom_mask)
        else:
            logger.info(f"No custom mask given, leaving self.custom is None")
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

    def write(self, outfile: Union[str, Path]):
        """Write the output PDB file to a file.
        :param outfile: the output file
        :type outfile: Union[str, Path]
        """
        if self._operated:
            with open(outfile, 'w') as f:
                f.write(self.pdbout)
        else:
            logger.error(f'{self} has not been operated, skipping.')

    @property
    def success(self) -> bool:
        return self.pdbout is not None
