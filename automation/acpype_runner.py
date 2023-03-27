import acpype
from acpype.topol import ACTopol
import io
from util import StringStream, ligand_from_rcsb_model, ligand_from_ligand_expo, find_ligand
from typing import Union
from runner import Runner
from subprocess_runner import SubprocessRunner
import parmed
from pathlib import Path
import openbabel

import logging

logging.basicConfig()
logger = logging.getLogger(__name__)


def add_H4pH_ligand_obabel(inp: Union[str, Path, io.IOBase, StringStream], pH: float, format: str='mol', name: str = ''):
    """Adjust the H given input file, using openbabel, better use it with sdf or mol2 format as input.
    """
    inp = StringStream(inp, name)
    obabel_runner = SubprocessRunner(f'obabel -i{format} -omol -p {pH} --partialcharge eem', stdin=inp.read(), check=True)
    obabel_runner.run()
    return obabel_runner.stdout


class AcpypeRunner(Runner):
    def __init__(self, mol_or_sdf_inp: Union[str, Path, io.IOBase, StringStream], name: str = '',
                 pH: float=7.0, charge_type: str='bcc', atom_type: str='gaff2', result_dir: Union[str, Path]=None, out_topol: str='', **kwargs):
        self.molin = add_H4pH_ligand_obabel(StringStream(mol_or_sdf_inp, name), pH, name=name)
        self._molecule = None
        self.name = name
        self.pH = pH
        # ["gas", "bcc", "user"]
        assert charge_type in ["gas", "bcc", "user"], f'charge_type must be one of ["gas", "bcc", "user"], but got {charge_type}'
        self.charge_type = charge_type
        # ["gaff", "amber", "gaff2", "amber2"]
        assert atom_type in ["gaff", "amber", "gaff2", "amber2"], f'atom_type must be one of ["gaff", "amber", "gaff2", "amber2"], but got {atom_type}'
        self.atom_type = atom_type
        # ["gmx", "cns", "charmm", 'all']
        assert out_topol in ["gmx", "cns", "charmm", 'all', ''], f'out_topol must be one of ["gmx", "cns", "charmm", "all", ""], but got {out_topol}'
        self.out_topol = out_topol
        self.kwargs = kwargs
        super().__init__(f'See acpype for more detail', False)
        # result stored here
        self.result_dir = result_dir
        # set corresponding working dir for acpype
        self.create_result_dir()


    def __repr__(self) -> str:
        return f'<AcpypeRunner: {self.name} pH={self.pH} charge_type={self.charge_type} atom_type={self.atom_type} temp_dir={self.result_dir} out_topol={self.out_topol} kwargs={self.kwargs}>'

    def create_result_dir(self):
        if not self.result_dir:
            self.result_dir = Path(f'acpype_{self.name}')
        else:
            self.result_dir = Path(self.result_dir)
        logger.info(f"Using temporary directory: {self.result_dir}")
        if self.result_dir.is_dir():
            logger.warning(f'{self.result_dir} already exists, will overwriting files.')
        else:
            self.result_dir.mkdir(parents=True)
        self.kwargs['basename'] = str(self.result_dir.absolute() / self.result_dir.name)
        return self.result_dir


    def run(self, ignore_opreated: bool=False):
        """This function is modified from acpype.run() funtion.
        Most of the parameters are fixed as we has an already fixed pdb file.
        """
        if self._operated:
            if ignore_opreated:
                logging.warning(f'{self} has already been operated, running again.')
            else:
                logging.error(f'{self} has already been operated, skipping.')
                return
        self._operated = True
        logger.info(f"Using Openbabel (version {openbabel.__version__}) to correct the molecule to pH={self.pH}")
        
        molpath = self.result_dir / (self.result_dir.name + '.mol')
        with open(molpath, 'w') as f:
            f.write(self.molin)
        logger.info(f"Using acpype package with version: {acpype.__version__}")
        self._molecule = ACTopol(str(molpath.absolute()), chargeType=self.charge_type, atomType=self.atom_type,
                           outTopol=self.out_topol, **self.kwargs)
        self._molecule.createACTopol()
        self._molecule.createMolTopol()

    
    @property
    def success(self) -> bool:
        # which files are used in tleap?
        pkl_exist = (self.result_dir / (self.result_dir.name + '.pkl')).is_file()
        frcmod_exist = (self.result_dir / (self.result_dir.name + '_AC.frcmod')).is_file()
        inpcrd_exist = (self.result_dir / (self.result_dir.name + '_AC.inpcrd')).is_file()
        prmtop_exist = (self.result_dir / (self.result_dir.name + '_AC.prmtop')).is_file()
        return pkl_exist and frcmod_exist and inpcrd_exist and prmtop_exist
        

