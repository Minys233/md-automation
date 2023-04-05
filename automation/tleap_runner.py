from subprocess_runner import SubprocessRunner
from typing import Union
from pathlib import Path
import parmed
import logging

DEFAULT_FF_TEMPLATE = \
"""
### Default forcefields ###
# Protein force field
source leaprc.protein.ff14SB
# DNA forcefield
source leaprc.DNA.OL15
# RNA forcefield
# source leaprc.RNA.OL3
# Carbonhydrate forcefield
source leaprc.GLYCAM_06j-1
# Lipid forcefield
source leaprc.lipid21
# Genral Amber forcefield
source leaprc.gaff2

"""

DEFAULT_IMPLICIT_SOLVENT_FF_TEMPLATE = \
"""
### Default implicit solvent forcefields ###
# Protein force field
source leaprc.protein.ff14SBonlysc
"""


DEFAULT_WATER_MODEL = \
"""
# using which watermodel 
source leaprc.water.tip3p
# loadAmberParams frcmod.tip3pfb
"""


"""
# load file
prot = loadpdb 10GS_for_tleap.pdb
check prot
saveamberparm prot 10GS_run.prmtop 10GS_run.inpcrd
quit
"""



class UnitNameGenerator:
    def __init__(self, prefix: str="") -> None:
        self.name_set = set()
        self.prefix = prefix
        self.current_idx = 1
    def get_name(self):
        name = f"{self.prefix}{self.current_idx}"
        self.current_idx += 1
        return name


def find_disulfide(parm: parmed.Structure):
    residues = [res for res in parm.residues if res.name in ['CYS', 'CYX']]
    cys_cys_resid_set = set()
    cys_cys_atomidx_set = set()
    for residue in residues:
        for atom in residue.atoms:
            if 'SG' in atom.name:
                for partner in atom.bond_partners:
                    if (partner.residue.name.startswith('CY') and partner.name.startswith('SG')):
                        # use tuple for hashing
                        cys_cys_resid_set.add(tuple(sorted((atom.residue.idx, partner.residue.idx))))
                        cys_cys_atomidx_set.add(tuple(sorted((atom.idx, partner.idx))))
    return sorted(cys_cys_resid_set), sorted(cys_cys_atomidx_set)

def preamble(forcefields: list=[], cmds: list=[]):
    if not forcefields:
        return DEFAULT_FF_TEMPLATE + '\n'.join(cmds)
    else:
        return '\n'.join([f"source leaprc.{ff}" for ff in forcefields] + cmds) + '\n'


ProteinNameGenerator = UnitNameGenerator("ProteinUnit")
SmallMoleculeGenerator = UnitNameGenerator("SmallMoleculeUnit")


def load_file(fpath: Union[str, Path], varname: str='', unit_registry: set=set()):
    fpath = Path(fpath).absolute()
    if fpath.suffix == '.pdb':
        cmd = f"loadpdb {str(fpath)}\n"
    elif fpath.suffix == '.prepi':
        cmd =  f"loadamberprep {fpath.name}\n"
    elif fpath.suffix == '.frcmod':
        cmd =  f"loadamberparams {fpath.name}\n"
    elif fpath.suffix == '.mol2':
        cmd =  f"loadmol2 {fpath.name}\n"
    elif fpath.suffix == '.mol3':
        cmd =  f"loadmol3 {fpath.name}\n"
    elif fpath.suffix == '.off' or fpath.suffix == '.lib':
        cmd =  f"loadoff {fpath.name}\n"
    else:
        raise ValueError(f"Unknown file type: {fpath.suffix}")
    if varname:
        cmd = f"{varname} = " + cmd
        unit_registry.add(varname)
    return cmd


def tleap_small_molecule(libfile: Union[str, Path], frcmodfile: Union[str, Path],
                         mol2file: Union[str, Path], cmds_before: list=[], cmds_after: list=[], unit_registry: set=set()):
    tleap_input = '\n'.join(cmds_before)
    tleap_input += load_file(libfile)
    tleap_input += load_file(frcmodfile)
    small_mol_name = SmallMoleculeGenerator.get_name()
    # small_mol_name = mol2file.name
    tleap_input += load_file(mol2file, small_mol_name, unit_registry=unit_registry)
    tleap_input += '\n'.join(cmds_after)
    return tleap_input


def set_box(unit: str, box: list):
    if box:
        return f"setbox {unit} {box[0]} {box[1]} {box[2]}\n"
    else:
        return ""


def tleap_standard_protein(pdbfile: str, box: list=[], extra_bonds: list=[],
                           cmds_before: list=[], cmds_after: list=[], unit_registry: set=set()):
    parm = parmed.read_PDB(pdbfile)
    tleap_input = '\n'.join(cmds_before)
    prot_name = ProteinNameGenerator.get_name()
    # prot_name = pdbfile.name
    tleap_input += load_file(pdbfile, prot_name, unit_registry=unit_registry)
    # Link disulfide bonds
    sslist, atomidx = find_disulfide(parm)
    if sslist:
        for resid1, resid2 in sslist:
            tleap_input += f"bond {prot_name}.{resid1+1}.SG {prot_name}.{resid2+1}.SG\n"
    if extra_bonds:
        for resid1, atomname1, resid2, atomname2 in extra_bonds:
            tleap_input += f"bond {prot_name}.{resid1+1}.{atomname1} {prot_name}.{resid2+1}.{atomname2}\n"
    tleap_input += set_box(prot_name, box)
    # Others? Gap? there should not have any other gaps
    tleap_input += '\n'.join(cmds_after)
    tleap_input += f'check {prot_name}\n'
    return tleap_input

def save_all(nameprefix: str, unit_registry: set):
    save_cmds = []
    for unit in unit_registry:
        save_cmds.append(f"saveAmberParm {unit} {nameprefix}.prmtop {nameprefix}.inpcrd")
    return '\n'.join(save_cmds)
    

class tLEaPRunner(SubprocessRunner):
    def __init__(self, working_dir: Union[str, Path], input_file: Union[str, Path], 
                 output_prefix: str="", forcefields: list=["protein.ff14SB", "DNA.OL15", "GLYCAM_06j-1", "lipid21", "gaff2"]):
        self.input_file = Path(input_file).absolute()
        self.tleap_input = ""
        self.working_dir = Path(working_dir).absolute()
        self.output_prefix = output_prefix
        self.forcefields = forcefields
        self.unit_registry = set()
        super().__init__('tleap -f', check=False)
    
    def run(self):
        self.tleap_input += f'logfile {str(self.working_dir / self.input_file.name)}.tleap.log\n'
        self.tleap_input += preamble(self.forcefields)
        self.tleap_input += tleap_standard_protein(str(self.input_file), unit_registry=self.unit_registry)
        self.tleap_input += save_all(str(self.working_dir / self.output_prefix), self.unit_registry)
        with open(self.working_dir / f"{self.output_prefix}.tleap.in", 'w') as f:
            f.write(self.tleap_input)
        self.append_cmd(f"{self.working_dir / f'{self.output_prefix}.tleap.in'}")
        super().run()
    
        

if __name__ == '__main__':
    print(preamble())
    print(tleap_standard_protein('automation/test_data/10GS_for_tleap.pdb', box=[]))
    print("saveAmberParm ProteinUnit1 automation/test_data/10GS.prmtop automation/test_data/10GS.inpcrd")


# leap_template = """
# {force_fields}
# {more_force_fields}
# x = loadpdb {input_pdb}
# {box_info}
# {more_leap_cmds}
# set default PBradii mbondi3
# set default nocenter on
# saveAmberParm x {prmtop} {rst7}
# quit
# """