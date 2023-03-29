# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from collections import OrderedDict
import logging
from pathlib import Path
import numpy as np
import pandas as pd

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#     === ===============================
#     ASH Neutral ASP
#     CYX SS-bonded CYS
#     CYM Negative CYS
#     GLH Neutral GLU
#     HIP Positive HIS
#     HID Neutral HIS, proton HD1 present
#     HIE Neutral HIS, proton HE2 present
#     LYN Neutral LYS
#     TYM Negative TYR
#     AR0 Neutral ARG
#     === ===============================

#     ========= ======= =========
#     Charge +1 Neutral Charge -1
#     ========= ======= =========
#     -         ASH     ASP
#     -         CYS     CYM
#     -         GLH     GLU
#     HIP       HID/HIE -
#     LYS       LYN     -
#     -         TYR     TYM
#     ARG       AR0     -
#     ========= ======= =========
    
import logging
from io import StringIO, IOBase
from pdb2pqr import io, pdb
from pdb2pqr.main import setup_molecule, drop_water, is_repairable
import pdb2pqr.forcefield as forcefield
import pdb2pqr.hydrogens as hydrogens
import pdb2pqr.debump as debump

import propka
import propka.input as pk_in
from propka.parameters import Parameters
from propka.molecular_container import MolecularContainer
import argparse
import pandas as pd
from typing import Union
from util import StringStream

from runner import Runner


class PDB2PQRRunner(Runner):
    def __init__(self, pdbin: Union[str, IOBase, StringStream], drop_water: bool=False,
                 ff: str="parse", ffout: str="amber", neutraln: bool=False, neutralc: bool=False,
                 debump: bool=True, pH: float=7.0, opt: bool=True):
        # pdbstr, argdrop_water, ff, ffout, neutraln, neutralc, argdebump, ph, opt
        self.pdbin = StringStream(pdbin).read()
        self.drop_water = drop_water
        self.ff = ff
        self.ffout = ffout
        self.neutraln = neutraln
        self.neutralc = neutralc
        self.debump = debump
        self.pH = pH
        self.opt = opt
        # stores the result
        self.pka_table = None
        self.biomolecule = None
        self.modified = None
        self.pdbout = None
        super().__init__('pdb2pqr', False)
    
    @property
    def success(self):
        pass

    def run(self,ignore_opreated: bool=False):
        if self._operated:
            if ignore_opreated:
                logging.warning(f'{self} has already been operated, running again.')
            else:
                logging.error(f'{self} has already been operated, skipping.')
                return
        self._operated = True

        definition = io.get_definitions()
        pdblist, _ = pdb.read_pdb(StringIO(self.pdbin))
        if self.drop_water:
            pdblist = drop_water(pdblist)

        self.biomolecule, definition, _ = setup_molecule(pdblist, definition, None) # no ligand assumed
        self.biomolecule.set_termini(self.neutraln, self.neutralc)
        self.biomolecule.update_bonds()

        forcefield_ = forcefield.Forcefield(self.ff, definition, None, None)# args.userff, args.usernames
        hydrogen_handler = hydrogens.create_handler()
        debumper = debump.Debump(self.biomolecule)
        if is_repairable(self.biomolecule, False):# args.ligand is not None):
            self.biomolecule.repair_heavy()
        self.biomolecule.update_ss_bridges()
        if self.debump:
            debumper.debump_biomolecule()
        self.biomolecule.remove_hydrogens()
        pka_df = self._run_propka(self.biomolecule, keep_protons=False, chains=None)
        self.biomolecule.apply_pka_values(
            forcefield_.name, self.pH,
            {f"{row['res_name']} {row['res_num']} {row['chain_id']}": row["pKa"]
            for row in pka_df if row["group_label"].startswith(row["res_name"])},)
        self.biomolecule.add_hydrogens()
        if self.debump:
            debumper.debump_biomolecule()
        hydrogen_routines = hydrogens.HydrogenRoutines(debumper, hydrogen_handler)
        if self.opt:
            hydrogen_routines.set_optimizeable_hydrogens()
            hydrogen_routines.initialize_full_optimization()
        else:
            hydrogen_routines.initialize_wat_optimization()
        hydrogen_routines.optimize_hydrogens()
        hydrogen_routines.cleanup()

        self.biomolecule.set_states()
        matched_atoms, missing_atoms = self.biomolecule.apply_force_field(forcefield_)
        name_scheme = forcefield.Forcefield(self.ffout, definition, None)
        self.biomolecule.apply_name_scheme(name_scheme)
        lines = io.print_biomolecule_atoms(atomlist=matched_atoms, pdbfile=True)
        lines[-1] += '\n'
        missing_lines = io.print_biomolecule_atoms(atomlist=missing_atoms, pdbfile=True)
        missing_lines[-1] += '\n'
        self._set_pka_table(pka_df)
        if len(missing_lines) > 1:
            self.pdbout = ''.join(lines + missing_lines)
        else:
            self.pdbout = ''.join(lines)
        
        self.modified_residues()


    def _run_propka(self, biomolecule, keep_protons: bool=False, chains: list=None):
        """Run a PROPKA calculation.
        :param biomolecule:  A biomolecule object
        :type  biomolecule:  :class:`.MolecularContainer`
        :param keep_protons:  Keep protons in the output PDB file
        :type  keep_protons:  :class:`bool`
        :param chains:  List of chains to titrate
        :type  chains:  :class:`list`
        :returns:  A list of PROPKA results
        """
        parameter_file = str((Path(propka.__file__).parent / 'propka.cfg').absolute())
        args = argparse.Namespace(display_coupled_residues=True,
                                  keep_protons=keep_protons,
                                  chains=chains,
                                  titrate_only=None,
                                  parameters=parameter_file,
                                  protonate_all=False,
                                  grid=(0.0, 14.0, 0.1))
        lines = io.print_biomolecule_atoms(atomlist=biomolecule.atoms, pdbfile=True)
        with StringIO() as fpdb:
            fpdb.writelines(lines)
            parameters = pk_in.read_parameter_file(args.parameters, Parameters())
            molecule = MolecularContainer(parameters, args)
            molecule = pk_in.read_molecule_file("placeholder.pdb", molecule, fpdb)

        molecule.calculate_pka()
        conformation = molecule.conformations["AVR"]
        rows = []
        for group in conformation.groups:
            row_dict = OrderedDict()
            atom = group.atom
            row_dict["res_num"] = atom.res_num
            row_dict["ins_code"] = atom.icode
            row_dict["res_name"] = atom.res_name
            row_dict["chain_id"] = atom.chain_id
            row_dict["group_label"] = group.label
            row_dict["group_type"] = getattr(group, "type", None)
            row_dict["pKa"] = group.pka_value
            row_dict["model_pKa"] = group.model_pka
            row_dict["buried"] = group.buried
            if group.coupled_titrating_group:
                row_dict["coupled_group"] = group.coupled_titrating_group.label
            else:
                row_dict["coupled_group"] = None
            rows.append(row_dict)
        return rows


    def _set_pka_table(self, pkadf: list):
        df = {"chainID": [], "resSeq": [], "iCode": [], "stdResName": [], "protonatedResName": [], "pKa": [], "buried": []}
        for res in self.biomolecule.residues:
            # find pka and buried in pkadf
            for pkadict in pkadf:
                if pkadict["chain_id"] == res.chain_id and pkadict["res_num"] == res.res_seq and pkadict["ins_code"].strip() == res.ins_code:
                    df['chainID'].append(res.chain_id)
                    df["resSeq"].append(res.res_seq)
                    df["iCode"].append(res.ins_code)
                    df["stdResName"].append(res.name)
                    df["protonatedResName"].append(res.ffname if hasattr(res, 'ffname') else "N/A")
                    df["pKa"].append(pkadict["pKa"])
                    df["buried"].append(pkadict["buried"])
                    break
        self.pka_table = pd.DataFrame(df)
    

    def modified_residues(self, tol: float=1):
        if self.modified is not None:
            return self.modified
        modified = self.pka_table[
            (self.pka_table['stdResName'] != self.pka_table['protonatedResName'].str[-3:]) &\
            (self.pka_table['protonatedResName'] != 'N/A')
        ]
        logger.warning(f"These residues have been protonated:")
        logger.warning(f"  chainID\tresSeq\tiCode\tstdResName\tprotonatedResName\tpKa\tburied")
        for _, dr in modified.iterrows():
            logger.warning(f"  {dr['chainID']}\t{dr['resSeq']}\t{dr['iCode']}\t{dr['stdResName']}\t{dr['protonatedResName']}\t{dr['pKa']:.1f}\t{dr['buried']}")
        dubious = modified[abs(modified.pKa - self.pH) < 1]
        if len(dubious):
            logger.warning(f"Dubious protonation state: the pKa of {len(dubious)} residues is within {tol:.1f} units of pH {self.pH:.1f}:")
            logger.warning(f"  chainID\tresSeq\tiCode\tstdResName\tprotonatedResName\tpKa\tburied")
            for _, dr in dubious.iterrows():
                logger.warning(f"  {dr['chainID']}\t{dr['resSeq']}\t{dr['iCode']}\t{dr['stdResName']}\t{dr['protonatedResName']}\t{dr['pKa']:.1f}\t{dr['buried']}")
        self.modified = modified
        return modified
