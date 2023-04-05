from collections import Counter
import parmed
from parmed.formats.pdb import _is_hetatm
import io
import logging
from runner import Runner
from typing import Union, List
from util import StringStream, ligand_from_rcsb_model, protein_from_rcsb
from pdbfixer_runner import PDBFixerRunner
from pdb4amber_runner import PDB4AmberRunner
from pdb2pqr_runner import PDB2PQRRunner
from acpype_runner import AcpypeRunner
from tleap_runner import tLEaPRunner
from pdb4amber import residue
import numpy as np
from pathlib import Path
logger = logging.getLogger(__name__)
logging.basicConfig()

WATER_RESIDUE = {'WAT', 'HOH'}
# Protein residues + RNA/DNA residues + metal ions
PROTEIN_RESIDUE = set(residue.RESPROT) | set(residue.RESNA) | (set(residue.RESSOLV) - WATER_RESIDUE)
IONS_RESIDUE = set(residue.RESSOLV) - WATER_RESIDUE

def selection_with_number(parm: parmed.Structure, sele_arr: np.ndarray):
    sele = parm[sele_arr]
    for idx, atom in enumerate(sele.atoms, 1):
        atom.number = idx
    return sele

def combine_structrure_with_number(*args):
    atom_nums = [i for i in range(1, sum([len(s.atoms) for s in args]) + 1)]
    structure = parmed.Structure()
    for idx, a in enumerate(structure.atoms):
        a.number = atom_nums[idx]
    return structure

def group_structure_by_bond(structure: parmed.Structure):
    """Group residues by linkage between residues, with Union-Find algorithm."""
    # initialize the Union-Find algorithm
    uf = dict()
    for idx, residue in enumerate(structure.residues):
        uf[idx] = [idx]
    # find the connected residues
    for residue in structure.residues:
        for atom in residue:
            for partner in atom.bond_partners:
                if partner.residue.idx != residue.idx:
                    # find the root of the partner
                    root = partner.residue.idx
                    while root in uf and root != uf[root][0]:
                        root = uf[root][0]
                    # find the root of the residue
                    root2 = residue.idx
                    while root2 in uf and root2 != uf[root2][0]:
                        root2 = uf[root2][0]
                    # merge the two groups
                    if root != root2:
                        uf[root2] += uf[root]
                        uf[root] = [root2]
    # return the groups
    groups = []
    for residue in structure.residues:
        if residue.idx in uf and residue.idx == uf[residue.idx][0]:
            groups.append(uf[residue.idx])
    groups = [sorted(i) for i in groups]
    return sorted(groups, key=lambda x: len(x))


def find_ligand(pdb: Union[str, io.IOBase, 'StringStream', parmed.Structure], cutoff: int):
    if isinstance(pdb, parmed.Structure):
        structure = pdb
    else:
        structure = parmed.formats.PDBFile.parse(StringStream(pdb))
    group_list = group_structure_by_bond(structure)
    logger.warning(f"Found following groups with (length: frequency): {Counter([len(i) for i in group_list])}")
    water_group, protein_group, ligand_group = [], [], []
    # filter out the water to one group
    # filter out the ions to protein group
    for g in group_list:
        if len(g) == 1:
            if structure.residues[g[0]].name in WATER_RESIDUE:
                water_group.append(g)
            elif structure.residues[g[0]].name in IONS_RESIDUE:
                protein_group.append(g)
            else:
                ligand_group.append(g)
        elif len(g) <= cutoff:
            ligand_group.append(g)
        else:
            protein_group.append(g)
    # water selection, merge all water groups
    water_group = [i for j in water_group for i in j]
    water_sele_arr = np.zeros(len(structure.atoms), dtype=bool)
    for idx, atom in enumerate(structure.atoms):
        if atom.residue.idx in water_group:
            water_sele_arr[idx] = True
    water_sele = [selection_with_number(structure, water_sele_arr)]
    # protein selection, not merge the groups
    protein_sele = []
    for g in protein_group:
        sele_arr = np.zeros(len(structure.atoms), dtype=bool)
        for idx, atom in enumerate(structure.atoms):
            if atom.residue.idx in g:
                sele_arr[idx] = True
        protein_sele.append(selection_with_number(structure, sele_arr))
    # ligand selection, not merge the groups
    ligand_sele = []
    for g in ligand_group:
        sele_arr = np.zeros(len(structure.atoms), dtype=bool)
        for idx, atom in enumerate(structure.atoms):
            if atom.residue.idx in g:
                sele_arr[idx] = True
        ligand_sele.append(selection_with_number(structure, sele_arr))
    return protein_sele, ligand_sele, water_sele

_DEFAULT_HEADER_PREFIX = (
    # sect2.html Title section
    "HEADER", "OBSLTE", "TITLE ", "SPLIT ", "CAVEAT", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "REVDAT", "SPRSDE", "JRNL  ", "REMARK",
    # sect3.html Primary structure section
    "SEQRES", "DBREF ", "DBREF1", "DBREF2", "SEQADV", "MODRES", #  no longer used
    # sect4.html Heterogen section
    "HET   ", "HETNAM", "HETSYN", "FORMUL", # no longer used
    # sect5.html Secondary structure section
    "HELIX", "SHEET", "TURN", # no longer used
    # sect8.html Crystallographic and Coordinate Transformation Section
    "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3", "MTRIX1", "MTRIX2", "MTRIX3",
    # sect6.html The connectivity annotation section
    "SSBOND", "LINK  ", "CISPEP",
    # sect7.html Miscellaneous Features Section
    "SITE  ",
)

_DEFAULT_MAINBODY_PREFIX = (
    # sect9.html Coordinate section
    "MODEL ", "ATOM  ", "TER   ", "HETATM", "ENDMDL", # "ANISOU", # no longer used
    # sect10.html Connectivity section
    "CONECT",
    # sect11.html Bookkeeping section
    "END   ", # "MASTER", # no longer used
)


class StructureChecker:
    def __init__(self, working_dir: Union[str, Path], pdbid: str, pdbin: Union[str, io.IOBase, StringStream]=None,
                 mutations: str="", replace_nonstandard: bool=True, ligand_cutoff: int=1, pH: float=7.0, charge_type: str='gas'):
        self.working_dir = Path(working_dir)
        self.working_dir.mkdir(parents=True, exist_ok=True)
        self.pdbid = pdbid.upper()
        logging.warning(f"Working directory: {self.working_dir}")
        # read the pdb file
        if pdbin:
            self._pdbin_handle = StringStream(pdbin)
        elif self.pdbid:
            self._pdbin_handle = StringStream(protein_from_rcsb(pdbid))
        # write the pdb file to working directory
        self._pdbin_handle.seek(0)
        with open(self.working_dir / f"{pdbid}_input.pdb", "w") as f:
            f.write(self._pdbin_handle.read())
        self._pdbin_handle.seek(0)
        # parse the pdb file
        self.header_lines = []
        self.mainbody_lines = []
        self.other_lines = []
        self._parse()
        logger.warning("===================================")
        logger.warning("===== Checking Input PDB File =====")
        logger.warning("===================================")
        self.check_structure(self._pdbin_handle.read())
        self._pdbin_handle.seek(0)
        # some parameters
        self.mutations = mutations
        self.replace_nonstandard = replace_nonstandard
        self.ligand_cutoff = ligand_cutoff
        self.pH = pH
        self.charge_type = charge_type
        logger.warning("===================================")
        logger.warning("============= Fixing ==============")
        logger.warning("===================================")
        self.parm = self._fix_pdb()
        logger.warning("===================================")
        logger.warning("=========== Seperating ============")
        logger.warning("===================================")
        self.seperation = None
        self.heterogen_mol = dict()
        self._seperate()
        logger.warning("===================================")
        logger.warning("============ PDB2PQR ==============")
        logger.warning("===================================")
        self.pdb2pqr()
        # logger.warning("===================================")
        # logger.warning("=========== PDB4Amber =============")
        # logger.warning("===================================")
        # self.pdb4amber()
        logger.warning("===================================")
        logger.warning("============= tLEaP ===============")
        logger.warning("===================================")
        self.tleap()
        logger.warning("===================================")
        logger.warning("============= ACPYPE ==============")
        logger.warning("===================================")
        self.acpype()
        # logger.warning("===================================")
        # logger.warning("===== Checking Output PDB File ====")
        # logger.warning("===================================")
        # pdb2check = self.parmed2pdb(self.parm)
        # self.check_structure(pdb2check)
    
    def parmed2pdb(self, parm: parmed.Structure) -> str:
        pdbstream = io.StringIO()
        parm.save(pdbstream, format='pdb', renumber=True, altlocs='first',
                  write_anisou=False, charmm=False, use_hetatoms=True,
                  standard_resnames=False, increase_tercount=False,
                  write_links=False, conect=True)
        pdbstr = pdbstream.getvalue()
        if self.header_lines:
            return ''.join(self.header_lines) + '\n'.join(self._keep_mainbody_lines(pdbstr))
        else:
            return pdbstr

    def _parse(self):
        """Parse the PDB file."""
        self._pdbin_handle.seek(0)
        for line in self._pdbin_handle:
            if line.startswith(_DEFAULT_HEADER_PREFIX):
                self.header_lines.append(line)
            elif line.startswith(_DEFAULT_MAINBODY_PREFIX):
                self.mainbody_lines.append(line)
            else:
                logger.debug(f'Skipping line with unused prefix: {line.strip()}')
                self.other_lines.append(line)
        self._pdbin_handle.seek(0)
    

    def _keep_mainbody_lines(self, lines: Union[str, io.IOBase, StringStream]) -> str:
        """Keep only the mainbody lines."""
        if isinstance(lines, str):
            lines = lines.split('\n')
        else:
            lines.seek(0)
            lines = lines.read().split('\n')
        main_lines = [l for l in lines if l.startswith(_DEFAULT_MAINBODY_PREFIX)]
        return main_lines


    def _fix_pdb(self):
        """Fix the PDB file."""
        self._pdbin_handle.seek(0)
        pdbfixer = PDBFixerRunner(self._pdbin_handle.read(), remove_het=False, keep_water=True, mutations=self.mutations,
                                  replace_nonstandard=self.replace_nonstandard, addH=False)
        pdbfixer.run()
        assert pdbfixer.success, "PDBFixer failed to fix the PDB file."
        fixed_pdbstr = ''.join(self.header_lines) + '\n'.join(self._keep_mainbody_lines(pdbfixer.pdbout))
        return parmed.read_PDB(StringStream(fixed_pdbstr))


    def _seperate(self):
        """Seperate the PDB file into protein, heterogen, and water."""
        # new seperate method
        chain_set = set([r.chain for r in self.parm.residues])
        chain_dict = {}
        logger.warning(f"A bonding group is defined as a set of atoms that form a single molecule (connected graph)")
        logger.warning(f"Bonding groups that (length <= {self.ligand_cutoff}) and (not water) and (not ions) are considered as ligands")
        for chain in chain_set:
            sele_arr = np.array([a.residue.chain == chain for a in self.parm.atoms])
            chain_sele = selection_with_number(self.parm, sele_arr)
            chain_pro, chain_lig, chain_wat = find_ligand(chain_sele, cutoff=self.ligand_cutoff)
            logger.warning(f"Parsing chain with id = {chain}")
            logger.warning(f"  Protein: {chain_pro}; Ligand: {chain_lig}; Water: {chain_wat}")
            # maybe covalent ligand
            for item in chain_pro:
                het_res_mask = np.array([_is_hetatm(r.name) for r in item.residues])
                if het_res_mask.sum():
                    het_res_idx = np.where(het_res_mask)[0]
                    covalent_lig_res = [item.residues[r] for r in het_res_idx]
                    logger.warning(f"  Found covalent ligand in chain {chain}: {covalent_lig_res}")
                    for cr in covalent_lig_res:
                        logger.warning(f"    For residue {cr.name}:")
                        for atom in cr.atoms:
                            for other in atom.bond_partners:
                                if atom.residue.name != other.residue.name:
                                    # same chain but different residue
                                    logger.warning(f"    Bond: (NAME <{atom.name}> ATMNUM <{atom.number}> @ RESNAME <{atom.residue.name}> RESNUM <{atom.residue.number}>) <-> (NAME <{other.name}> ATMNUM <{other.number}> @ RESNAME <{other.residue.name}> RESNUM <{other.residue.number}>)")
                                    # TODO: maybe we should keed the records
                else:
                    logger.warning(f"  Chain {chain} has no covalent ligand (all residues are standard residues)")
            # non-covalent ligand

            if not chain_lig:
                logger.warning(f"  Chain {chain} has no ligands")
                continue
            for ligand in chain_lig:
                logger.warning(f"  Ligand with name {[r.name for r in ligand.residues]}")
                if len(ligand.residues) == 1:
                    logger.warning(f"  Small molecule, downloading from RCSB")
                    self.heterogen_mol[ligand] = ligand_from_rcsb_model(self.pdbid, ligand.residues[0].name, ligand.residues[0].chain, 'sdf')
                else:
                    logger.warning(f"  Ligand has residues > 1, skip downloading")

            chain_dict[chain] = [chain_pro, chain_lig, chain_wat]
        self.seperation = chain_dict

    def check_structure(self, pdbin: Union[str, io.IOBase, StringStream]):
        if not isinstance(pdbin, str):
            pdbin = pdbin.read()
        structure = parmed.formats.PDBFile.parse(StringStream(pdbin))
        if structure is None:
            logger.error("Structure is not valid.")
            return
        logger.warning("Checking the structure...")
        # Check the number of MODELs
        num_model = structure.get_coordinates().shape[0]
        if num_model > 1:
            logger.warning(f"  1. Num. MODELs: {num_model} > 1. Only the first model will be used.")
        else:
            logger.warning(f"  1. Num. MODELs: {num_model} = 1. OK.")
        # Check MODRES records
        modres = []
        for line in pdbin.split('\n'):
            if line.startswith('MODRES'):
                code, res, chain, seqnum, inser, stdres, comment = line[7:11], line[12:15], line[16], int(line[18:22]), line[22], line[24:27], line[29:].strip()
                modres.append((code, res, chain, seqnum, inser, stdres, comment))
        if len(modres):
            logger.warning(f"  2. MODRES records found: ")
            logger.warning("    index\tidCode\tresName\tchainID\tseqNum\tiCode\tstdRes\tcomment")
            for i, m in enumerate(modres):
                logger.warning(f"    {i}\t{m[0]}\t{m[1]}\t{m[2]}\t{m[3]}\t{m[4]}\t{m[5]}\t{m[6]}")
        else:
            logger.warning("  2. No MODRES records found. OK.")
        # Check altlocs
        altlocs = [atom.altloc for atom in structure.atoms]
        altloc_counter = Counter(altlocs)
        if len(altloc_counter) > 1:
            logger.warning(f"  3. Altlocs found: ")
            logger.warning(f"    {altloc_counter}")
            logger.warning(f"    (Using the first altloc bydefault)")
        else:
            logger.warning("  3. No altlocs found. OK")
        
        # Check for missing residues
        # TODO: implement it in parmed context
        from pdbfixer import PDBFixer
        pdb2fix = self.parmed2pdb(structure)
        fixer = PDBFixer(pdbfile=io.StringIO(pdb2fix))
        fixer.findMissingResidues()
        if len(fixer.missingResidues):
            logger.warning(f"  4. Found missing residues:")
            logger.warning("    index\tchainID\tresSeq\tresName")
            cids = [c.id for c in fixer.topology.chains()]
            for idx, ((cid, rid), rn) in enumerate(fixer.missingResidues.items()):
                logger.warning(f"    {idx}\t{cids[cid]}\t{rid+1}\t{rn}")
        else:
            logger.warning("  4. No missing residues found. OK.")

        # Check for nonstandard residues
        fixer.findNonstandardResidues()
        if len(fixer.nonstandardResidues):
            logger.warning(f"  5. Found nonstandard protein residues: ")
            logger.warning("    index\tchainID\tresSeq\tresName\tstdRes")
            for idx, (nstd, std) in enumerate(fixer.nonstandardResidues):
                logger.warning(f"    {idx}\t{nstd.chain.id}\t{nstd.id}\t{nstd.name}\t{std}")
        else:
            logger.warning("  5. No nonstandard residues found. OK.")
        # Check for heterogens (ligands, ions, etc.)
        def _het2delete(fixer) -> list:
            het_residues, waters = [], []
            for residue in fixer.topology.residues():
                if residue.name not in PROTEIN_RESIDUE:
                    if residue.name in WATER_RESIDUE:
                        waters.append((residue.name, residue.chain.id, residue.id))
                    else:
                        het_residues.append((residue.name, residue.chain.id, residue.id))
            return het_residues, waters
    
        het_res, wat_res = _het2delete(fixer)
        if len(het_res) or len(wat_res):
            logger.warning(f"  6. Found heterogens: ")

            if len(het_res):
                logger.warning(f"  6.1 Heterogens: ")
                logger.warning(f"    resName\tchainID\tresSeq")
                for h in het_res:
                    logger.warning(f"    {h[0]}\t{h[1]}\t{h[2]}")
            else:
                logger.warning("  6.1 No heterogens found. OK.")
            if len(wat_res):
                logger.warning(f"  6.2 Waters: {len(wat_res)} water residues found.")    
            else:
                logger.warning("  6.2 No waters found. OK.")        
        else:
            logger.warning("  6. No heterogens found. OK.")

        # Check for missing atoms
        fixer.findMissingAtoms()
        if len(fixer.missingAtoms) or len(fixer.missingTerminals):
            logger.warning(f"  7. Found missing heavy atoms: ")
            if len(fixer.missingAtoms):
                logger.warning(f"  7.1 Missing atoms:")
                logger.warning(f"    index\tchainID\tresSeq\tresName")
                for idx, (res, matm) in enumerate(fixer.missingAtoms.items()):
                    logger.warning(f"    {idx}\t{res.chain.id}\t{res.id}\t{matm}")
            else:
                logger.warning("  7.1 No missing atoms found. OK.")
            if len(fixer.missingTerminals):
                logger.warning(f"  7.2 Missing terminal atoms:")
                logger.warning(f"    Missing terminal atoms:")
                logger.warning(f"    index\tchainID\tresSeq\tresName")
                for idx, (res, matm) in enumerate(fixer.missingTerminals.items()):
                    logger.warning(f"    {idx}\t{res.chain.id}\t{res.id}\t{matm}")
            else:
                logger.warning("  7.2 No missing terminal atoms found. OK.")
        else:
            logger.warning("  7. No missing heavy atoms found. OK.")
        # check for gaps
        from pdb4amber import AmberPDBFixer
        amberfixer = AmberPDBFixer(parm=structure)
        gaplist = amberfixer.find_gaps()
        if gaplist:
            logger.warning(f"  8. Found gaps in the protein backbone:")
            logger.warning("    Gap(s) found:")
            logger.warning("    index\tgap(A)\tresName1\tresSeq1\tresName2\tresSeq2")
            for idx, (d, resname0, resid0, resname1, resid1) in enumerate(gaplist):
                # convert to 1-based
                logger.warning(f"    {d:.3f}\t{resname0}\t{resid0 + 1}\t{resname1}\t{resid1 + 1}")
        else:
            logger.warning(f"  8. No gaps found in the protein backbone. OK.")
        # TODO: Check SSLINK records

    def pdb2pqr(self):
        # self.seperation
        newpro = {chain: [] for chain in self.seperation.keys()}
        for idx, (chainid, contents) in enumerate(self.seperation.items()):
            logger.warning(f"Run PDB2PQR for chain {chainid}: protein part")
            for struc in contents[0]:
                prot_group_pdbstr = self.parmed2pdb(struc)
                pdb2pqr_runner = PDB2PQRRunner(prot_group_pdbstr, drop_water=False, ff='amber', ffout='amber',
                                               neutraln=False, neutralc=False, debump=True, pH=self.pH, opt=True)
                pdb2pqr_runner.run()
                # TODO: maybe need to add header lines and CONECT records
                # But we currently assume that the input pdb here all has standard residues without covalent bonds
                newpro[chainid].append(parmed.read_PDB(io.StringIO(pdb2pqr_runner.pdbout)))

        for chainid in newpro.keys():
            self.seperation[chainid][0] = newpro[chainid]


    def acpype(self):
        for ligand, ligand_molstr in self.heterogen_mol.items():
            result_dir = self.working_dir / f"{self.pdbid}_ligand_{ligand.residues[0].chain}_{ligand.residues[0].name}"
            acpype_runner = AcpypeRunner(ligand_molstr, name='', pH=self.pH,
                                         charge_type=self.charge_type, atom_type='gaff2',
                                         result_dir=result_dir, out_topol='', level=20)
            acpype_runner.run()

    
    def tleap(self):
        for chainid, contents in self.seperation.items():
            for pro in contents[0]:
                pdbname = f"{self.pdbid}_protein_{pro.residues[0].chain}.pdb"
                pdbpath = self.working_dir / pdbname
                with open(pdbpath, 'w') as f:
                    f.write(self.parmed2pdb(pro))
                tleap_runner = tLEaPRunner(self.working_dir, pdbpath, pdbpath.name)
                tleap_runner.run()
            
            for wat in contents[2]:
                pdbname = f"{self.pdbid}_water_{wat.residues[0].chain}.pdb"
                pdbpath = self.working_dir / pdbname
                with open(pdbpath, 'w') as f:
                    f.write(self.parmed2pdb(wat))
                tleap_runner = tLEaPRunner(self.working_dir, pdbpath, pdbpath.name)
                tleap_runner.run()



if __name__ == "__main__":
    sep = StructureChecker('/home/yaosen/10GS-SC', '10GS',)
    