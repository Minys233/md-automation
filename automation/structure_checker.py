from collections import Counter
import parmed
import io
import logging
from runner import Runner
from typing import Union, List
from util import StringStream
from pdbfixer_runner import PDBFixerRunner
from pdb4amber_runner import PDB4AmberRunner
from pdb4amber import residue
import numpy as np
logger = logging.getLogger(__name__)
logging.basicConfig()

WATER_RESIDUE = {'WAT', 'HOH'}
# Protein residues + RNA/DNA residues + metal ions
PROTEIN_RESIDUE = set(residue.RESPROT) | set(residue.RESNA) | (set(residue.RESSOLV) - WATER_RESIDUE)
IONS_RESIDUE = set(residue.RESSOLV) - WATER_RESIDUE

def selection_with_number(parm: parmed.Structure, sele_arr: np.ndarray):
    sele = parm[sele_arr]
    nums = [a.number for idx, a in enumerate(parm.atoms) if sele_arr[idx]]
    for idx, a in enumerate(sele.atoms):
        a.number = nums[idx]
    return sele

def combine_structrure_with_number(*args):
    atom_nums = []
    for s in args:
        atom_nums += [a.number for a in s.atoms]
    structure = parmed.Structure()
    for s in args:
        structure += s
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
    # structure.copy(parmed.Structure)
    if isinstance(pdb, parmed.Structure):
        structure = pdb
    else:
        structure = parmed.formats.PDBFile.parse(StringStream(pdb))
    group_list = group_structure_by_bond(structure)
    logger.warning(f"Found following groups with (length: frequency): {Counter([len(i) for i in group_list])}")
    logger.warning(f"Groups with length <= {cutoff} and not water or ions are considered as ligands.")
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
    def __init__(self, pdbin: Union[str, io.IOBase, StringStream], pdbid: str="", mutations: str="", replace_nonstandard: bool=True, ligand_cutoff: int=5):
        self._pdbin_handle = StringStream(pdbin)
        self.pdbid = pdbid
        self.header_lines = []
        self.mainbody_lines = []
        self.other_lines = []
        self._parse()

        logger.warning("===================================")
        logger.warning("===== Checking Input PDB File =====")
        logger.warning("===================================")
        self.check_structure(self._pdbin_handle.read())
        self._pdbin_handle.seek(0)
        # logger.warning("===================================")
        # logger.warning("=========== Checking End ==========")
        # logger.warning("===================================")
        self.mutations = mutations
        self.replace_nonstandard = replace_nonstandard
        self.ligand_cutoff = ligand_cutoff
        logger.warning("===================================")
        logger.warning("============= Fixing ==============")
        logger.warning("===================================")
        self.parm = self._fix_pdb()
        logger.warning("===================================")
        logger.warning("=========== Seperating ============")
        logger.warning("===================================")
        self.protein, self.heterogen, self.water = None, None, None
        self.protein_groups, self.heterogen_groups, self.water_groups = None, None, None
        self._seperate()
        logger.warning("===================================")
        logger.warning("=========== Seperating ============")
        logger.warning("===================================")



        logger.warning("===================================")
        logger.warning("===== Checking Output PDB File ====")
        logger.warning("===================================")
        pdb2check = self.parmed2pdb(self.parm)
        self.check_structure(pdb2check)
    
    @property
    def protein_pdbstr(self):
        return self.parmed2pdb(self.protein)

    @property
    def water_pdbstr(self):
        return self.parmed2pdb(self.water)
    
    @property
    def heterogen_pdbstr(self):
        return self.parmed2pdb(self.heterogen)

    def parmed2pdb(self, parm: parmed.Structure) -> str:
        pdbstream = io.StringIO()
        parm.save(pdbstream, format='pdb', renumber=False, altlocs='first',
                  write_anisou=False, charmm=False, use_hetatoms=True,
                  standard_resnames=False, increase_tercount=True,
                  write_links=True, conect=True)
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

        prot_mask_arr = np.zeros(len(self.parm.atoms), dtype=bool)
        het_mask_arr = np.zeros(len(self.parm.atoms), dtype=bool)
        water_mask_arr = np.zeros(len(self.parm.atoms), dtype=bool)

        for idx, atom in enumerate(self.parm.atoms):
            if atom.residue.name in PROTEIN_RESIDUE:
                prot_mask_arr[idx] = True
            elif atom.residue.name in WATER_RESIDUE:
                water_mask_arr[idx] = True
            else:
                het_mask_arr[idx] = True
        # selected atoms' number is -1, need to fix it here to align with the original PDB file
        # so we use the selection array ourselves and set the atom.number manually.
        protein = selection_with_number(self.parm, prot_mask_arr)
        heterogen = selection_with_number(self.parm, het_mask_arr)
        water = selection_with_number(self.parm, water_mask_arr)

        logger.warning(f"Find ligands (residue name mask) with residue names: ")
        logger.warning(f"These ligands are saved to self.heterogen along with self.protein and self.water.")
        logger.warning(f"  index\tresName\tresSeq\tchainID")
        for idx, res in enumerate(heterogen.residues):
            logger.warning(f"  {idx}\t{res.name}\t{res.number}\t{res.chain}")
        self.protein = protein
        self.heterogen = heterogen
        self.water = water
        
        logger.warning(f"Find ligands (topology Union-Find) with residue names: ")
        logger.warning(f"These ligands are saved to self.protein_groups along with self.heterogen_groups and self.water_groups.")
        protein_groups, ligand_groups, water_groups = find_ligand(self.parm, cutoff=self.ligand_cutoff)
        for idx, ligand in enumerate(ligand_groups):
            logger.warning(f"  Group {idx}:")
            logger.warning(f"  index\tresName\tresSeq\tchainID")
            for iidx, res in enumerate(ligand.residues):
                logger.warning(f"  {iidx}\t{res.name}\t{res.number}\t{res.chain}")
        self.protein_groups = protein_groups
        self.heterogen_groups = ligand_groups
        self.water_groups = water_groups

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
        # TODO: Check SSLINK records

        
if __name__ == "__main__":
    sep = StructureChecker('/home/yaosen/md-automation/automation/test_data/1UCY.pdb')