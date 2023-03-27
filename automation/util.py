import logging
import io
from typing import Union
from pathlib import Path
import parmed
from bs4 import BeautifulSoup
import urllib.request

logging.basicConfig()
logger = logging.getLogger(__name__)

# maybe shutdown soon! but we kept a local copy haha
LIGAND_EXPO_BASE = 'http://ligand-expo.rcsb.org/files/{}/{}/isdf/'
MODELS_RCSB_BASE = 'https://models.rcsb.org/v1/{}/ligand?auth_asym_id={}&label_comp_id={}&model_nums=1&encoding={}'

def ligand_from_ligand_expo(pdbid: str, resname: str, chain: str) -> str:
    url = LIGAND_EXPO_BASE.format(resname[0].upper(), resname.upper())
    logger.info(f"Fetching ligand data from {url}")
    html_doc = urllib.request.urlopen(url).read().decode('utf-8')
    soup = BeautifulSoup(html_doc, 'html.parser')
    ligand = None
    for link in soup.body.find_all('a'):
        href = link.get('href')
        
        # PDB ID _ Component ID _ Model No. _ Chain ID _ Residue No. _ mmCIF Sequence No. _ mmCIF Asym ID _ Disorder Flag _
        elem = href.split('_')
        if len(elem) < 2:
            logger.info(f"Skipping {href} with elem {elem}")
            continue
        _pdbid, _resname, _model, _chain,  = elem[0], elem[1], elem[2], elem[3]
        if pdbid.lower() == _pdbid.lower() and resname.lower() == _resname.lower() and chain.lower() == _chain.lower():
            ligand = url + href
    if ligand is None:
        logger.error(f"Cannot find ligand ({pdbid}-{resname}-{chain}) at {url}")
        return ""
    return urllib.request.urlopen(ligand).read().decode('utf-8')

def ligand_from_rcsb_model(pdbid: str, resname: str, chain: str, format: str='mol2'):
    url = MODELS_RCSB_BASE.format(pdbid, chain, resname, format)
    logger.info(f"Fetching ligand data from {url}")
    ligand_str = urllib.request.urlopen(url).read().decode('utf-8')
    # check success
    success = False
    for l in ligand_str.split('\n'):
        if not l.startswith('#'):
            success = True
            break
    if success:
        return ligand_str
    else:
        logger.error(f"Cannot find ligand ({pdbid}-{resname}-{chain}), try to check the URL: {url}")
        return ""

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


def find_ligand(pdb: Union[str, io.IOBase, 'StringStream'], cutoff: int):
    # structure.copy(parmed.Structure)
    structure = parmed.formats.PDBFile.parse(StringStream(pdb))
    water_mask = ':' + ','.join(parmed.residue.WATER_NAMES)
    structure.strip(water_mask)
    group_list = group_structure_by_bond(structure)
    print(f"Found following groups with lengths: {[len(i) for i in group_list]}")
    print(f"Groups with length <= {cutoff} are considered as ligands.")
    ligands = []
    for g in group_list:
        sele = structure[':' + ','.join([str(i+1) for i in g])]
        # resname_set = set([r.name for r in sele.residues])
        if len(g) <= cutoff: # and not resname_set.intersection(parmed.residue.WATER_NAMES):
            print(f"Found a ligand with N_res = {len(g)}: {list(sele.residues)}")
            ligands.append(sele)
    return ligands

class StringStream:
    def __init__(self, input: Union[str, Path, io.IOBase, 'StringStream'], name: str = ''):
        """StringStream is a wrapper for a string from multiple sources that can be used as a file handle.

        :param input: text input, from a single string, text file at some path, or an opened IOBase object
        :type input: Union[str, Path, io.IOBase, StringStream]
        :param name: the name of the stream, defaults to ''
        :type name: str, optional
        :raises TypeError: if input is not a string, Path, or IOBase
        """
        if isinstance(input, str) and '\n' in input:
            self._handle = io.StringIO(input)
        elif isinstance(input, str) and '\n' not in input:
            with open(input, 'r') as f:
                self._handle = io.StringIO(f.read())
        elif isinstance(input, Path):
            self._handle = io.StringIO(input.read_text())
        elif isinstance(input, (io.IOBase, StringStream)):
            self._handle = io.StringIO(input.read())
        else:
            raise TypeError(f'input must be str, Path, IOBase or StringStream, not {type(input)}')
        self.name = name

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name})'
    
    def __iter__(self):
        return iter(self._handle)

    def __next__(self):
        return self._handle.__next__()

    def read(self):
        return self._handle.read()
    
    def readline(self):
        return self._handle.readline()
    
    def readlines(self):
        return self._handle.readlines()
    
    def write(self, s):
        return self._handle.write(s)
    
    def seek(self, offset):
        return self._handle.seek(offset)
    
    def tell(self):
        return self._handle.tell()

