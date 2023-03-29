import logging
import io
from typing import Union
from pathlib import Path
import parmed
from bs4 import BeautifulSoup
import urllib.request

logging.basicConfig()
logger = logging.getLogger(__name__)

# https://www.wwpdb.org/documentation/file-format-content/format33
PDB_TITLE_PREFIX = (
    # sect2.html Title section
    "HEADER", "OBSLTE", "TITLE ", "SPLIT ", "CAVEAT", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "REVDAT", "SPRSDE", "JRNL  ", "REMARK"
)
PDB_MAINBODY_PREFIX = (
    # sect3.html Primary structure section
    "SEQRES", "DBREF ", # "DBREF1", "DBREF2", "SEQADV", "MODRES", #  no longer used
    # sect4.html Heterogen section
    "HET   ", "HETNAM", "HETSYN", "FORMUL", # no longer used
    # sect5.html Secondary structure section
    # "HELIX", "SHEET", "TURN", # no longer used
    # sect8.html Crystallographic and Coordinate Transformation Section
    "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3", "MTRIX1", "MTRIX2", "MTRIX3",
    # sect9.html Coordinate section
    "MODEL ", "ATOM  ", "TER   ", "HETATM", "ENDMDL", # "ANISOU", # no longer used
    # sect10.html Connectivity section
    "CONECT",
    # sect11.html Bookkeeping section
    "END   ", # "MASTER", # no longer used
)
PDB_OTHER_PREFIX = (
    # sect6.html The connectivity annotation section
    "SSBOND", "LINK  ", "CISPEP",
    # sect7.html Miscellaneous Features Section
    "SITE  ",
)



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
        self.name = input.name if isinstance(input, StringStream) and not name else name

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

