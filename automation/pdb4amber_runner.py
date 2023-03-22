from subprocess_runner import SubprocessRunner
from typing import Union, List
from pathlib import Path
import io


class PDB4AmberRunner(SubprocessRunner):
    def __init__(self, input: Union[str, Path, io.IOBase], name: str = '', **kwargs):
        """PDB4AmberRunner is a wrapper for pdb4amber"""
        super().__init__('pdb4amber', input, True, **kwargs)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._cmd}, {self._stdin}, {self._check}, {self._kwargs})'
