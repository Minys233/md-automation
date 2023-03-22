import logging
import io
from typing import Union
from pathlib import Path

logging.basicConfig()
logger = logging.getLogger(__name__)


class StringStream:
    def __init__(self, input: Union[str, Path, io.IOBase], name: str = ''):
        """StringStream is a wrapper for a string from multiple sources that can be used as a file handle.

        :param input: text input, from a single string, text file at some path, or an opened IOBase object
        :type input: Union[str, Path, io.IOBase]
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
        elif isinstance(input, io.IOBase):
            self._handle = io.StringIO(input.read())
        else:
            raise TypeError(f'input must be str, Path, or IOBase, not {type(input)}')
        self.name = name

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name})'
    
    def __iter__(self):
        return iter(self._handle)

    def __next__(self):
        return self._handle.__next__()

    def read(self, size=-1):
        return self._handle.read(size=size)
    
    def readline(self, size=-1):
        return self._handle.readline(size=size)
    
    def readlines(self, hint=-1):
        return self._handle.readlines(hint)
    
    def write(self, s):
        return self._handle.write(s)
    
    def seek(self, offset):
        return self._handle.seek(offset)
    
    def tell(self):
        return self._handle.tell()

