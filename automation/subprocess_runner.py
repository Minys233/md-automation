from subprocess import run
from typing import List, Union
import logging
from runner import Runner

logger = logging.getLogger(__name__)


class SubprocessRunner:
    _USED_KWARGS = ['capture_output', 'encoding', 'shell', 'input', 'errors', 'check']

    def __init__(self, cmd: Union[List, str], stdin: str="", check: bool=False, **kwargs):
        """Run a subprocess and capture its output when we ask for it.

        :param cmd: command to be executed
        :type cmd: Union[List, str]
        :param stdin: data to be passed to stdin, defaults to ""
        :type stdin: str, optional
        :param check: check if the subporocess return 0, defaults to False
        :type check: bool, optional
        """
        if isinstance(cmd, str):
            cmd = cmd.split()
        self._cmd = cmd
        self._stdin = stdin
        self._check = check
        self._proc = None
        self._operated = False
        self._kwargs = kwargs
    
    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._cmd}, {self._stdin}, {self._kwargs})'

    @property
    def operated(self) -> bool:
        return self._operated

    @property
    def cmd(self) -> str:
        return self._cmd
    
    @property
    def check(self) -> bool:
        return self._check
    
    @property
    def success(self) -> bool:
        return self._proc.returncode == 0 if self._proc is not None else False

    @property
    def proc(self):
        return self._proc

    @property
    def stderr(self) -> str:
        return self._proc.stderr if self._proc is not None else ''
    @property
    def stdout(self) -> str:
        return self._proc.stdout if self._proc is not None else ''

    def append_cmd(self, cmd: Union[List, str]) -> None:
        if isinstance(cmd, str):
            cmd = cmd.split()
        self._cmd += cmd

    def run(self, ignore_opreated: bool=False) -> None:
        if self._operated:
            if ignore_opreated:
                logging.warning(f'{self} has already been operated, running again.')
            else:
                logging.error(f'{self} has already been operated, skipping.')
                return
        self._operated = True
        kwargs= {k: v for k, v in self._kwargs.items() if k not in self._USED_KWARGS}
        self._proc = run(self._cmd, capture_output=True, encoding='utf-8', shell=False, errors='strict', input=self._stdin, check=self._check, **kwargs)



