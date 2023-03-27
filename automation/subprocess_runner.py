from subprocess import run
from typing import List, Union
import logging
from runner import Runner
from pathlib import Path

logger = logging.getLogger(__name__)


class SubprocessRunner(Runner):
    _USED_KWARGS = ['capture_output', 'encoding', 'shell', 'input', 'errors', 'check']

    def __init__(self, cmd: Union[List, str], stdin: str="", check: bool=False, **kwargs):
        """Run a subprocess and capture its output when we ask for it.

        :param cmd: command to be executed
        :type cmd: Union[List, str]
        :param stdin: data to be passed to stdin, defaults to ""
        :type stdin: str, optional
        :param check: check if the subporocess return 0, defaults to False
        :type check: bool, optional
        :param kwargs: keyword arguments to be passed to subprocess.run
        :type kwargs: dict
        """
        if isinstance(cmd, str):
            cmd = cmd.split()
        self._stdin = stdin
        self._check = check
        self._proc = None
        self._kwargs = kwargs
        super().__init__(cmd, operated=False)
    
    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._cmd}, {self._stdin}, {self._kwargs})'

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
        if self._proc is None:
            raise RuntimeError('Subprocess has not been run yet.')
        return self._proc.stderr
    
    @property
    def stdout(self) -> str:
        if self._proc is None:
            raise RuntimeError('Subprocess has not been run yet.')
        return self._proc.stdout
    
    def append_cmd(self, cmd: Union[List, str]) -> None:
        if isinstance(cmd, str):
            cmd = cmd.split()
        self._cmd += cmd

    def write(self, saveto: Union[str, Path]) -> None:
        if self._proc is None:
            raise RuntimeError('Subprocess has not been run yet.')
        with open(saveto, 'w') as f:
            f.write(self._proc.stdout)

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



