import logging
from abc import abstractmethod
from typing import Any

logging.basicConfig()
logger = logging.getLogger(__name__)


class Runner:
    """Running a command or task and capture its output when we ask for it.

    """
    def __init__(self, cmd: Any, operated: bool):
        self._cmd = cmd
        self._operated = operated

    @abstractmethod
    def __repr__(self) -> str:
        pass

    @abstractmethod
    def run(self, *args, **kwargs):
        pass
    
    @property
    @abstractmethod
    def success(self) -> bool:
        pass

    @property
    def operated(self) -> bool:
        return self._operated

    @property
    def cmd(self) -> str:
        return self._cmd
