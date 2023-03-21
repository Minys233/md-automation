import logging
from abc import ABC, abstractmethod

logging.basicConfig()
logger = logging.getLogger(__name__)


class Runner(ABC):
    """Running a command or task and capture its output when we ask for it.

    """
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __repr__(self) -> str:
        pass

    @abstractmethod
    def run(self, *args, **kwargs):
        pass

    @property
    @abstractmethod
    def operated(self) -> bool:
        pass

    @property
    @abstractmethod
    def cmd(self) -> str:
        pass
    
    @property
    @abstractmethod
    def success(self) -> bool:
        pass
