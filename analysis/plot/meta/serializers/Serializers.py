import typing
from abc import ABC, abstractmethod
from enum import Enum


class FILE_TYPE(Enum):
    MATRIX_TABLE = "MATRIX_TABLE"
    HAIL_TABLE = "HAIL_TABLE"
    DATAFRAME = "DATAFRAME"
    PICKLE = "PICKLE"
    PLOT = "PLOT"
    REGRESSION_GLM = "REGRESSION_GLM"
    DATAFRAME_TEX = "DATAFRAME_TEX"


class Serializer(ABC):
    def __init__(self, *args, **kwargs):
        pass

    @classmethod
    @abstractmethod
    def read(cls, path: str) -> typing.Any:
        """Read the data at path into a valid instance of the provided type."""
        pass

    @classmethod
    @abstractmethod
    def write(cls, data: typing.Any, path: str):
        """Serialize the data into a file or directory at path."""
        pass

    @classmethod
    @abstractmethod
    def provides(cls) -> FILE_TYPE:
        """Functionality provided by class."""
        return None
