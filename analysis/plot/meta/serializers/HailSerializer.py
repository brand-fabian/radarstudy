import os
import typing

import hail as hl

from .Serializers import FILE_TYPE, Serializer


class MatrixTableSerializer(Serializer):
    @classmethod
    def provides(cls):
        return FILE_TYPE.MATRIX_TABLE

    @classmethod
    def read(cls, path: str):
        return hl.read_matrix_table(os.path.abspath(path))

    @classmethod
    def write(cls, data: typing.Any, path: str):
        data.write(os.path.abspath(path), overwrite=True)


class TableSerializer(Serializer):
    @classmethod
    def provides(cls):
        return FILE_TYPE.HAIL_TABLE

    @classmethod
    def read(cls, path: str):
        return hl.read_table(os.path.abspath(path))

    @classmethod
    def write(cls, data: typing.Any, path: str):
        data.write(os.path.abspath(path), overwrite=True)
