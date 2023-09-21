import os
import pickle
import typing

import pandas

from .Serializers import FILE_TYPE, Serializer


class DataFrameSerializer(Serializer):
    @classmethod
    def provides(cls):
        return FILE_TYPE.DATAFRAME

    @classmethod
    def read(cls, path: str):
        return pandas.read_csv(os.path.abspath(path))

    @classmethod
    def write(cls, data: typing.Any, path: str):
        data.to_csv(os.path.abspath(path), index=False)


class PickleSerializer(Serializer):
    @classmethod
    def provides(cls):
        return FILE_TYPE.PICKLE

    @classmethod
    def read(cls, path: str):
        with open(os.path.abspath(path), "rb") as pickle_f:
            return pickle.load(pickle_f)

    @classmethod
    def write(cls, data: typing.Any, path: str):
        with open(os.path.abspath(path), "wb") as pickle_f:
            pickle.dump(data, pickle_f)


class DataFrameSerializerTex(DataFrameSerializer):
    @classmethod
    def provides(cls):
        return FILE_TYPE.DATAFRAME_TEX

    @classmethod
    def read(cls, path: str):
        return super(DataFrameSerializerTex, cls).read(path)

    @classmethod
    def write(cls, data: pandas.DataFrame, path: str):
        tex_path = os.path.abspath(os.path.splitext(path)[0] + ".tex")
        data.to_latex(tex_path, index=False)
        super(DataFrameSerializerTex, cls).write(data, path)
