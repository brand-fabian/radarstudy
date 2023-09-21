import logging
import os
import typing

import pandas
import plotly

from .Serializers import FILE_TYPE, Serializer

logger = logging.getLogger(__name__)


class PlotSerializerException(Exception):
    pass


class PlotSerializer(Serializer):
    @classmethod
    def provides(cls):
        return FILE_TYPE.PLOT

    @classmethod
    def read(cls, path: str):
        return None

    @classmethod
    def write(cls, data: typing.Any, path: str):
        if len(data) != 2:
            raise PlotSerializerException("missing input data to serialize plot.")
        fig, table = data
        base_path = os.path.abspath(path)
        table.to_csv("{}.csv".format(base_path), index=False)
        table.to_latex("{}.tex".format(base_path), index=False)
        fig.write_html("{}.html".format(base_path))
        fig.write_image("{}.png".format(base_path))
        fig.write_image("{}.svg".format(base_path))
        fig.write_image("{}.pdf".format(base_path))


class RegressionSerializer(Serializer):
    @classmethod
    def provides(cls):
        return FILE_TYPE.REGRESSION_GLM

    @classmethod
    def read(cls, path: str):
        return None

    @classmethod
    def write(cls, data: typing.Any, path: str):
        table, fit = data
        summary = fit.summary()
        base_path = os.path.abspath(path)
        logger.debug("writing glm output for {}".format(base_path))
        with open("{}.glm.txt".format(base_path), "w") as out_f:
            out_f.write(summary.as_text())

        with open("{}.glm.tex".format(base_path), "w") as out_f:
            out_f.write(summary.as_latex())

        table.to_csv("{}.csv".format(base_path), index=False)
