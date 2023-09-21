import hashlib
import inspect
import logging
import os
from glob import glob
from itertools import chain

import pandas

from meta.checkpoint import NamedCheckpoint
from meta.config import CACHE_DIR
from meta.plots import AnalysisFactory
from meta.serializers import FILE_TYPE

logger = logging.getLogger(__name__)
FOLDER = "simulations"
VERSION = ".version"


class SimulationAnalysisFactory(AnalysisFactory):
    @classmethod
    def version(cls, ppv: float):
        return hashlib.md5(
            (inspect.getsource(cls) + ";ppv={ppv}".format(ppv=ppv)).encode("utf-8")
        ).hexdigest()

    @classmethod
    def check_version(cls, ppv: float):
        if os.path.isfile(os.path.join(CACHE_DIR, FOLDER, VERSION)):
            old = None
            with open(os.path.join(CACHE_DIR, FOLDER, VERSION), "r") as vers_f:
                old = vers_f.read().strip()

            if old is None:
                return

            if old != cls.version(ppv):
                for f in glob(os.path.join(CACHE_DIR, FOLDER, "*")):
                    try:
                        os.remove(f)
                    except OSError as err:
                        logger.error("error deleting {}: {}".format(f, err))

        with open(os.path.join(CACHE_DIR, FOLDER, VERSION), "w") as vers_f:
            vers_f.write(cls.version(ppv))

    @classmethod
    def get(cls, iteration: int, ppv: float, matching_order={}):
        os.makedirs(os.path.join(CACHE_DIR, FOLDER), exist_ok=True)
        cls.check_version(ppv)

        class SimulationAnalysis:
            @staticmethod
            def _simulate_ppv(
                msdn_df: pandas.DataFrame, ppv: float = ppv, cohort_col: str = "cohort"
            ):
                return pandas.concat(
                    [
                        msdn_df[msdn_df[cohort_col] == c].sample(frac=ppv)
                        for c in msdn_df[cohort_col].unique()
                    ]
                )

            @staticmethod
            @AnalysisFactory.STATS_FACTORY.register(
                name="{}/simulation_{}".format(FOLDER, iteration),
                index=0,
                group_col="cohort",
                value_col="msdn_clusters",
                order=matching_order,
            )
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR,
                name="{}/simulation_{}.pickle".format(FOLDER, iteration),
                output_type=FILE_TYPE.PICKLE,
            )
            def simulation(msdn_df: pandas.DataFrame, age_data: pandas.DataFrame):
                df = SimulationAnalysis._simulate_ppv(msdn_df)

                MsdnAnalysis = AnalysisFactory.instance(
                    "MsdnAnalysisFactory", matching_order=matching_order
                )

                _, reg_fit_pairs = MsdnAnalysis.regression_msdn_pairs(df, age_data)
                fig_df, reg_fit_clusters = MsdnAnalysis.regression_msdn_clusters(
                    df, age_data
                )

                data = []
                for row, p_value, group in chain(
                    zip(
                        reg_fit_pairs.params.items(),
                        reg_fit_pairs.pvalues,
                        ["msdn_pairs"] * len(reg_fit_pairs.params),
                    ),
                    zip(
                        reg_fit_clusters.params.items(),
                        reg_fit_clusters.pvalues,
                        ["msdn_clusters"] * len(reg_fit_clusters.params),
                    ),
                ):
                    data.append([iteration, group, *row, p_value])

                return fig_df, pandas.DataFrame(
                    data,
                    columns=["iteration", "group", "coef", "coef_value", "p_value"],
                )

        return SimulationAnalysis
