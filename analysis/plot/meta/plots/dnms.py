import typing

import pandas
import plotly.express as px
import statsmodels.api as sm
import statsmodels.formula.api as smf

from meta.checkpoint import NamedCheckpoint
from meta.config import CACHE_DIR, LEVELS_BASE
from meta.helpers import fill_missing
from meta.plots import AnalysisFactory
from meta.plotting import get_labels, update_legend
from meta.serializers import FILE_TYPE


class DnmAnalysisFactory(AnalysisFactory):
    @classmethod
    def get(cls, matching_order={}):
        class DnmAnalysis:
            @staticmethod
            def _get_dnms(
                dnm_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
            ):
                dnm_data = []
                for sample, group in dnm_df.groupby(["s", "cohort"]):
                    dnm_data.append(
                        [
                            *sample,
                            len(group),
                        ]
                    )
                dnm_per_sample_df = pandas.DataFrame(
                    dnm_data, columns=["sample", "cohort", "count"]
                )
                return dnm_per_sample_df

            @staticmethod
            @AnalysisFactory.STATS_FACTORY.register(
                name="dnms", index=1, order=matching_order
            )
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="dnms", output_type=FILE_TYPE.PLOT
            )
            def plot_dnms(
                dnm_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                lang=AnalysisFactory.LANGUAGE,
            ):
                fig_df = DnmAnalysis._get_dnms(dnm_df, age_data)

                fig = px.box(
                    fig_df,
                    y="count",
                    color="cohort",
                    hover_data=["sample"],
                    labels=get_labels(lang=lang, labels={"count": "DNMs"}),
                )
                update_legend(fig)
                return fig, fig_df

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="glm_dnms", output_type=FILE_TYPE.REGRESSION_GLM
            )
            def regression_dnms(dnm_df: pandas.DataFrame, age_data: pandas.DataFrame):
                fig_df = DnmAnalysis._get_dnms(dnm_df, age_data)
                l = LEVELS_BASE

                reg_fit = smf.glm(
                    formula="count ~ C(cohort, levels=l)",
                    data=fig_df,
                    family=sm.families.NegativeBinomial(),
                ).fit()

                return fig_df, reg_fit

        return DnmAnalysis
