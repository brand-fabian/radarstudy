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


class MsdnAnalysisFactory(AnalysisFactory):
    @classmethod
    def get(cls, matching_order={}):
        class MsdnAnalysis:
            @staticmethod
            def _get_msdns(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                what: typing.Literal["pairs", "clusters"],
            ):
                group_size = {
                    "pairs": lambda g: len(g),
                    "clusters": lambda g: len(g.msdn_id.unique()),
                }

                col_name = {
                    "pairs": "msdn_pairs",
                    "clusters": "msdn_clusters",
                }

                data = []
                for meta, group in msdn_df.groupby(["s", "cohort"]):
                    data.append([*meta, group_size[what](group)])

                fig_df = pandas.DataFrame(
                    data, columns=["sample", "cohort", col_name[what]]
                )
                fig_df = fill_missing(fig_df, {col_name[what]: 0}, metadata=age_data)

                return fig_df

            @staticmethod
            @AnalysisFactory.STATS_FACTORY.register(
                name="msdn_pairs", index=1, value_col="msdn_pairs", order=matching_order
            )
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="msdn_pairs", output_type=FILE_TYPE.PLOT
            )
            def plot_msdn_pairs(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                lang=AnalysisFactory.LANGUAGE,
            ):
                fig_df = MsdnAnalysis._get_msdns(msdn_df, age_data, "pairs")

                fig = px.box(
                    fig_df, y="msdn_pairs", color="cohort", labels=get_labels(lang=lang)
                )
                update_legend(fig)
                return fig, fig_df

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR,
                name="glm_msdn_pairs",
                output_type=FILE_TYPE.REGRESSION_GLM,
            )
            def regression_msdn_pairs(
                msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
            ):
                fig_df = MsdnAnalysis._get_msdns(msdn_df, age_data, "pairs")

                l = LEVELS_BASE
                reg_fit = smf.glm(
                    formula="msdn_pairs ~ C(cohort, levels=l)",
                    data=fig_df,
                    family=sm.families.NegativeBinomial(),
                ).fit()

                return fig_df, reg_fit

            @staticmethod
            @AnalysisFactory.STATS_FACTORY.register(
                name="msdn_clusters",
                index=1,
                value_col="msdn_clusters",
                order=matching_order,
            )
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="msdn_clusters", output_type=FILE_TYPE.PLOT
            )
            def plot_msdn_clusters(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                lang=AnalysisFactory.LANGUAGE,
            ):
                fig_df = MsdnAnalysis._get_msdns(msdn_df, age_data, "clusters")

                fig = px.box(
                    fig_df,
                    y="msdn_clusters",
                    color="cohort",
                    labels=get_labels(lang=lang),
                )
                update_legend(fig)
                return fig, fig_df

            @staticmethod
            @AnalysisFactory.STATS_FACTORY.register(
                name="msdn_cluster_violin",
                index=1,
                value_col="msdn_clusters",
                order=matching_order,
            )
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="msdn_clusters_violin", output_type=FILE_TYPE.PLOT
            )
            def plot_msdn_clusters_violin(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                lang=AnalysisFactory.LANGUAGE,
            ):
                fig_df = MsdnAnalysis._get_msdns(msdn_df, age_data, "clusters")

                fig = px.violin(
                    fig_df,
                    y="msdn_clusters",
                    points="all",
                    box=True,
                    color="cohort",
                    labels=get_labels(lang=lang),
                )
                fig.update_traces(spanmode="hard", bandwidth=0.6)
                update_legend(fig)

                return fig, fig_df

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR,
                name="glm_msdn_clusters",
                output_type=FILE_TYPE.REGRESSION_GLM,
            )
            def regression_msdn_clusters(
                msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
            ):
                fig_df = MsdnAnalysis._get_msdns(msdn_df, age_data, "clusters")

                l = LEVELS_BASE
                reg_fit = smf.glm(
                    formula="msdn_clusters ~ C(cohort, levels=l)",
                    data=fig_df,
                    family=sm.families.NegativeBinomial(),
                ).fit()

                return fig_df, reg_fit

        return MsdnAnalysis
