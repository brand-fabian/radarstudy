import typing

import numpy as np
import pandas
import plotly.express as px
import statsmodels.api as sm
import statsmodels.formula.api as smf

from meta.checkpoint import NamedCheckpoint
from meta.config import CACHE_DIR, LEVELS_BASE
from meta.helpers import fill_missing
from meta.plots import AnalysisFactory
from meta.plotting import get_label, get_labels, merge_figures, update_legend, DEFAULT_PLOTLY_STYLES
from meta.serializers import FILE_TYPE


class MsdnByDosageAnalysisFactory(AnalysisFactory):
    @classmethod
    def get(cls, matchin_order={}):
        class MsdnByDosageAnalysis:
            @staticmethod
            def _get_msdns_by_dosage(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                what: typing.Literal["pairs", "clusters"],
            ):
                group_size = {
                    "pairs": lambda g: len(g),
                    "clusters": lambda g: len(g.msdn_id.unique()),
                }

                age_data_index = age_data.set_index("rid", drop=False)
                data = []
                for meta, group in msdn_df.groupby(["s", "cohort"]):
                    exposure = (
                        age_data_index.loc[meta[0]].sum_dosage_gonads
                        if meta[0] in age_data_index.index
                        else 0
                    )
                    exposure_status = (
                        age_data_index.loc[meta[0]].exposure_status
                        if meta[0] in age_data_index.index
                        and pandas.notna(age_data_index.loc[meta[0]].exposure_status)
                        else "UNKNOWN"
                    )
                    data.append(
                        [
                            *meta,
                            exposure,
                            exposure_status,
                            "{}::{}".format(exposure_status, meta[1]),
                            group_size[what](group),
                        ]
                    )

                fig_df = pandas.DataFrame(
                    data,
                    columns=[
                        "sample",
                        "cohort",
                        "sum_dosage_gonads",
                        "exposure_status",
                        "exposure_cohort",
                        "count",
                    ],
                )
                fig_df = fill_missing(
                    fig_df,
                    col_defaults={
                        "exposure_status": "NO_EXPOSURE",
                        "exposure_cohort": "UNKNOWN",
                        "count": 0,
                        "sum_dosage_gonads": 0,
                    },
                    metadata=age_data,
                )
                fig_df["exposure_cohort"] = (
                    fig_df["cohort"] + "::" + fig_df["exposure_status"]
                )
                return fig_df

            @staticmethod
            def _regression_msdns_by_dosage(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                what: typing.Literal["pairs", "clusters"],
            ):
                fig_df = MsdnByDosageAnalysis._get_msdns_by_dosage(
                    msdn_df,
                    age_data,
                    what=what,
                )

                reg_fit = smf.glm(
                    formula="count ~ sum_dosage_gonads",
                    data=fig_df,
                    family=sm.families.NegativeBinomial(),
                ).fit()

                return fig_df, reg_fit

            @staticmethod
            def _regression_msdns_by_dosage_and_cohort(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                what: typing.Literal["pairs", "clusters"],
            ):
                fig_df = MsdnByDosageAnalysis._get_msdns_by_dosage(
                    msdn_df, age_data, what=what
                )

                l = LEVELS_BASE
                reg_fit = smf.glm(
                    formula="count ~ sum_dosage_gonads : C(cohort, levels=l)",
                    data=fig_df,
                    family=sm.families.NegativeBinomial(),
                ).fit()

                return fig_df, reg_fit

            @staticmethod
            def _plot_msdns_by_dosage(
                fig_df: pandas.DataFrame,
                fit,
                lang=AnalysisFactory.LANGUAGE,
                y_label: typing.Literal["msdn_pairs", "msdn_clusters"] = "msdn_pairs",
            ):
                desc = fig_df.groupby("cohort").describe()

                dosage_values = {
                    cohort: np.linspace(
                        desc["sum_dosage_gonads"]["min"].loc[cohort],
                        desc["sum_dosage_gonads"]["max"].loc[cohort],
                        100,
                    )
                    for cohort in ["RADAR", "CRU"]
                }

                predict_df = pandas.concat(
                    list(
                        map(
                            lambda x: pandas.DataFrame(
                                map(lambda y: [x[0], y], x[1]),
                                columns=["cohort", "sum_dosage_gonads"],
                            ),
                            dosage_values.items(),
                        )
                    )
                )

                predict_df["count"] = fit.predict(predict_df)

                line = px.line(
                    predict_df,
                    x="sum_dosage_gonads",
                    y="count",
                    color="cohort",
                    range_x=(0, 1500),
                    category_orders={"cohort": ["INOVA", "CRU", "RADAR", "PILOT"]},
                    color_discrete_sequence=[
                        px.colors.qualitative.G10[0],
                        px.colors.qualitative.G10[2],
                        px.colors.qualitative.G10[1],
                        px.colors.qualitative.G10[3],
                    ]
                )
                line.data[0].name = "GLM Radar"
                line.data[1].name = "GLM CRU"

                fig = px.strip(
                    fig_df,
                    x="sum_dosage_gonads",
                    y="count",
                    color="cohort",
                    range_x=(0, 1500),
                    orientation="h",
                    stripmode="overlay",
                    # trendline="ols",
                    category_orders={"cohort": ["INOVA", "CRU", "RADAR", "PILOT"]},
                    labels=get_labels(
                        lang=lang, labels={"count": get_label(y_label, lang=lang)}
                    ),
                    color_discrete_sequence=[
                        px.colors.qualitative.G10[0],
                        px.colors.qualitative.G10[2],
                        px.colors.qualitative.G10[1],
                        px.colors.qualitative.G10[3],
                    ]
                )
                fig.update_traces(jitter=1, marker={ "size": 8 })
                for trace, symbol in zip(fig.data, [0, 2, 4, 3]):
                    trace["marker"]["symbol"] = symbol
                fig = merge_figures(fig, line)
                update_legend(fig)
                return fig, fig_df

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="pairs_by_dosage", output_type=FILE_TYPE.PLOT
            )
            def plot_msdn_pairs_by_dosage(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                lang=AnalysisFactory.LANGUAGE,
            ):
                fig_df = MsdnByDosageAnalysis._get_msdns_by_dosage(
                    msdn_df,
                    age_data,
                    what="pairs",
                )
                fig_df = fig_df[fig_df.sum_dosage_gonads > 0]

                _, fit = MsdnByDosageAnalysis._regression_msdns_by_dosage(
                    msdn_df, age_data, what="pairs"
                )

                return MsdnByDosageAnalysis._plot_msdns_by_dosage(
                    fig_df,
                    fit,
                    lang=lang,
                    y_label="msdn_pairs",
                )

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR,
                name="pairs_by_dosage_and_cohort",
                output_type=FILE_TYPE.PLOT,
            )
            def plot_msdn_pairs_by_dosage_and_cohort(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                lang=AnalysisFactory.LANGUAGE,
            ):
                fig_df = MsdnByDosageAnalysis._get_msdns_by_dosage(
                    msdn_df,
                    age_data,
                    what="pairs",
                )
                fig_df = fig_df[fig_df.sum_dosage_gonads > 0]

                _, fit = MsdnByDosageAnalysis._regression_msdns_by_dosage_and_cohort(
                    msdn_df, age_data, what="pairs"
                )

                return MsdnByDosageAnalysis._plot_msdns_by_dosage(
                    fig_df,
                    fit,
                    lang=lang,
                    y_label="msdn_pairs",
                )

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR,
                name="glm_pairs_by_dosage",
                output_type=FILE_TYPE.REGRESSION_GLM,
            )
            def regression_msdn_pairs_by_dosage(
                msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
            ):
                return MsdnByDosageAnalysis._regression_msdns_by_dosage(
                    msdn_df,
                    age_data,
                    what="pairs",
                )

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR,
                name="glm_pairs_by_dosage_and_cohort",
                output_type=FILE_TYPE.REGRESSION_GLM,
            )
            def regression_msdn_pairs_by_dosage_and_cohort(
                msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
            ):
                return MsdnByDosageAnalysis._regression_msdns_by_dosage_and_cohort(
                    msdn_df,
                    age_data,
                    what="pairs",
                )

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="clusters_by_dosage", output_type=FILE_TYPE.PLOT
            )
            def plot_msdn_clusters_by_dosage(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                lang=AnalysisFactory.LANGUAGE,
            ):
                fig_df = MsdnByDosageAnalysis._get_msdns_by_dosage(
                    msdn_df, age_data, what="clusters"
                )
                fig_df = fig_df[fig_df.sum_dosage_gonads > 0]

                _, fit = MsdnByDosageAnalysis._regression_msdns_by_dosage(
                    msdn_df,
                    age_data,
                    what="clusters",
                )

                return MsdnByDosageAnalysis._plot_msdns_by_dosage(
                    fig_df,
                    fit,
                    lang=lang,
                    y_label="msdn_clusters",
                )

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR,
                name="clusters_by_dosage_and_cohort",
                output_type=FILE_TYPE.PLOT,
            )
            def plot_msdn_clusters_by_dosage_and_cohort(
                msdn_df: pandas.DataFrame,
                age_data: pandas.DataFrame,
                lang=AnalysisFactory.LANGUAGE,
            ):
                fig_df = MsdnByDosageAnalysis._get_msdns_by_dosage(
                    msdn_df, age_data, what="clusters"
                )
                fig_df = fig_df[fig_df.sum_dosage_gonads > 0]

                _, fit = MsdnByDosageAnalysis._regression_msdns_by_dosage_and_cohort(
                    msdn_df,
                    age_data,
                    what="clusters",
                )

                return MsdnByDosageAnalysis._plot_msdns_by_dosage(
                    fig_df,
                    fit,
                    lang=lang,
                    y_label="msdn_clusters",
                )

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR,
                name="glm_clusters_by_dosage",
                output_type=FILE_TYPE.REGRESSION_GLM,
            )
            def regression_msdn_clusters_by_dosage(
                msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
            ):
                return MsdnByDosageAnalysis._regression_msdns_by_dosage(
                    msdn_df,
                    age_data,
                    what="clusters",
                )

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR,
                name="glm_clusters_by_dosage_and_cohort",
                output_type=FILE_TYPE.REGRESSION_GLM,
            )
            def regression_msdn_clusters_by_dosage_and_cohort(
                msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
            ):
                return MsdnByDosageAnalysis._regression_msdns_by_dosage_and_cohort(
                    msdn_df,
                    age_data,
                    what="clusters",
                )

        return MsdnByDosageAnalysis
