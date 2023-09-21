import typing

import pandas
import plotly.express as px
import numpy as np

from meta.checkpoint import NamedCheckpoint
from meta.config import CACHE_DIR, LEVELS_BASE, EXPOSED_COHORTS
from meta.helpers import fill_missing
from meta.plots import AnalysisFactory
from meta.plotting import get_labels, update_legend, get_label
from meta.serializers import FILE_TYPE


class MetadataAnalysisFactory(AnalysisFactory):
    @classmethod
    def get(cls, matching_order={}):
        class MetadataAnalysis:
            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR,
                name="age_distribution",
                output_type=FILE_TYPE.PLOT,
            )
            def plot_age_distribution(age_data, lang=AnalysisFactory.LANGUAGE):
                fig = px.scatter(
                    age_data,
                    x="father_age",
                    y="mother_age",
                    color="cohort",
                    symbol="cohort",
                    labels=get_labels(lang=lang),
                )
                update_legend(fig)
                return fig, age_data

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="service_duration", output_type=FILE_TYPE.PLOT
            )
            def plot_service_duration(age_data, lang=AnalysisFactory.LANGUAGE):
                age_data["service_duration"] = age_data["service_end"] - age_data["service_begin"]
                age_data["service_duration"] = age_data["service_duration"].clip(lower=1) + 0.5
                max_svc_dur = int(age_data[age_data.cohort == "RADAR"]["service_duration"].max())
                age_data = age_data[age_data.cohort == "RADAR"]
                fig = px.histogram(
                    age_data,
                    x="service_duration",
                    color="cohort",
                    nbins=int(max_svc_dur / 1.5),
                    text_auto=True,
                    category_orders={"service_duration": list(map(str, range(max_svc_dur)))},
                    color_discrete_sequence=px.colors.qualitative.G10[1:],
                    range_x=(0, max_svc_dur + 1),
                    labels=get_labels(lang=lang),
                )
                update_legend(fig)
                return fig, age_data

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="dosage_by_service_duration", output_type=FILE_TYPE.PLOT
            )
            def plot_dosage_by_service_duration(age_data, lang=AnalysisFactory.LANGUAGE):
                age_data["service_duration"] = age_data["service_end"] - age_data["service_begin"]
                age_data = age_data[age_data.cohort == "RADAR"]
                fig = px.scatter(
                    age_data,
                    y="sum_dosage_gonads",
                    x="service_duration",
                    color="cohort",
                    symbol="cohort",
                    log_y=True,
                    trendline="ols",
                    labels=get_labels(lang=lang),
                )
                update_legend(fig)
                return fig, age_data

            @staticmethod
            def _get_dosage_histogram(
                df: pandas.DataFrame,
                dosage_col: str = "sum_dosage_gonads",
                nbins: int = 25,
            ):
                data = []

                real_min = min(filter(lambda x: x > 0, df[dosage_col]))
                bins = np.append(
                    np.insert(np.linspace(real_min, df[dosage_col].max(), nbins), 0, 0),
                    df[dosage_col].max(),
                )

                for meta, group in df.groupby("cohort"):
                    for value, bin in zip(*np.histogram(group[dosage_col], bins=bins)):
                        data.append([
                            meta,
                            bin,
                            value,
                        ])
                return pandas.DataFrame(data, columns=["cohort", "sum_dosage_gonads", "count"])

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="dosage_histogram", output_type=FILE_TYPE.PLOT
            )
            def plot_dosage_histogram(
                age_data,
                lang=AnalysisFactory.LANGUAGE,
                dosage_col: str = "sum_dosage_gonads",
                exposed_cohorts: typing.List[str] = EXPOSED_COHORTS,
                nbins: int = 25,
            ):
                control_cohorts = set(age_data.cohort.unique()) - set(exposed_cohorts)
                fig_df = MetadataAnalysis._get_dosage_histogram(age_data[~age_data.cohort.isin(control_cohorts)], nbins=nbins, dosage_col=dosage_col)
                fig = px.bar(
                    fig_df,
                    x="sum_dosage_gonads",
                    y="count",
                    color="cohort",
                    barmode="group",
                    labels=get_labels(lang=lang, labels={"count": get_label("no_children", lang=lang)}),
                )
                update_legend(fig)

                # Add mean values for each exposed cohorts as vertical line
                for cohort in exposed_cohorts:
                    mean_exposure = age_data[age_data.cohort == "RADAR"].sum_dosage_gonads.mean()
                    fig.add_vline(
                        x=mean_exposure,
                        line_color="black",
                        line_width=2,
                        annotation_text=get_label("mean_exposure_annotation", lang=lang).format(cohort.title(), mean_exposure),
                        annotation_textangle=270,
                        annotation_position="top left",
                        opacity=1,
                    )
                return fig, fig_df

            @staticmethod
            def _get_parental_age_histogram(
                age_data: pandas.DataFrame,
                age_col: str = "father_age",
                bins: np.ndarray = np.arange(10, 60, 5),
            ) -> pandas.DataFrame:
                data = []

                for meta, group in age_data.groupby("cohort"):
                    for value, bin in zip(*np.histogram(group[age_col], bins=bins)):
                        data.append([
                            meta,
                            bin,
                            value / len(group)
                        ])
                return pandas.DataFrame(data, columns=["cohort", "age", "count"])

            @staticmethod
            @NamedCheckpoint.checkpoint(
                prefix=CACHE_DIR, name="parental_age_histogram", output_type=FILE_TYPE.PLOT
            )
            def plot_parental_age_histogram(
                age_data: pandas.DataFrame,
                lang=AnalysisFactory.LANGUAGE,
            ):
                father_df = MetadataAnalysis._get_parental_age_histogram(age_data)
                mother_df = MetadataAnalysis._get_parental_age_histogram(age_data, age_col="mother_age")

                fig = px.bar(
                    father_df,
                    x="age",
                    y="count",
                    color="cohort",
                    barmode="group",
                    labels=get_labels(lang=lang, labels={"count": get_label("parental_histogram_y", lang=lang)}),
                )

                mat_age = px.bar(
                    mother_df,
                    x="age",
                    y="count",
                    color="cohort",
                    barmode="group",
                )

                for trace in mat_age.data:
                    trace.y = trace.y * -1
                    trace.showlegend = False
                    fig.add_trace(trace)

                update_legend(fig)

                fig.update_layout(
                    yaxis={
                        "tickformat": ".0%",
                        "tickvals": np.arange(-1, 1, 0.1),
                        "ticktext": list(map(lambda x: str(abs(x)) + "%", np.arange(-100, 100, 10))),
                    }
                )
                return fig, pandas.concat([father_df, mother_df])
        return MetadataAnalysis