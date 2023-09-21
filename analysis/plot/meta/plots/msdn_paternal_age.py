import logging
import os
import typing

import numpy as np
import pandas
import plotly.express as px
import statsmodels.api as sm
import statsmodels.formula.api as smf

from meta.checkpoint import NamedCheckpoint
from meta.helpers import fill_missing
from meta.plotting import get_label, get_labels, merge_figures, update_legend
from meta.serializers import FILE_TYPE
from meta.config import LEVELS_BASE, CATEGORY_ORDER_BASE

CACHE_DIR = os.path.abspath(os.getenv("CACHE_DIR", "./output"))
logger = logging.getLogger(__name__)


def _get_x_by_paternal_age(
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
        data.append(
            [
                *meta,
                age_data_index.loc[meta[0]].father_age
                if meta[0] in age_data_index.index
                else 0,
                group_size[what](group),
            ]
        )
    fig_df = pandas.DataFrame(data, columns=["sample", "cohort", "father_age", "count"])
    fig_df = fill_missing(fig_df, {"father_age": 0, "count": 0}, metadata=age_data)
    fig_df.father_age = list(
        map(
            lambda x: age_data_index.loc[x].father_age
            if x in age_data_index.index
            else 0,
            fig_df["sample"],
        )
    )
    return fig_df


def _get_regression_fit(fig_df: pandas.DataFrame, fit):
    # Alternatively, retrieve cohorts from the df: fig_df.cohorts.unique(),
    # but this does not guarantee the correct ordering of the regressions in
    # the final table
    cohorts = list(set(CATEGORY_ORDER_BASE) & set(fig_df.cohort.unique()))
    desc = fig_df.groupby("cohort").describe()

    ages = {
        cohort.lower(): np.linspace(
            desc["father_age"]["min"].loc[cohort],
            desc["father_age"]["max"].loc[cohort],
            100,
        )
        for cohort in cohorts
    }

    predict_df = pandas.concat(
        [
            pandas.DataFrame(
                map(lambda x: [cohort, x], ages[cohort.lower()]),
                columns=["cohort", "father_age"],
            )
            for cohort in cohorts
        ]
    )
    predict_df["count"] = fit.predict(predict_df)
    return predict_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="pairs_by_paternal_age", output_type=FILE_TYPE.PLOT
)
def plot_msdn_pairs_by_paternal_age(
    msdn_df: pandas.DataFrame,
    age_data: pandas.DataFrame,
    lang: typing.Literal["en", "de"] = "en",
):
    fig_df = _get_x_by_paternal_age(msdn_df, age_data, what="pairs")
    fig_df = fig_df[fig_df.father_age > 0]

    _, fit = regression_msdn_clusters_by_paternal_age(msdn_df, age_data)
    predict_df = _get_regression_fit(fig_df, fit)

    line = px.line(
        predict_df,
        x="father_age",
        y="count",
        color="cohort",
        category_orders={"cohort": CATEGORY_ORDER_BASE},
    )
    line.data[0].name = "GLM Inova"
    line.data[1].name = "GLM Radar"
    line.data[2].name = "GLM CRU"
    if len(line.data) > 3:
        line.data[3].name = "GLM Pilot"

    fig = px.scatter(
        fig_df,
        x="father_age",
        y="count",
        color="cohort",
        symbol="cohort",
        hover_data=["sample"],
        labels=get_labels(
            lang=lang, labels={"count": get_label("msdn_pairs", lang=lang)}
        ),
    )
    fig = merge_figures(fig, line)
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="clusters_by_paternal_age", output_type=FILE_TYPE.PLOT
)
def plot_msdn_clusters_by_paternal_age(
    msdn_df: pandas.DataFrame,
    age_data: pandas.DataFrame,
    lang: typing.Literal["en", "de"] = "en",
):
    fig_df = _get_x_by_paternal_age(msdn_df, age_data, what="clusters")
    fig_df = fig_df[fig_df.father_age > 0]

    _, fit = regression_msdn_clusters_by_paternal_age(msdn_df, age_data)
    predict_df = _get_regression_fit(fig_df, fit)

    line = px.line(
        predict_df,
        x="father_age",
        y="count",
        color="cohort",
        category_orders={"cohort": CATEGORY_ORDER_BASE},
    )
    line.data[0].name = "GLM Inova"
    line.data[1].name = "GLM Radar"
    line.data[2].name = "GLM CRU"
    if len(line.data) > 3:
        line.data[3].name = "GLM Pilot"

    fig = px.scatter(
        fig_df,
        x="father_age",
        y="count",
        color="cohort",
        symbol="cohort",
        hover_data=["sample"],
        labels=get_labels(
            lang=lang, labels={"count": get_label("msdn_clusters", lang=lang)}
        ),
    )
    fig = merge_figures(fig, line)
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="glm_clusters_by_paternal_age",
    output_type=FILE_TYPE.REGRESSION_GLM,
)
def regression_msdn_clusters_by_paternal_age(
    msdn_df: pandas.DataFrame,
    age_data: pandas.DataFrame,
    levels: typing.List[str] = LEVELS_BASE,
):
    fig_df = _get_x_by_paternal_age(msdn_df, age_data, what="clusters")
    reg_fit = smf.glm(
        formula="count ~ C(cohort, levels=levels) + father_age",
        data=fig_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return fig_df, reg_fit


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="glm_pairs_by_paternal_age",
    output_type=FILE_TYPE.REGRESSION_GLM,
)
def regression_msdn_pairs_by_paternal_age(
    msdn_df: pandas.DataFrame,
    age_data: pandas.DataFrame,
    levels: typing.List[str] = LEVELS_BASE,
):
    fig_df = _get_x_by_paternal_age(msdn_df, age_data, what="pairs")
    reg_fit = smf.glm(
        formula="count ~ C(cohort, levels=levels) + father_age",
        data=fig_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return fig_df, reg_fit
