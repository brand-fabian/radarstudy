import argparse
import logging
import os
import pickle
import sys
import typing
from itertools import chain

import hail as hl
import numpy as np
import pandas
import patsy
import plotly.express as px
import statsmodels.api as sm
import statsmodels.formula.api as smf

from meta.arguments import (
    add_base_arguments,
    configure_logger,
    dump_call_information,
    hail_init,
    load_metadata_from_file,
)
from meta.checkpoint import NamedCheckpoint
from meta.config import CACHE_DIR, LEVELS_BASE, LOG_LEVEL, LEVELS_SUBGROUPS
from meta.helpers import downsample_age_match, fill_missing
from meta.loaders import GraphtyperLoader, MsdnLoader
from meta.plots import (
    AnalysisFactory,
    plot_msdn_clusters_by_paternal_age,
    plot_msdn_pairs_by_paternal_age,
    regression_msdn_clusters_by_paternal_age,
    regression_msdn_pairs_by_paternal_age,
)
from meta.plotting import (
    apply_styles,
    get_label,
    get_labels,
    merge_figures,
    update_legend,
)
from meta.serializers import FILE_TYPE
from meta.statistics import StatisticsFactory

# Statistics Factory used to get the main output statistics
StatsFactory = StatisticsFactory()

###############################################################################
# Logging                                                                     #
###############################################################################
logger = configure_logger()

###############################################################################
# Arguments                                                                   #
###############################################################################
parser = argparse.ArgumentParser(
    prog="PLOT-MSDN",
    description="Plot descriptive statistics for msdns.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "ht",
    metavar="HT",
    type=str,
    help="Input hail table (msdn_ht) containing multisite de novo mutations for both cohorts.",
)
parser = add_base_arguments(parser)

args = parser.parse_args()
logger.setLevel(LOG_LEVEL[args.verbose])
logger.debug("Arguments: {}".format(vars(args)))

logger.info("Writing cache output to {}".format(CACHE_DIR))
hail_init(args)

if not os.path.isdir(args.ht):
    logger.error(
        "Input path {} does not exist or is invalid matrix table".format(args.ht)
    )
    sys.exit(1)

if not args.load_cache:
    NamedCheckpoint.disable()

msdn_ht = MsdnLoader.load_input(args.ht)

if args.graphtyper is not None:
    graphtyper_df = GraphtyperLoader.load_graphtyper(
        args.graphtyper,
    )
    graphtyper_ht = GraphtyperLoader.load_graphtyper_hail(graphtyper_df)
    msdn_ht = GraphtyperLoader.annotate_graphtyper_msdn(
        graphtyper_ht,
        msdn_ht,
    )

msdn_df = MsdnLoader.load_dataframe(msdn_ht)
age_data: typing.Optional[pandas.DataFrame] = load_metadata_from_file(args)

logger.info(
    "Read {} msdns from {} samples.".format(
        len(msdn_df),
        len(age_data),
    )
)

if args.factor > 0:
    logger.info("Applying {}:1 downsampling".format(args.factor))
    msdn_df, age_data, matching = MsdnLoader.load_downsample_age_match(
        msdn_df, age_data, args.factor
    )
    matching_order = [
        *matching.keys(),
        *chain(*map(lambda x: list(x), matching.values())),
    ]
    control_matching = {s: r for r in matching.keys() for s in matching[r]}
    logger.info("Kept {} msdns from {} samples".format(len(msdn_df), len(age_data)))
else:
    matching = None
    matching_order = None
    control_matching = None

if not args.apply_control_matching:
    # Overwrite previous control_matching value based on option.
    control_matching = None

if args.apply_graphtyper_filter:
    msdn_all_df = msdn_df.copy(deep=True)
    msdn_df = MsdnLoader.load_graphtyper_filter(msdn_df)

# Apply plotly styling for all plots
apply_styles()


###############################################################################
# Checkpoints                                                                 #
###############################################################################
# These scriptlets/functions generate named output files in the cache dir that
# feature plots and datatables used to generate the data.
@StatsFactory.register(
    name="msdn_pairs", index=1, value_col="msdn_pairs", order=matching_order
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="msdn_pairs", output_type=FILE_TYPE.PLOT
)
def plot_msdn_pairs(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame, lang=args.language
):
    data = []
    for meta, group in msdn_df.groupby(["s", "cohort"]):
        data.append([*meta, len(group)])
    fig_df = pandas.DataFrame(data, columns=["sample", "cohort", "msdn_pairs"])
    fig_df = fill_missing(fig_df, {"msdn_pairs": 0}, metadata=age_data)

    fig = px.box(fig_df, y="msdn_pairs", color="cohort", labels=get_labels(lang=lang))
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="glm_msdn_pairs", output_type=FILE_TYPE.REGRESSION_GLM
)
def regression_msdn_pairs(msdn_df: pandas.DataFrame, age_data: pandas.DataFrame):
    data = []
    for meta, group in msdn_df.groupby(["s", "cohort"]):
        data.append([*meta, len(group)])
    fig_df = pandas.DataFrame(data, columns=["sample", "cohort", "msdn_pairs"])
    fig_df = fill_missing(fig_df, {"msdn_pairs": 0}, metadata=age_data)

    l = LEVELS_BASE
    reg_fit = smf.glm(
        formula="msdn_pairs ~ C(cohort, levels=l)",
        data=fig_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return fig_df, reg_fit


@StatsFactory.register(
    name="msdn_clusters", index=1, value_col="msdn_clusters", order=matching_order
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="msdn_clusters", output_type=FILE_TYPE.PLOT
)
def plot_msdn_clusters(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame, lang=args.language
):
    data = []
    for meta, group in msdn_df.groupby(["s", "cohort"]):
        data.append([*meta, len(group.msdn_id.unique())])
    fig_df = pandas.DataFrame(data, columns=["sample", "cohort", "msdn_clusters"])
    fig_df = fill_missing(fig_df, {"msdn_clusters": 0}, metadata=age_data)

    fig = px.box(
        fig_df, y="msdn_clusters", color="cohort", labels=get_labels(lang=lang)
    )
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="glm_msdn_clusters", output_type=FILE_TYPE.REGRESSION_GLM
)
def regression_msdn_clusters(msdn_df: pandas.DataFrame, age_data: pandas.DataFrame):
    data = []
    for meta, group in msdn_df.groupby(["s", "cohort"]):
        data.append([*meta, len(group.msdn_id.unique())])
    fig_df = pandas.DataFrame(data, columns=["sample", "cohort", "msdn_clusters"])
    fig_df = fill_missing(fig_df, {"msdn_clusters": 0}, metadata=age_data)

    l = LEVELS_BASE
    reg_fit = smf.glm(
        formula="msdn_clusters ~ C(cohort, levels=l)",
        data=fig_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return fig_df, reg_fit


@StatsFactory.register(
    name="cluster_size_histogram",
    index=1,
    value_col="cluster_size",
    order=matching_order,
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="cluster_size_histogram", output_type=FILE_TYPE.PLOT
)
def plot_msdn_cluster_size_histogram(msdn_df: pandas.DataFrame, lang=args.language):
    data = []
    for meta, group in msdn_df.groupby(["s", "cohort", "msdn_id"]):
        data.append(
            [
                *meta,
                group.iloc[0].locus,
                len(group) + 1,
            ]
        )
    fig_df = pandas.DataFrame(
        data, columns=["sample", "cohort", "cluster", "locus", "cluster_size"]
    )
    fig_df.locus = fig_df.locus.astype(str)

    fig = px.histogram(
        fig_df,
        x="cluster_size",
        color="cohort",
        barmode="group",
        log_y=True,
        labels=get_labels(
            lang=lang, labels={"y": get_label("msdn_clusters", lang=lang)}
        ),
    )
    update_legend(fig)
    return fig, fig_df


@StatsFactory.register(
    name="cluster_size_distribution",
    index=1,
    value_col="cluster_size",
    order=matching_order,
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="cluster_size_distribution", output_type=FILE_TYPE.PLOT
)
def plot_msdn_cluster_size_distribution(
    cluster_size_df: pandas.DataFrame, lang=args.language
):
    fig = px.box(
        cluster_size_df,
        y="cluster_size",
        color="cohort",
        hover_data=["sample", "locus", "cluster"],
        labels=get_labels(lang=lang),
    )
    update_legend(fig)
    return fig, cluster_size_df


@StatsFactory.register(name="cluster_size_by_cohort", index=1, order=matching_order)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="cluster_size_by_cohort", output_type=FILE_TYPE.PLOT
)
def plot_msdn_cluster_size_by_cohort(
    msdn_df: pandas.DataFrame,
    age_data: pandas.DataFrame,
    lang=args.language,
):
    data = []
    for meta, group in msdn_df.groupby(["s", "cohort", "msdn_id"]):
        data.append(
            [
                *meta,
                group.iloc[0].locus,
                len(group) + 1,
            ]
        )
    fig_df = pandas.DataFrame(
        data, columns=["sample", "cohort", "cluster", "locus", "cluster_size"]
    )

    data = []
    for meta, group in fig_df.groupby(["cohort", "cluster_size"]):
        divisor = len(msdn_df[msdn_df.cohort == meta[0]].msdn_id.unique())
        data.append([*meta, len(group) / divisor])
    fig_df = pandas.DataFrame(data, columns=["cohort", "cluster_size", "count"])

    fig = px.bar(
        fig_df,
        x="cluster_size",
        y="count",
        color="cohort",
        barmode="group",
        log_y=False,
        labels=get_labels(
            lang=lang, labels={"count": get_label("cluster_size_by_cohort", lang=lang)}
        ),
    )
    update_legend(fig)
    return fig, fig_df


def _get_msdn_pairs_by_exposure_status(msdn_df, age_data, control_matching=None):
    age_data_index = age_data.set_index("rid", drop=False)
    data = []
    for meta, group in msdn_df.groupby(["s", "cohort"]):
        matched_case_id = meta[0]
        if control_matching is not None and meta[0] in control_matching:
            # Retrieve status from matched cases to match control cohort to those cases
            matched_case_id = control_matching[meta[0]]
        exposure_status = (
            age_data_index.loc[matched_case_id].exposure_status
            if matched_case_id in age_data_index.index
            and pandas.notna(age_data_index.loc[matched_case_id].exposure_status)
            else "UNKNOWN"
        )

        data.append(
            [
                *meta,
                exposure_status,
                "{}::{}".format(exposure_status, meta[1]),
                len(group),
            ]
        )
    fig_df = pandas.DataFrame(
        data,
        columns=["sample", "cohort", "exposure_status", "exposure_cohort", "count"],
    )
    fig_df = fill_missing(
        fig_df,
        col_defaults={
            "exposure_status": "NO_EXPOSURE",
            "exposure_cohort": "UNKNOWN",
            "count": 0,
        },
        metadata=age_data,
    )

    fig_df["exposure_status"] = fig_df.apply(
        lambda row: age_data.loc[row["sample"]].exposure_status
        if row["sample"] in age_data.index
        else row["exposure_status"],
        axis=1,
    )
    fig_df = fig_df[fig_df["exposure_status"] != "UNKNOWN"]
    fig_df["exposure_cohort"] = fig_df["cohort"] + "::" + fig_df["exposure_status"]
    return fig_df


@StatsFactory.register(
    name="msdn_pairs_by_exposure_status",
    index=1,
    group_col="exposure_cohort",
    order=matching_order,
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="pairs_by_exposure_status", output_type=FILE_TYPE.PLOT
)
def plot_msdn_pairs_by_exposure_status(
    msdn_df: pandas.DataFrame,
    age_data: pandas.DataFrame,
    control_matching: typing.Optional[typing.Dict[str, str]] = control_matching,
    lang=args.language,
):
    fig_df = _get_msdn_pairs_by_exposure_status(
        msdn_df, age_data, control_matching=control_matching
    )
    fig_df = fig_df[fig_df.exposure_cohort != "CRU::NO_EXPOSURE"]

    fig = px.box(
        fig_df,
        y="count",
        color="exposure_cohort",
        hover_data=["sample"],
        labels=get_labels(
            lang=lang, labels={"count": get_label("msdn_pairs", lang=lang)}
        ),
    )
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="glm_pairs_by_exposure_status",
    output_type=FILE_TYPE.REGRESSION_GLM,
)
def regression_msdn_pairs_by_exposure_status(
    msdn_df, age_data, control_matching=control_matching
):
    fig_df = _get_msdn_pairs_by_exposure_status(
        msdn_df, age_data, control_matching=control_matching
    )
    fig_df = fig_df[fig_df.exposure_cohort != "CRU::NO_EXPOSURE"]
    l = LEVELS_SUBGROUPS
    reg_fit = smf.glm(
        formula="count ~ C(exposure_cohort, levels=l)",
        data=fig_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return fig_df, reg_fit


def _get_msdn_clusters_by_exposure_status(
    msdn_df, age_data, control_matching=control_matching
):
    age_data_index = age_data.set_index("rid", drop=False)
    data = []
    for meta, group in msdn_df.groupby(["s", "cohort"]):
        matched_case_id = meta[0]
        if control_matching is not None and meta[0] in control_matching:
            # Retrieve status from matched cases to match control cohort to those cases
            matched_case_id = control_matching[meta[0]]
        exposure_status = (
            age_data_index.loc[matched_case_id].exposure_status
            if matched_case_id in age_data_index.index
            and pandas.notna(age_data_index.loc[matched_case_id].exposure_status)
            else "UNKNOWN"
        )
        data.append(
            [
                *meta,
                exposure_status,
                "{}::{}".format(exposure_status, meta[1]),
                len(group.msdn_id.unique()),
            ]
        )
    fig_df = pandas.DataFrame(
        data,
        columns=["sample", "cohort", "exposure_status", "exposure_cohort", "count"],
    )
    fig_df = fill_missing(
        fig_df,
        col_defaults={
            "exposure_status": "NO_EXPOSURE",
            "exposure_cohort": "UNKNOWN",
            "count": 0,
        },
        metadata=age_data,
    )
    fig_df["exposure_status"] = fig_df.apply(
        lambda row: age_data.loc[row["sample"]].exposure_status
        if row["sample"] in age_data.index
        else row["exposure_status"],
        axis=1,
    )
    fig_df = fig_df[fig_df["exposure_status"] != "UNKNOWN"]
    fig_df["exposure_cohort"] = fig_df["cohort"] + "::" + fig_df["exposure_status"]
    return fig_df


@StatsFactory.register(
    name="msdn_clusters_by_exposure_status",
    index=1,
    group_col="exposure_cohort",
    order=matching_order,
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="clusters_by_exposure_status", output_type=FILE_TYPE.PLOT
)
def plot_msdn_clusters_by_exposure_status(
    msdn_df: pandas.DataFrame,
    age_data: pandas.DataFrame,
    control_matching: typing.Optional[typing.Dict[str, str]] = control_matching,
    lang=args.language,
):
    fig_df = _get_msdn_clusters_by_exposure_status(
        msdn_df, age_data, control_matching=control_matching
    )
    fig_df = fig_df[fig_df.exposure_cohort != "CRU::NO_EXPOSURE"]

    fig = px.box(
        fig_df,
        y="count",
        color="exposure_cohort",
        hover_data=["sample"],
        labels=get_labels(
            lang=lang, labels={"count": get_label("msdn_clusters", lang=lang)}
        ),
    )
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="glm_clusters_by_exposure_status",
    output_type=FILE_TYPE.REGRESSION_GLM,
)
def regression_msdn_clusters_by_exposure_status(
    msdn_df, age_data, control_matching=control_matching
):
    fig_df = _get_msdn_clusters_by_exposure_status(
        msdn_df, age_data, control_matching=control_matching
    )
    fig_df = fig_df[fig_df.exposure_cohort != "CRU::NO_EXPOSURE"]

    l = LEVELS_SUBGROUPS
    reg_fit = smf.glm(
        formula="count ~ C(exposure_cohort, levels=l)",
        data=fig_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return fig_df, reg_fit


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="msdn_locations", output_type=FILE_TYPE.PLOT
)
def plot_msdn_locations(msdn_df: pandas.DataFrame, lang=args.language):
    msdn_df_coordinate = msdn_df
    msdn_df_coordinate["contig"] = list(map(lambda x: x.contig, msdn_df.locus))
    msdn_df_coordinate["position"] = list(map(lambda x: x.position, msdn_df.locus))
    msdn_df_coordinate = msdn_df_coordinate.sort_values(by=["contig", "position"])
    valid_contigs = set(["X", "Y", *map(str, range(1, 23))])
    msdn_df_coordinate = msdn_df_coordinate[
        msdn_df_coordinate.contig.isin(valid_contigs)
    ]

    fig = px.scatter(
        msdn_df_coordinate,
        x="position",
        y="contig",
        facet_row="cohort",
        color="cohort",
        labels=get_labels(lang=lang),
    )
    fig.update_yaxes(
        categoryorder="array", categoryarray=list([*map(str, range(1, 23)), "X", "Y"])
    )
    update_legend(fig)
    return fig, msdn_df_coordinate


def _get_clusters_by_service_duration(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
):
    age_data["service_duration"] = age_data["service_end"] - age_data["service_begin"]
    age_data_index = age_data.set_index("rid", drop=False)
    data = []
    for meta, group in msdn_df.groupby(["s", "cohort"]):
        service_duration = (
            age_data_index.loc[meta[0]].service_duration
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
                service_duration,
                exposure_status,
                len(group.msdn_id.unique()),
            ]
        )
    fig_df = pandas.DataFrame(
        data,
        columns=["sample", "cohort", "service_duration", "exposure_status", "count"],
    )
    fig_df = fill_missing(
        fig_df,
        col_defaults={
            "exposure_status": "NO_EXPOSURE",
            "service_duration": 0,
            "count": 0,
        },
        metadata=age_data,
    )

    return fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="clusters_by_service_duration", output_type=FILE_TYPE.PLOT
)
def plot_clusters_by_service_duration(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame, lang=args.language
):
    fig_df = _get_clusters_by_service_duration(msdn_df, age_data)
    fig_df = fig_df[fig_df.service_duration > 0]

    _, fit = regression_clusters_by_service_duration(msdn_df, age_data)
    desc = fig_df.describe()

    predict_df = pandas.DataFrame(
        np.linspace(
            desc["service_duration"]["min"], desc["service_duration"]["max"], 100
        ),
        columns=["service_duration"],
    )
    predict_df["count"] = fit.predict(predict_df)

    line = px.line(
        predict_df,
        x="service_duration",
        y="count",
    )

    fig = px.scatter(
        fig_df,
        x="service_duration",
        y="count",
        color="exposure_status",
        symbol="exposure_status",
        # trendline="ols",
        labels=get_labels(
            lang=lang, labels={"count": get_label("msdn_clusters", lang=lang)}
        ),
    )
    fig = merge_figures(fig, line)
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="glm_clusters_by_service_duration",
    output_type=FILE_TYPE.REGRESSION_GLM,
)
def regression_clusters_by_service_duration(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
):
    fig_df = _get_clusters_by_service_duration(msdn_df, age_data)

    reg_fit = smf.glm(
        formula="count ~ service_duration",
        data=fig_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return fig_df, reg_fit


def _get_pairs_by_service_duration(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
):
    age_data["service_duration"] = age_data["service_end"] - age_data["service_begin"]
    age_data_index = age_data.set_index("rid", drop=False)
    data = []
    for meta, group in msdn_df.groupby(["s", "cohort"]):
        service_duration = (
            age_data_index.loc[meta[0]].service_duration
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
                service_duration,
                exposure_status,
                len(group),
            ]
        )
    fig_df = pandas.DataFrame(
        data,
        columns=["sample", "cohort", "service_duration", "exposure_status", "count"],
    )
    fig_df = fill_missing(
        fig_df,
        col_defaults={
            "exposure_status": "NO_EXPOSURE",
            "service_duration": 0,
            "count": 0,
        },
        metadata=age_data,
    )
    return fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="pairs_by_service_duration", output_type=FILE_TYPE.PLOT
)
def plot_pairs_by_service_duration(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame, lang=args.language
):
    fig_df = _get_pairs_by_service_duration(msdn_df, age_data)
    fig_df = fig_df[fig_df.service_duration > 0]

    _, fit = regression_pairs_by_service_duration(msdn_df, age_data)
    desc = fig_df.describe()

    predict_df = pandas.DataFrame(
        np.linspace(
            desc["service_duration"]["min"], desc["service_duration"]["max"], 100
        ),
        columns=["service_duration"],
    )
    predict_df["count"] = fit.predict(predict_df)
    desc = fig_df.describe()

    predict_df = pandas.DataFrame(
        np.linspace(
            desc["service_duration"]["min"], desc["service_duration"]["max"], 100
        ),
        columns=["service_duration"],
    )
    predict_df["count"] = fit.predict(predict_df)

    line = px.line(
        predict_df,
        x="service_duration",
        y="count",
    )

    fig = px.scatter(
        fig_df,
        x="service_duration",
        y="count",
        color="exposure_status",
        symbol="exposure_status",
        # trendline="ols",
        labels=get_labels(
            lang=lang, labels={"count": get_label("msdn_pairs", lang=lang)}
        ),
    )
    fig = merge_figures(fig, line)
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="glm_pairs_by_service_duration",
    output_type=FILE_TYPE.REGRESSION_GLM,
)
def regression_pairs_by_service_duration(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
):
    fig_df = _get_pairs_by_service_duration(msdn_df, age_data)

    reg_fit = smf.glm(
        formula="count ~ service_duration",
        data=fig_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return fig_df, reg_fit


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="mixedlm_cluster_size", output_type=FILE_TYPE.REGRESSION_GLM
)
def regression_cluster_size(cluster_size_df: pandas.DataFrame):
    # smf.mixedlm doesnt work with levels...
    levels = LEVELS_BASE
    formula = "cluster_size ~ C(cohort, levels=levels)"

    y, X = patsy.dmatrices(formula, cluster_size_df, return_type="dataframe")

    mlm = sm.MixedLM(y, X, groups=cluster_size_df["sample"]).fit()

    return cluster_size_df, mlm


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="glm_cluster_size", output_type=FILE_TYPE.REGRESSION_GLM
)
def regression_cluster_size_glm(cluster_size_df: pandas.DataFrame):
    l = LEVELS_BASE
    reg_fit = smf.glm(
        "cluster_size ~ C(cohort, levels=l)",
        data=cluster_size_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return cluster_size_df, reg_fit


def _get_pairs_by_exposition_phase(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
):
    age_data_index = age_data.set_index("rid", drop=False)
    data = []
    for meta, group in msdn_df.groupby(["cohort", "s"]):
        phase = age_data_index.loc[meta[1]].phase if meta[0] == "RADAR" else "INOVA"
        phase = "UNKNOWN" if phase == "UNKOWN" else phase
        data.append(
            [
                *meta,
                phase,
                len(group),
            ]
        )
    fig_df = pandas.DataFrame(data, columns=["cohort", "sample", "phase", "count"])
    fig_df = fill_missing(
        fig_df, col_defaults={"count": 0, "phase": "UNKNOWN"}, metadata=age_data
    )
    fig_df.phase = np.where(fig_df.cohort == "INOVA", "INOVA", fig_df.phase)
    fig_df = fig_df[fig_df.phase != "UNKNOWN"]
    return fig_df


@StatsFactory.register(
    name="pairs_by_exposition_phase",
    index=1,
    group_col="phase",
    order=matching_order,
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="pairs_by_exposition_phase", output_type=FILE_TYPE.PLOT
)
def plot_pairs_by_exposition_phase(msdn_df, age_data, lang=args.language):
    fig_df = _get_pairs_by_exposition_phase(msdn_df, age_data)

    fig = px.box(
        fig_df,
        color="phase",
        y="count",
        labels=get_labels(
            lang=lang, labels={"count": get_label("msdn_pairs", lang=lang)}
        ),
    )
    update_legend(fig)
    return fig, fig_df


def _get_clusters_by_exposition_phase(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
):
    age_data_index = age_data.set_index("rid", drop=False)
    data = []
    for meta, group in msdn_df.groupby(["cohort", "s"]):
        phase = age_data_index.loc[meta[1]].phase if meta[0] == "RADAR" else "INOVA"
        phase = "UNKNOWN" if phase == "UNKOWN" else phase
        data.append(
            [
                *meta,
                phase,
                len(group.msdn_id.unique()),
            ]
        )
    fig_df = pandas.DataFrame(data, columns=["cohort", "sample", "phase", "count"])
    fig_df = fill_missing(
        fig_df, col_defaults={"count": 0, "phase": "UNKNOWN"}, metadata=age_data
    )
    fig_df.phase = np.where(fig_df.cohort == "INOVA", "INOVA", fig_df.phase)
    fig_df = fig_df[fig_df.phase != "UNKNOWN"]
    return fig_df


@StatsFactory.register(
    name="clusters_by_exposition_phase",
    index=1,
    group_col="phase",
    order=matching_order,
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="clusters_by_exposition_phase", output_type=FILE_TYPE.PLOT
)
def plot_clusters_by_exposition_phase(msdn_df, age_data, lang=args.language):
    fig_df = _get_clusters_by_exposition_phase(msdn_df, age_data)

    fig = px.box(
        fig_df,
        color="phase",
        y="count",
        labels=get_labels(
            lang=lang, labels={"count": get_label("msdn_clusters", lang=lang)}
        ),
    )
    update_legend(fig)
    return fig, fig_df


@StatsFactory.register(
    name="cohort_contingency",
    test="chi2",
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="cohort_contingency", output_type=FILE_TYPE.DATAFRAME_TEX
)
def test_cohort_contingency(
    msdn_df: pandas.DataFrame, age_data: pandas.DataFrame
) -> pandas.DataFrame:
    data = []
    for meta, group in msdn_df.groupby(["cohort", "s"]):
        data.append([*meta, len(group)])

    sample_df = pandas.DataFrame(data, columns=["cohort", "sample", "count"])
    sample_df = fill_missing(sample_df, col_defaults={"count": 0}, metadata=age_data)

    contingency = {"has_msdn": {}, "no_msdn": {}}
    for meta, group in sample_df.groupby("cohort"):
        contingency["has_msdn"][meta] = len(group[group["count"] > 0])
        contingency["no_msdn"][meta] = len(group[group["count"] == 0])

    data = [[k, *v.values()] for k, v in contingency.items()]
    return pandas.DataFrame(
        data, columns=["type", *sample_df.cohort.unique()]
    ).set_index("type")


@StatsFactory.register(
    name="msdn_validation",
    index=1,
    order=matching_order,
    value_col="ppv",
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="msdn_validation",
    output_type=FILE_TYPE.PLOT,
)
def plot_msdn_validation(
    msdn_df: pandas.DataFrame,
    lang=args.language,
    expected_evidence=set(["CONFIRMATION", "PARTIAL_CONFIRMATION"]),
):
    data = []
    for meta, group in msdn_df.groupby(["cohort", "s"]):
        data.append(
            [
                *meta,
                len(group[group["graphtyper_evidence"].isin(expected_evidence)])
                / len(group),
            ]
        )

    fig_df = pandas.DataFrame(data, columns=["cohort", "sample", "ppv"])

    fig = px.box(
        fig_df,
        color="cohort",
        y="ppv",
        labels=get_labels(lang=lang),
    )
    update_legend(fig)

    return fig, fig_df


###############################################################################
# Main script                                                                 #
###############################################################################
NamedCheckpoint.enable()

with dump_call_information(
    os.path.join(CACHE_DIR, "args.json"),
    args,
):
    AnalysisFactory.STATS_FACTORY = StatsFactory
    AnalysisFactory.LANGUAGE = args.language
    MsdnAnalysis = AnalysisFactory.instance("MsdnAnalysisFactory")
    MsdnByDosageAnalysis = AnalysisFactory.instance("MsdnByDosageAnalysisFactory")

    MsdnAnalysis.plot_msdn_pairs(msdn_df, age_data)
    MsdnAnalysis.plot_msdn_clusters(msdn_df, age_data)
    MsdnAnalysis.plot_msdn_clusters_violin(msdn_df, age_data)
    _, cluster_size_df = plot_msdn_cluster_size_histogram(msdn_df)
    plot_msdn_cluster_size_distribution(cluster_size_df)
    plot_msdn_cluster_size_by_cohort(msdn_df, age_data)
    plot_msdn_pairs_by_paternal_age(msdn_df, age_data, lang=args.language)
    plot_msdn_clusters_by_paternal_age(msdn_df, age_data, lang=args.language)
    plot_msdn_pairs_by_exposure_status(msdn_df, age_data)
    plot_msdn_clusters_by_exposure_status(msdn_df, age_data)
    MsdnByDosageAnalysis.plot_msdn_pairs_by_dosage(msdn_df, age_data)
    MsdnByDosageAnalysis.plot_msdn_clusters_by_dosage(msdn_df, age_data)
    MsdnByDosageAnalysis.plot_msdn_pairs_by_dosage_and_cohort(msdn_df, age_data)
    MsdnByDosageAnalysis.plot_msdn_clusters_by_dosage_and_cohort(msdn_df, age_data)
    plot_msdn_locations(msdn_df)
    plot_clusters_by_service_duration(msdn_df, age_data)
    plot_pairs_by_service_duration(msdn_df, age_data)
    plot_pairs_by_exposition_phase(msdn_df, age_data)
    plot_clusters_by_exposition_phase(msdn_df, age_data)
    MsdnAnalysis.regression_msdn_pairs(msdn_df, age_data)
    MsdnAnalysis.regression_msdn_clusters(msdn_df, age_data)
    MsdnByDosageAnalysis.regression_msdn_pairs_by_dosage(msdn_df, age_data)
    MsdnByDosageAnalysis.regression_msdn_clusters_by_dosage(msdn_df, age_data)
    MsdnByDosageAnalysis.regression_msdn_pairs_by_dosage_and_cohort(msdn_df, age_data)
    MsdnByDosageAnalysis.regression_msdn_clusters_by_dosage_and_cohort(msdn_df, age_data)
    regression_pairs_by_service_duration(msdn_df, age_data)
    regression_clusters_by_service_duration(msdn_df, age_data)
    regression_cluster_size(cluster_size_df)
    regression_cluster_size_glm(cluster_size_df)
    regression_msdn_pairs_by_exposure_status(msdn_df, age_data)
    regression_msdn_clusters_by_exposure_status(msdn_df, age_data)
    regression_msdn_clusters_by_paternal_age(msdn_df, age_data)
    regression_msdn_pairs_by_paternal_age(msdn_df, age_data)
    test_cohort_contingency(msdn_df, age_data)
    if (
        args.graphtyper is not None
        and args.apply_graphtyper_filter
        and msdn_all_df is not None
    ):
        plot_msdn_validation(msdn_all_df)
    elif args.graphtyper is not None:
        plot_msdn_validation(msdn_df)

    StatsFactory.get().to_dataframe().to_csv(
        os.path.join(CACHE_DIR, "msdn_statistics.csv"),
        index=False,
    )
