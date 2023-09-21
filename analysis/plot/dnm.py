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
from meta.config import CACHE_DIR, LOG_LEVEL, LEVELS_BASE
from meta.helpers import impute_parental_age
from meta.loaders import DnmLoader, GraphtyperLoader, MsdnLoader
from meta.plots import AnalysisFactory
from meta.plotting import apply_styles, get_labels, merge_figures, update_legend
from meta.serializers import FILE_TYPE
from meta.statistics import StatisticsFactory

# Statistics Factory used to get the main output statistics
StatsFactory = StatisticsFactory()

###############################################################################
# Logging                                                                     #
###############################################################################
logger = configure_logger()

###############################################################################
# Arguments and script preamble                                               #
###############################################################################
parser = argparse.ArgumentParser(
    prog="PLOT-DNM",
    description="Plot descriptive statistics for dnms.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "mt",
    metavar="MT",
    type=str,
    help="Input matrix table (refined_dnm) containing de novo mutations for both cohorts.",
)
parser.add_argument(
    "--msdn-ht",
    type=str,
    default=None,
    help="Input hail table (msdn_ht) containing multisite de novo mutations for both cohorts. Required for --isolated",
)
parser.add_argument(
    "--isolated",
    default=False,
    action="store_true",
    help="If true, compute isolated dnm rates (i.e. exlude msdn variants).",
)
parser = add_base_arguments(parser)

args = parser.parse_args()
logger.setLevel(LOG_LEVEL[args.verbose])
logger.debug("Arguments: {}".format(vars(args)))

logger.info("Writing cache output to {}".format(CACHE_DIR))

hail_init(args)

if args.isolated and args.msdn_ht is None:
    logger.error("--msdn-ht is required for computing isolated dnms.")
    sys.exit(1)
if args.isolated and not os.path.isdir(args.msdn_ht):
    logger.error("could not find msdn hail table at {}".format(args.msdn_ht))
    sys.exit(1)

if not os.path.isdir(args.mt):
    logger.error(
        "Input path {} does not exist or is invalid matrix table".format(args.mt)
    )
    sys.exit(1)

if not args.load_cache:
    NamedCheckpoint.disable()

dnm_ht = DnmLoader.load_input(args.mt)

if args.graphtyper is not None:
    logger.info("loading graphtyper results data")
    graphtyper_df = GraphtyperLoader.load_graphtyper(
        args.graphtyper,
    )
    graphtyper_ht = GraphtyperLoader.load_graphtyper_hail(graphtyper_df)
    dnm_ht = GraphtyperLoader.annotate_graphtyper_dnm(
        graphtyper_ht,
        dnm_ht,
    )

dnm_df = DnmLoader.load_dataframe(dnm_ht, has_graphtyper=args.graphtyper is not None)

if args.isolated:
    logger.info("computing isolated de novo mutations")
    msdn_ht = MsdnLoader.load_input(args.msdn_ht)
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
    dnm_df = DnmLoader.load_isolated(dnm_df, msdn_df)

age_data = load_metadata_from_file(args)

logger.info(
    "Read {} dnms from {} samples.".format(
        len(dnm_df),
        len(age_data),
    )
)

if args.factor > 0:
    ##
    # Matching of case to control cohort. Besides eliminating confounders (i.e. age)
    # from the dataset, we also use this to inform some stratification in the statistics
    # layer to properly compare the two cohorts (by matching sub-cohorts for exposed
    # and non-exposed individuals) in subsequent statitsic tests.
    #
    # This can be overriden by the apply_control_matching flag.
    logger.info("Applying {}:1 downsampling".format(args.factor))
    dnm_df, age_data, matching = DnmLoader.load_downsample_age_match(
        dnm_df, age_data, args.factor
    )
    matching_order = [
        *matching.keys(),
        *chain(*map(lambda x: list(x), matching.values())),
    ]
    control_matching = {s: r for r in matching.keys() for s in matching[r]}
    logger.info(
        "Kept {} dnms from {} samples".format(
            len(dnm_df),
            len(age_data),
        )
    )
else:
    matching = None
    matching_order = None
    control_matching = None
    age_data = impute_parental_age(age_data)

if not args.apply_control_matching:
    # Overwrite previous control_matching value based on option.
    control_matching = None

if args.apply_graphtyper_filter:
    logger.info("discarding all variants not found by graphtyper")
    dnm_all_df = dnm_df.copy(deep=True)
    dnm_df = DnmLoader.load_graphtyper_filter(dnm_df)

apply_styles()

###############################################################################
# Checkpoints                                                                 #
###############################################################################
# These scriptlets/functions generate named output files in the cache dir that
# feature plots and datatables used to generate the data.
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="substitutions", output_type=FILE_TYPE.PLOT
)
def plot_substitutions(dnm_df, lang=args.language):
    subst_data = []
    for subst_type, group in dnm_df.groupby(["str", "s", "cohort"]):
        subst_data.append([*subst_type, len(group)])
    subst_df = pandas.DataFrame(
        subst_data, columns=["substitution", "sample", "cohort", "count"]
    )

    fig = px.box(
        subst_df,
        y="count",
        points="outliers",
        x="substitution",
        color="cohort",
        hover_data=["sample", "count", "cohort"],
        labels=get_labels(lang=lang, labels={"count": "DNMs"})
        # boxmean=True,
    )

    update_legend(fig)
    return fig, subst_df


@StatsFactory.register(name="dnms", index=1, order=matching)
@NamedCheckpoint.checkpoint(prefix=CACHE_DIR, name="dnms", output_type=FILE_TYPE.PLOT)
def plot_dnms(dnm_df, age_data, lang=args.language):
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

    fig = px.box(
        dnm_per_sample_df,
        y="count",
        color="cohort",
        hover_data=["sample"],
        labels=get_labels(lang=lang, labels={"count": "DNMs"}),
    )

    update_legend(fig)
    return fig, dnm_per_sample_df


@StatsFactory.register(
    name="dnms_by_exposure", index=1, group_col="exposure_cohort", order=matching_order
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="dnms_by_exposure", output_type=FILE_TYPE.PLOT
)
def plot_dnms_by_exposure(
    dnm_df, age_data, control_matching=control_matching, lang=args.language
):
    age_data_index = age_data.set_index("rid", drop=False)
    dnm_data = []
    for sample, group in dnm_df.groupby(["s", "cohort"]):
        matched_case_id = sample[0]
        if control_matching is not None and sample[0] in control_matching:
            # Retrieve status from matched cases to match control cohort to those cases
            matched_case_id = control_matching[sample[0]]
        exposure_status = (
            age_data_index.loc[matched_case_id].exposure_status
            if matched_case_id in age_data_index.index
            and pandas.notna(age_data_index.loc[matched_case_id].exposure_status)
            else "UNKNOWN"
        )

        dnm_data.append(
            [
                *sample,
                exposure_status,
                len(group),
            ]
        )
    dnm_per_sample_df = pandas.DataFrame(
        dnm_data, columns=["sample", "cohort", "exposure_status", "count"]
    )
    dnm_per_sample_df = dnm_per_sample_df[
        dnm_per_sample_df["exposure_status"] != "UNKNOWN"
    ]
    dnm_per_sample_df["exposure_cohort"] = (
        dnm_per_sample_df["cohort"] + "::" + dnm_per_sample_df["exposure_status"]
    )

    fig = px.box(
        dnm_per_sample_df,
        y="count",
        color="exposure_cohort",
        hover_data=["sample"],
        labels=get_labels(lang=lang, labels={"count": "DNMs"}),
    )

    update_legend(fig)
    return fig, dnm_per_sample_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="dnms_by_dosage", output_type=FILE_TYPE.PLOT
)
def plot_dnms_by_dosage(dnm_df, age_data, lang=args.language):
    age_data_index = age_data.set_index("rid", drop=False)
    dnm_data = []
    for sample, group in dnm_df.groupby(["s", "cohort"]):
        dnm_data.append(
            [
                *sample,
                age_data_index.loc[sample[0]].sum_dosage_gonads
                if sample[0] in age_data_index.index
                else 0,
                len(group),
            ]
        )
    exposure_gonads_df = pandas.DataFrame(
        dnm_data,
        columns=["sample", "cohort", "dosage", "count"],
    )
    exposure_gonads_df = exposure_gonads_df[exposure_gonads_df.dosage > 0]

    fig = px.scatter(
        exposure_gonads_df,
        x="dosage",
        y="count",
        color="cohort",
        symbol="cohort",
        labels=get_labels(lang=lang, labels={"count": "DNMs"}),
    )
    update_legend(fig)
    return fig, exposure_gonads_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="paternal_age_effect", output_type=FILE_TYPE.PLOT
)
def plot_paternal_age_effect(dnm_df, age_data, lang=args.language):
    age_data_index = age_data.set_index("rid", drop=False)
    dnm_data = []
    for sample, group in dnm_df.groupby(["s", "cohort"]):
        dnm_data.append(
            [
                *sample,
                age_data_index.loc[sample[0]].father_age
                if sample[0] in age_data_index.index
                else 0,
                len(group),
            ]
        )
    father_age_df = pandas.DataFrame(
        dnm_data, columns=["sample", "cohort", "father_age", "count"]
    )

    # Compute trend line based on statistical model that will be used
    # later.
    _, fit = regression_paternal_age_effect(father_age_df)
    desc = father_age_df.groupby("cohort").describe()

    ages = {
        "radar": np.linspace(
            desc["father_age"]["min"].loc["RADAR"],
            desc["father_age"]["max"].loc["RADAR"],
            100,
        ),
        "inova": np.linspace(
            desc["father_age"]["min"].loc["INOVA"],
            desc["father_age"]["max"].loc["INOVA"],
            100,
        ),
        "cru": np.linspace(
            desc["father_age"]["min"].loc["CRU"],
            desc["father_age"]["max"].loc["CRU"],
            100,
        ),
    }

    predict_df = pandas.concat(
        [
            pandas.DataFrame(
                map(lambda x: ["INOVA", x], ages["inova"]),
                columns=["cohort", "father_age"],
            ),
            pandas.DataFrame(
                map(lambda x: ["RADAR", x], ages["radar"]),
                columns=["cohort", "father_age"],
            ),
            pandas.DataFrame(
                map(lambda x: ["CRU", x], ages["cru"]),
                columns=["cohort", "father_age"],
            ),
        ]
    )

    predict_df["count"] = fit.predict(predict_df)

    line = px.line(
        predict_df,
        x="father_age",
        y="count",
        color="cohort",
        category_orders={"cohort": ["INOVA", "RADAR", "CRU"]},
    )
    line.data[0].name = "GLM Inova"
    line.data[1].name = "GLM Radar"
    line.data[2].name = "GLM CRU"

    fig = px.scatter(
        father_age_df,
        x="father_age",
        y="count",
        color="cohort",
        symbol="cohort",
        hover_data=["sample", "count", "father_age"],
        labels=get_labels(lang=lang, labels={"count": "DNMs"}),
    )
    fig = merge_figures(fig, line)
    update_legend(fig)
    return fig, father_age_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="glm_paternal_age_effect",
    output_type=FILE_TYPE.REGRESSION_GLM,
)
def regression_paternal_age_effect(father_age_df: pandas.DataFrame):
    l = LEVELS_BASE
    reg_fit = smf.glm(
        formula="count ~ C(cohort, levels=l) + father_age + C(cohort, levels=l) * father_age",
        data=father_age_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return father_age_df, reg_fit


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="maternal_age_effect", output_type=FILE_TYPE.PLOT
)
def plot_maternal_age_effect(dnm_df, age_data, lang=args.language):
    age_data_index = age_data.set_index("rid", drop=False)
    dnm_data = []
    for sample, group in dnm_df.groupby(["s", "cohort"]):
        dnm_data.append(
            [
                *sample,
                age_data_index.loc[sample[0]].mother_age
                if sample[0] in age_data_index.index
                else 0,
                len(group),
            ]
        )
    mother_age_df = pandas.DataFrame(
        dnm_data, columns=["sample", "cohort", "mother_age", "count"]
    )
    mother_age_df = mother_age_df[mother_age_df.mother_age > 15]

    fig = px.scatter(
        mother_age_df,
        x="mother_age",
        y="count",
        color="cohort",
        symbol="cohort",
        hover_data=["sample", "count", "mother_age"],
        trendline="ols",
        labels=get_labels(lang=lang, labels={"count": "DNMs"}),
    )
    update_legend(fig)
    return fig, mother_age_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="alternate_allele_frequencies", output_type=FILE_TYPE.PLOT
)
def plot_alternate_allele_frequencies(dnm_df, lang=args.language):
    aaf_hist_data = []
    aaf_grouped_df = dnm_df
    aaf_grouped_df["af_groups"] = pandas.cut(
        aaf_grouped_df["AF"],
        bins=list(x / 100 for x in range(30, 70, 1)) + [1.1],
        labels=list(x / 100 for x in range(30, 70, 1)),
    )
    for desc, group in aaf_grouped_df.groupby(["cohort", "s", "af_groups"]):
        aaf_hist_data.append(
            [
                *desc,
                len(group),
            ]
        )
    aaf_hist_df = pandas.DataFrame(
        aaf_hist_data,
        columns=[
            "cohort",
            "sample",
            "AF",
            "count",
        ],
    )

    fig = px.line(
        aaf_hist_df,
        x="AF",
        y="count",
        hover_data=["sample", "cohort"],
        color="sample",
        labels=get_labels(lang=lang, labels={"count": "DNMs"}),
    )
    update_legend(fig)
    return fig, aaf_hist_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="dnms_by_service_duration", output_type=FILE_TYPE.PLOT
)
def plot_dnms_by_service_duration(dnm_df, age_data, lang=args.language):
    age_data["service_duration"] = age_data["service_end"] - age_data["service_begin"]
    age_data_index = age_data.set_index("rid", drop=False)
    fig_data = []
    for sample, group in dnm_df.groupby(["s", "cohort"]):
        fig_data.append(
            [
                *sample,
                age_data_index.loc[sample[0]].service_duration
                if sample[0] in age_data_index.index
                else 0,
                len(group),
            ]
        )
    fig_df = pandas.DataFrame(
        fig_data, columns=["sample", "cohort", "service_duration", "count"]
    )
    fig_df = fig_df[fig_df.service_duration > 0]

    fig = px.scatter(
        fig_df,
        x="service_duration",
        y="count",
        color="cohort",
        symbol="cohort",
        hover_data=["sample", "count", "service_duration"],
        trendline="ols",
        labels=get_labels(lang=lang, labels={"count": "DNMs"}),
    )
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="glm_dnms_by_service_duration",
    output_type=FILE_TYPE.REGRESSION_GLM,
)
def regression_dnm_by_service_duration(
    dnm_df: pandas.DataFrame, age_data: pandas.DataFrame
):
    age_data["service_duration"] = age_data["service_end"] - age_data["service_begin"]
    age_data_index = age_data.set_index("rid", drop=False)
    fig_data = []
    for sample, group in dnm_df.groupby(["s", "cohort"]):
        fig_data.append(
            [
                *sample,
                age_data_index.loc[sample[0]].service_duration
                if sample[0] in age_data_index.index
                else 0,
                len(group),
            ]
        )
    fig_df = pandas.DataFrame(
        fig_data, columns=["sample", "cohort", "service_duration", "count"]
    )

    fit = smf.glm(
        "count ~ service_duration", data=fig_df, family=sm.families.NegativeBinomial()
    ).fit()

    return fig_df, fit


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="glm_dnms_by_dosage", output_type=FILE_TYPE.REGRESSION_GLM
)
def regression_dnms_by_dosage(dnm_df: pandas.DataFrame, age_data: pandas.DataFrame):
    age_data_index = age_data.set_index("rid", drop=False)
    dnm_data = []
    for sample, group in dnm_df.groupby(["s", "cohort"]):
        dnm_data.append(
            [
                *sample,
                age_data_index.loc[sample[0]].sum_dosage_gonads
                if sample[0] in age_data_index.index
                else 0,
                len(group),
            ]
        )
    exposure_gonads_df = pandas.DataFrame(
        dnm_data,
        columns=["sample", "cohort", "dosage", "count"],
    )

    l = LEVELS_BASE
    fit = smf.glm(
        "count ~ dosage + C(cohort, levels=l)",
        data=exposure_gonads_df,
        family=sm.families.NegativeBinomial(),
    ).fit()

    return exposure_gonads_df, fit


@StatsFactory.register(name="dnms_by_phase", index=1, group_col="phase", order=matching)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="dnms_by_phase", output_type=FILE_TYPE.PLOT
)
def plot_dnms_by_phase(
    dnm_df: pandas.DataFrame, age_data: pandas.DataFrame, lang=args.language
):
    age_data_index = age_data.set_index("rid", drop=False)
    data = []
    for meta, group in dnm_df.groupby(["cohort", "s"]):
        data.append(
            [
                *meta,
                age_data_index.loc[meta[1]].phase if meta[0] == "RADAR" else "INOVA",
                len(group),
            ]
        )

    data = filter(lambda x: x[2] != "UNKNOWN", data)
    fig_df = pandas.DataFrame(data, columns=["cohort", "sample", "phase", "count"])

    fig = px.box(
        fig_df,
        color="phase",
        y="count",
        labels=get_labels(lang=lang, labels={"count": "DNMs"}),
    )
    update_legend(fig)

    return fig, fig_df


@StatsFactory.register(
    name="dnm_validation",
    index=1,
    order=matching_order,
    value_col="ppv",
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="dnm_validation",
    output_type=FILE_TYPE.PLOT,
)
def plot_dnm_validation(
    dnm_df: pandas.DataFrame,
    lang=args.language,
):
    data = []
    for meta, group in dnm_df.groupby(["cohort", "s"]):
        data.append(
            [
                *meta,
                len(group[group["graphtyper_evidence"] == "CONFIRMATION"]) / len(group),
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
    DnmAnalysis = AnalysisFactory.instance("DnmAnalysisFactory", matching_order=matching)
    MetadataAnalysis = AnalysisFactory.instance("MetadataAnalysisFactory", matching_order=matching)

    plot_substitutions(dnm_df)
    DnmAnalysis.plot_dnms(dnm_df, age_data)
    DnmAnalysis.regression_dnms(dnm_df, age_data)
    plot_dnms_by_exposure(dnm_df, age_data)
    plot_dnms_by_dosage(dnm_df, age_data)
    _, father_age_df = plot_paternal_age_effect(dnm_df, age_data)
    regression_paternal_age_effect(father_age_df)
    plot_maternal_age_effect(dnm_df, age_data)
    MetadataAnalysis.plot_age_distribution(age_data)
    plot_alternate_allele_frequencies(dnm_df)
    MetadataAnalysis.plot_service_duration(age_data)
    MetadataAnalysis.plot_dosage_by_service_duration(age_data)
    plot_dnms_by_service_duration(dnm_df, age_data)
    plot_dnms_by_phase(dnm_df, age_data)
    regression_dnm_by_service_duration(dnm_df, age_data)
    regression_dnms_by_dosage(dnm_df, age_data)
    MetadataAnalysis.plot_dosage_histogram(age_data)
    MetadataAnalysis.plot_parental_age_histogram(age_data)
    if (
        args.graphtyper is not None
        and args.apply_graphtyper_filter
        and dnm_all_df is not None
    ):
        plot_dnm_validation(dnm_all_df)
    elif args.graphtyper is not None:
        plot_dnm_validation(dnm_df)

    StatsFactory.get().to_dataframe().to_csv(
        os.path.join(CACHE_DIR, "dnm_statistics.csv"),
        index=False,
    )
