import argparse
import os
import sys
import typing
from itertools import chain

import hail as hl
import pandas
import plotly.express as px

from meta.arguments import (
    add_base_arguments,
    configure_logger,
    dump_call_information,
    hail_init,
    load_metadata_from_file,
)
from meta.checkpoint import NamedCheckpoint
from meta.config import CACHE_DIR, LOG_LEVEL
from meta.loaders import DnmLoader, GraphtyperLoader, MsdnLoader, PhasingLoader
from meta.plots import AnalysisFactory
from meta.plotting import apply_styles, get_label, get_labels, update_legend
from meta.serializers import FILE_TYPE

####################s###########################################################
# Logging                                                                     #
###############################################################################
logger = configure_logger()


###############################################################################
# Arguments                                                                   #
###############################################################################
parser = argparse.ArgumentParser(
    prog="PHASING",
    description="Plot phasing analysis.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "dnm_mt",
    metavar="DNM_MT",
    type=str,
    help="Input matrix table (refined_dnm) containing de novo mutations for both cohorts.",
)
parser.add_argument(
    "msdn_ht",
    metavar="MSDN_HT",
    type=str,
    help="Input hail table (msdn_ht) containing multisite de novo mutations for both cohorts.",
)
parser.add_argument(
    "-p",
    "--phasing",
    type=str,
    default=None,
    required=True,
    help=".pickle file containing phase information for dnms.",
)

parser = add_base_arguments(parser)
args = parser.parse_args()
logger.setLevel(LOG_LEVEL[args.verbose])
logger.debug("Arguments: {}".format(vars(args)))

logger.info("Writing cache output to {}".format(CACHE_DIR))

hail_init(args)

if not os.path.isfile(args.phasing):
    logger.error("Phasing file {} does not exist".format(args.phasing))
    sys.exit(1)

for p in [args.dnm_mt, args.msdn_ht]:
    if not os.path.isdir(p):
        logger.error(
            "Input path {} does not exist or is invalid matrix table".format(p)
        )
        sys.exit(1)

if not args.load_cache:
    NamedCheckpoint.disable()

msdn_ht = MsdnLoader.load_input(args.msdn_ht)
dnm_ht = DnmLoader.load_input(args.dnm_mt)

if args.graphtyper is not None:
    graphtyper_df = GraphtyperLoader.load_graphtyper(
        args.graphtyper,
    )
    graphtyper_ht = GraphtyperLoader.load_graphtyper_hail(graphtyper_df)
    dnm_ht = GraphtyperLoader.annotate_graphtyper_dnm(
        graphtyper_ht,
        dnm_ht,
    )
    msdn_ht = GraphtyperLoader.annotate_graphtyper_msdn(
        graphtyper_ht,
        msdn_ht,
    )

phase_ht = PhasingLoader.load_phasing(args.phasing)

age_data = load_metadata_from_file(args)

logger.info(
    "Read data from {} samples.".format(
        len(age_data),
    )
)

dnm_phase_df = PhasingLoader.annotate_dnm_phase(dnm_ht, phase_ht)
msdn_phase_df = PhasingLoader.annotate_msdn_phase(msdn_ht, phase_ht)

if args.factor > 0:
    logger.info("Applying {}:1 downsampling".format(args.factor))
    (
        dnm_phase_df,
        msdn_phase_df,
        age_data,
        matching,
    ) = PhasingLoader.load_downsample_age_match(
        dnm_phase_df,
        msdn_phase_df,
        age_data,
        args.factor,
    )

    matching_order = [
        *matching.keys(),
        *chain(*map(lambda x: list(x), matching.values())),
    ]
    control_matching = {s: r for r in matching.keys() for s in matching[r]}
    logger.info(
        "Kept {} msdns from {} samples".format(len(msdn_phase_df), len(age_data))
    )
else:
    matching = None
    matching_order = None
    control_matching = None

if not args.apply_control_matching:
    # Overwrite previous control_matching value based on option.
    control_matching = None

if args.apply_graphtyper_filter:
    dnm_phase_df = DnmLoader.load_graphtyper_filter(
        dnm_phase_df,
    )
    msdn_phase_df = MsdnLoader.load_graphtyper_filter(
        msdn_phase_df,
    )

apply_styles()
AnalysisFactory.create_stats_factory()
AnalysisFactory.LANGUAGE = args.language

###############################################################################
# Checkpoints                                                                 #
###############################################################################
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="dnm_phase_absolute", output_type=FILE_TYPE.PLOT
)
def plot_dnm_phase_absolute(dnm_phase_df: pandas.DataFrame, lang=args.language):
    data = []
    for meta, group in dnm_phase_df.groupby(["cohort", "s", "phasing"]):
        data.append(
            [
                *meta,
                len(group),
            ]
        )

    fig_df = pandas.DataFrame(data, columns=["cohort", "sample", "group", "count"])

    fig = px.box(
        fig_df,
        y="count",
        x="group",
        color="cohort",
        labels=get_labels(lang=lang, labels={"count": get_label("dnm", lang=lang)}),
    )
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="dnm_phase_per_sample",
    output_type=FILE_TYPE.PLOT,
)
def plot_dnm_phase_per_sample(dnm_phase_df, lang=args.language):
    data = []
    for meta, group in dnm_phase_df.groupby(["cohort", "s", "phasing"]):
        data.append(
            [
                *meta,
                len(group),
            ]
        )

    fig_df = pandas.DataFrame(data, columns=["cohort", "sample", "group", "count"])

    fig = px.bar(
        fig_df,
        x="sample",
        y="count",
        color="group",
        category_orders={"group": ["UNPHASED", "MATERNAL", "PATERNAL"]},
        facet_row="cohort",
        labels=get_labels(lang=lang, labels={"count": get_label("dnm", lang=lang)}),
    )
    update_legend(fig)
    return fig, fig_df


@AnalysisFactory.STATS_FACTORY.register(
    name="dnm_phased_percentage",
    index=1,
    group_col="cohort_vartype",
    value_col="percentage",
    order=matching_order,
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="dnm_phased_percentage",
    output_type=FILE_TYPE.PLOT,
)
def plot_dnm_phase_percentage(dnm_phase_df: pandas.DataFrame, lang=args.language):
    data = []
    for meta, group in dnm_phase_df.groupby(["cohort", "s"]):
        calls = {
            "phased": len(group[group.is_phased]) / len(group),  # Percentage phased
            "paternal": len(group[group.is_phased & group.is_paternal])
            / len(group[group.is_phased])
            if len(group[group.is_phased]) > 0
            else 0,  # Percentage paternally phased
            "maternal": (
                len(group[group.is_phased])
                - len(group[group.is_phased & group.is_paternal])
            )
            / len(group[group.is_phased])
            if len(group[group.is_phased]) > 0
            else 0,  # Percentage maternally phased
        }

        for k, v in calls.items():
            group_name = "{}::{}".format(meta[0], k.upper())
            data.append([*meta, k, v * 100, group_name])

    fig_df = pandas.DataFrame(
        data, columns=["cohort", "sample", "type", "percentage", "cohort_vartype"]
    )
    fig_df = fig_df[fig_df.percentage > 0]

    fig = px.box(
        fig_df,
        x="type",
        y="percentage",
        color="cohort",
        points="all",
        hover_data=["sample"],
        labels=get_labels(
            lang=lang, labels={"percentage": get_label("percentage_phased", lang=lang)}
        ),
    )
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="msdn_pairs_phased_percentage",
    output_type=FILE_TYPE.PLOT,
)
def plot_msdn_pairs_phase_percentage(
    msdn_phase_df: pandas.DataFrame, lang=args.language
):
    data = []
    for meta, group in msdn_phase_df.groupby(["cohort", "s"]):
        num_phased = len(group[group.phasing != "UNPHASED"])
        calls = {
            "phased": num_phased / len(group),
            "paternal": len(group[group.phasing == "PATERNAL"]) / num_phased
            if num_phased > 0
            else 0,
            "maternal": len(group[group.phasing == "MATERNAL"]) / num_phased
            if num_phased > 0
            else 0,
        }

        for k, v in calls.items():
            data.append([*meta, k, v * 100])

    fig_df = pandas.DataFrame(data, columns=["cohort", "sample", "type", "percentage"])
    fig_df = fig_df[fig_df.percentage > 0]

    fig = px.box(
        fig_df[fig_df.percentage > 0],
        x="cohort",
        y="percentage",
        color="type",
        hover_data=["sample"],
        labels=get_labels(
            lang=lang, labels={"percentage": get_label("percentage_phased", lang=lang)}
        ),
    )
    update_legend(fig)
    return fig, fig_df


def _collapse(vals):
    """Collapse phasing evidence for MSDN clusters."""
    current = "UNPHASED"
    for val in vals:
        if val == "UNPHASED":
            continue
        if val != current and current == "UNPHASED":
            current = val
        elif val != current and current != "UNPHASED":
            current = "INCONCLUSIVE"
    return current


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="msdn_clusters_phased_percentage",
    output_type=FILE_TYPE.PLOT,
)
def plot_msdn_clusters_phase_percentage(
    msdn_phase_df: pandas.DataFrame, lang=args.language
):
    msdn_phase_df["cluster_phasing"] = msdn_phase_df.msdn_id.apply(
        lambda x: _collapse(msdn_phase_df[msdn_phase_df.msdn_id == x].phasing)
    )

    data = []
    for meta, group in msdn_phase_df.groupby(["cohort", "s"]):
        num_phased = len(group[group.cluster_phasing != "UNPHASED"].msdn_id.unique())
        calls = {
            "phased": num_phased / len(group),
            "paternal": len(group[group.cluster_phasing == "PATERNAL"].msdn_id.unique())
            / num_phased
            if num_phased > 0
            else 0,
            "maternal": len(group[group.cluster_phasing == "MATERNAL"].msdn_id.unique())
            / num_phased
            if num_phased > 0
            else 0,
        }

        for k, v in calls.items():
            data.append([*meta, k, v * 100])

    fig_df = pandas.DataFrame(data, columns=["cohort", "sample", "type", "percentage"])
    fig_df = fig_df[fig_df.percentage > 0]

    fig = px.box(
        fig_df[fig_df.percentage > 0],
        x="cohort",
        y="percentage",
        color="type",
        hover_data=["sample"],
        labels=get_labels(
            lang=lang, labels={"percentage": get_label("percentage_phased", lang=lang)}
        ),
    )
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="msdn_clusters_absolute",
    output_type=FILE_TYPE.PLOT,
)
def plot_msdn_clusters_absolute(msdn_phase_df: pandas.DataFrame, lang=args.language):
    msdn_phase_df["cluster_phasing"] = msdn_phase_df.msdn_id.apply(
        lambda x: _collapse(msdn_phase_df[msdn_phase_df.msdn_id == x].phasing)
    )

    data = []
    for meta, group in msdn_phase_df.groupby(["cohort"]):
        # num_phased = len(group[group.cluster_phasing != "UNPHASED"].msdn_id.unique())
        calls = {
            "UNPHASED": len(
                group[group.cluster_phasing == "UNPHASED"].msdn_id.unique()
            ),  #  / len(group),
            "MATERNAL": len(
                group[group.cluster_phasing == "MATERNAL"].msdn_id.unique()
            ),  #  / num_phased if num_phased > 0 else 0,
            "PATERNAL": len(
                group[group.cluster_phasing == "PATERNAL"].msdn_id.unique()
            ),  #  / num_phased if num_phased > 0 else 0,
        }

        for k, v in calls.items():
            data.append([meta, k, v])  # * 100 ])

    fig_df = pandas.DataFrame(data, columns=["cohort", "phasing", "count"])
    fig = px.bar(
        fig_df,
        y="count",
        color="cohort",
        x="phasing",
        barmode="group",
        labels=get_labels(
            lang=lang, labels={"count": get_label("msdn_clusters", lang=lang)}
        ),
    )
    update_legend(fig)
    return fig, fig_df


@AnalysisFactory.STATS_FACTORY.register(
    name="phase_dnm_contingency_radar",
    test="chi2",
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="phase_dnm_contingency_radar",
    output_type=FILE_TYPE.DATAFRAME_TEX,
)
def test_dnm_contingency_radar(dnm_phase_df: pandas.DataFrame):
    return (
        dnm_phase_df.groupby(["cohort", "phasing"])
        .count()["locus"]
        .to_frame()
        .reset_index()
        .pivot(index="phasing", columns="cohort", values="locus")[["RADAR", "INOVA"]]
        .transpose()[["MATERNAL", "PATERNAL"]]
    )


@AnalysisFactory.STATS_FACTORY.register(
    name="phase_dnm_contingency_cru",
    test="chi2",
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="phase_dnm_contingency_cru",
    output_type=FILE_TYPE.DATAFRAME_TEX,
)
def test_dnm_contingency_cru(dnm_phase_df: pandas.DataFrame):
    return (
        dnm_phase_df.groupby(["cohort", "phasing"])
        .count()["locus"]
        .to_frame()
        .reset_index()
        .pivot(index="phasing", columns="cohort", values="locus")[["CRU", "INOVA"]]
        .transpose()[["MATERNAL", "PATERNAL"]]
    )


@AnalysisFactory.STATS_FACTORY.register(
    name="phase_msdn_contingency_radar",
    test="chi2",
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="phase_msdn_contingency_radar",
    output_type=FILE_TYPE.DATAFRAME_TEX,
)
def test_msdn_contingency_radar(msdn_phase_df: pandas.DataFrame):
    return (
        msdn_phase_df.groupby(["cohort", "phasing"])
        .count()["locus"]
        .to_frame()
        .reset_index()
        .pivot(index="phasing", columns="cohort", values="locus")[["RADAR", "INOVA"]]
        .transpose()[["MATERNAL", "PATERNAL"]]
    )


@AnalysisFactory.STATS_FACTORY.register(
    name="phase_msdn_contingency_cru",
    test="chi2",
)
@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="phase_msdn_contingency_cru",
    output_type=FILE_TYPE.DATAFRAME_TEX,
)
def test_msdn_contingency_cru(msdn_phase_df: pandas.DataFrame):
    return (
        msdn_phase_df.groupby(["cohort", "phasing"])
        .count()["locus"]
        .to_frame()
        .reset_index()
        .pivot(index="phasing", columns="cohort", values="locus")[["CRU", "INOVA"]]
        .transpose()[["MATERNAL", "PATERNAL"]]
    )


###############################################################################
# Main script                                                                 #
###############################################################################
NamedCheckpoint.enable()

with dump_call_information(
    os.path.join(CACHE_DIR, "args.json"),
    args,
):
    plot_dnm_phase_absolute(dnm_phase_df)
    plot_dnm_phase_per_sample(dnm_phase_df)
    plot_dnm_phase_percentage(dnm_phase_df)
    plot_msdn_pairs_phase_percentage(msdn_phase_df)
    plot_msdn_clusters_phase_percentage(msdn_phase_df)
    plot_msdn_clusters_absolute(msdn_phase_df)

    test_dnm_contingency_radar(dnm_phase_df)
    test_dnm_contingency_cru(dnm_phase_df)
    test_msdn_contingency_radar(msdn_phase_df)
    test_msdn_contingency_cru(msdn_phase_df)

    stats = AnalysisFactory.STATS_FACTORY.get()
    if len(stats) > 0:
        stats.to_dataframe().to_csv(
            os.path.join(CACHE_DIR, "phasing_statistics.csv"),
            index=False,
        )
