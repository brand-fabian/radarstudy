import argparse
import logging
import os
import pickle
import sys
import typing
from itertools import chain
from multiprocessing import Pool

import hail as hl
import pandas
import plotly.express as px
from tqdm import tqdm

from meta.arguments import (
    add_base_arguments,
    configure_logger,
    dump_call_information,
    hail_init,
    load_metadata_from_file,
)
from meta.checkpoint import NamedCheckpoint
from meta.config import CACHE_DIR, LOG_LEVEL
from meta.loaders import MsdnLoader
from meta.plots import AnalysisFactory
from meta.plotting import apply_styles, get_label, get_labels, update_legend
from meta.serializers import FILE_TYPE

###############################################################################
# Logging                                                                     #
###############################################################################
logger = configure_logger()


###############################################################################
# Arguments                                                                   #
###############################################################################
parser = argparse.ArgumentParser(
    prog="PLOT-PPV",
    description="Simulate the effect of reduced positive predictive values on the statistics.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "ht",
    metavar="HT",
    type=str,
    help="Input hail table (msdn_ht) containing multisite de novo mutations for both cohorts.",
)
parser.add_argument(
    "-n", "--num-simulations", type=int, default=1000, help="Number of simulation runs."
)
parser.add_argument(
    "-p",
    "--ppv",
    type=float,
    default=0.25,
    help="Positive predictive value to simulate",
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
msdn_df = MsdnLoader.load_dataframe(msdn_ht)

age_data = load_metadata_from_file(args)

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

apply_styles()
AnalysisFactory.create_stats_factory()
AnalysisFactory.LANGUAGE = args.language

###############################################################################
# Checkpoints                                                                 #
###############################################################################
def simulation(arguments):
    iteration, msdn_df, age_data = arguments
    SimulationAnalysisFactory = AnalysisFactory.instance(
        "SimulationAnalysisFactory",
        iteration,
        args.ppv,
        matching_order=matching_order,
    )

    return SimulationAnalysisFactory.simulation(
        msdn_df,
        age_data,
    )[1]


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="simulation_statistics_descriptors",
    output_type=FILE_TYPE.PLOT,
)
def plot_simulation_statistics_descriptors(
    simulation_df: pandas.DataFrame, lang=args.language
):
    data = []
    for _, row in simulation_df.iterrows():
        for g in ["group1", "group2"]:
            for t in ["mean", "median", "std"]:
                data.append(
                    [
                        row.name,
                        row[g],
                        t,
                        row["{}.{}".format(g, t)],
                    ]
                )

    fig_df = pandas.DataFrame(data, columns=["name", "cohort", "type", "value"])

    fig = px.box(
        fig_df,
        x="type",
        y="value",
        color="cohort",
        points=False,
        labels=get_labels(lang=lang, labels={"value": get_label("value", lang=lang)}),
    )
    update_legend(fig)
    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR, name="simulation_statistics_base", output_type=FILE_TYPE.PLOT
)
def plot_simulation_statistics_base(
    simulation_df: pandas.DataFrame, lang=args.language
):
    data = []
    for meta, group in simulation_df.groupby(["group1", "group2", "test"]):
        name = "::".join(meta[0:2])
        data.extend(zip([name] * len(group), [meta[2]] * len(group), group["p_value"]))
    fig_df = pandas.DataFrame(data, columns=["comparison", "test", "p_value"])

    fig = px.box(
        fig_df,
        x="comparison",
        y="p_value",
        color="test",
        points="all",
        labels=get_labels(lang=lang),
    )

    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="simulation_statistics_glm_pvalue",
    output_type=FILE_TYPE.PLOT,
)
def plot_simulation_statistics_glm_pvalue(
    simulation_df: pandas.DataFrame,
    lang=args.language,
):
    fig_df = simulation_df.copy()
    fig_df["coef"] = fig_df.coef.map(
        lambda x: "INOVA"
        if x == "Intercept"
        else "RADAR"
        if "RADAR" in x
        else "CRU"
        if "CRU" in x
        else "PILOT"
        if "PILOT" in x
        else None
    )

    fig = px.box(
        fig_df,
        y="p_value",
        color="coef",
        facet_row="group",
        labels=get_labels(lang=lang),
    )

    return fig, fig_df


@NamedCheckpoint.checkpoint(
    prefix=CACHE_DIR,
    name="simulation_statistics_glm_coefficients",
    output_type=FILE_TYPE.PLOT,
)
def plot_simulation_statistics_glm_coefficients(
    simulation_df: pandas.DataFrame,
    lang=args.language,
):
    fig_df = simulation_df.copy()
    fig_df["coef"] = fig_df.coef.map(
        lambda x: "INOVA"
        if x == "Intercept"
        else "RADAR"
        if "RADAR" in x
        else "CRU"
        if "CRU" in x
        else "PILOT"
        if "PILOT" in x
        else None
    )

    fig = px.box(
        fig_df,
        y="coef_value",
        color="coef",
        facet_row="group",
        labels=get_labels(lang=lang),
    )

    return fig, fig_df


###############################################################################
# Main script                                                                 #
###############################################################################
NamedCheckpoint.enable()

with dump_call_information(
    os.path.join(CACHE_DIR, "args.json"),
    args,
):
    simulations = None
    # with Pool(processes=args.threads) as pool:
    #     simulations = pandas.concat(
    #         tqdm(pool.imap(
    #             simulation,
    #             [(i, msdn_df, age_data) for i in range(args.num_simulations)],
    #         ), total=args.num_simulations)
    #     )
    simulations = pandas.concat(
        simulation((i, msdn_df, age_data)) for i in tqdm(range(args.num_simulations))
    )

    if simulations is None:
        logger.error("error simulating ppv data")
        sys.exit(1)

    plot_simulation_statistics_glm_pvalue(simulations)
    plot_simulation_statistics_glm_coefficients(simulations)

    simulations.to_csv(
        os.path.join(CACHE_DIR, "simulations.csv"),
        index=False,
    )

    stats = AnalysisFactory.STATS_FACTORY.get()
    if len(stats) > 0:
        sim_base_df = AnalysisFactory.STATS_FACTORY.get().to_dataframe()

        sim_base_df.to_csv(
            os.path.join(CACHE_DIR, "simulation_statistics.csv"),
            index=False,
        )

        plot_simulation_statistics_descriptors(sim_base_df)
        plot_simulation_statistics_base(sim_base_df)
    else:
        logger.info("no statistics in table for this run")
