import logging
import typing
from collections import namedtuple

import networkx as nx
import numpy as np
import pandas

logger = logging.getLogger(__name__)


###############################################################################
# Helper functions                                                            #
###############################################################################
def fill_missing(
    df: pandas.DataFrame,
    col_defaults: typing.Dict[str, typing.Any],
    metadata: pandas.DataFrame,
    sample_col="sample",
    cohort_col="cohort",
):
    """Add missing rows to the table based on the metadata presented.

    We assume that the cohort and sample column is named the same, both in
    the dataframe and the metadata.

    Rows are added with the columns set to their default values as defined in
    col_defaults for each missing sample present in the metadata but not in
    the actual data.

    Note
    ----
    The col_defaults parameter should include default values for all columns
    except `sample` and `cohort` (as denoted by the respective parameters).
    """

    def _resolve_defaults(
        sample: str,
        cohort: str,
        col_defaults: typing.Dict[str, typing.Any],
        md: pandas.DataFrame = metadata,
    ):
        return [sample, cohort,] + [
            md.loc[sample][key] if key in md.columns and sample in md.index else value
            for key, value in col_defaults.items()
        ]

    def _add_missing(source, cohort, samples, col_defaults) -> pandas.DataFrame:
        return pandas.concat(
            [
                source,
                pandas.DataFrame(
                    map(lambda s: _resolve_defaults(s, cohort, col_defaults), samples),
                    columns=[sample_col, cohort_col, *col_defaults.keys()],
                ),
            ]
        )

    samples = {cohort: set(group.rid) for cohort, group in metadata.groupby(cohort_col)}
    checked_cohorts = set()
    for cohort, group in df.groupby(cohort_col):
        checked_cohorts.add(cohort)
        missing_samples = samples[cohort] - set(group[sample_col])
        df = _add_missing(df, cohort, missing_samples, col_defaults)

    for cohort in set(samples.keys()) - checked_cohorts:
        missing_samples = set(metadata[metadata.cohort == cohort].rid)
        df = _add_missing(df, cohort, missing_samples, col_defaults)

    return df


DownsampleColumnInfo = namedtuple(
    "DownsampleColumnInfo",
    [
        "cohort_col",
        "case_cohort_name",
        "control_cohort_name",
        "sample_col",
        "father_age_col",
        "mother_age_col",
        "reference_sample_col",
    ],
)
default_col_info = DownsampleColumnInfo(
    cohort_col="cohort",
    case_cohort_name="RADAR",
    control_cohort_name="INOVA",
    sample_col="rid",
    father_age_col="father_age",
    mother_age_col="mother_age",
    reference_sample_col="s",
)
ParentalAge = namedtuple("ParentalAge", ["father", "mother"])
default_min_parental_age = ParentalAge(father=15, mother=15)


def downsample_age_match(
    metadata: pandas.DataFrame,
    factor: int = 2,
    reference_data: typing.Optional[pandas.DataFrame] = None,
    column_info: DownsampleColumnInfo = default_col_info,
    divider_char: str = "#",
    min_parental_age: ParentalAge = default_min_parental_age,
    include_other_cohorts: bool = False,
    match_other_cohorts: bool = True,
) -> typing.Tuple[pandas.DataFrame, typing.Dict[str, typing.Set[str]]]:
    """Downsample the metadata dataframe matching the cohorts by age.

    Apply the age matching to a given metadata dataframe. Matching is done by finding
    the minimum weight bipartite match for the two groups. For each sample of the case
    cohort we attempt to find `factor`-many matches in the control cohort.

    If reference data is provided, only control samples with rows in this dataframe
    are considered for the matching.

    To facilitate the factor, we use the character `divider_char` as suffix to the
    actual sample_id to included multiple cases for the case cohort. This character
    should not appear in _any_ sample id.

    Arguments
    ---------
    metadata : pandas.DataFrame
              Pandas DataFrame containing metadata to be sampled.
    factor : int
             Number of control samples to draw per case sample.
    reference_data : pandas.DataFrame
                     DataFrame to check for reference data when deciding which
                     samples to match.
    column_info : DownsampleColumnInfo
                  Column names to use to retrieve the different columns/values from
                  the metadata and reference_data tables.
    divider_char : str
                   Character to split sample id in matching nodes. Will be used
                   as a pattern like `{sample_id}{divider_char}{index}` to derive
                   multiple sample ids for factors > 1.
    min_parental_age: ParentalAge
                      Parental age cutoffs below which the respective age of the
                      parent is corrected/imputed with the age of the partner. A
                      row is removed from the case cohort if both mother and father
                      are below this threshold.
    include_other_cohorts: bool
                           (default = False) If this flag is set, all samples that are
                           part of other cohorts (i.e. not case and control cohorts),
                           as defined in the column_info object will be added to the
                           output table regardless of the matching.
    match_other_cohorts: bool
                         If True, other cohorts will be matched against the case cohort
                         as well.

    Returns
    -------
    pandas.DataFrame
        A dataframe sampled from the metadata table only containing matched samples.
    Dict[str, Set[str]]
        A mapping describing the matching of case to control samples.
    """

    metadata = impute_parental_age(
        metadata, column_info=column_info, min_parental_age=min_parental_age
    )

    cohorts = set(metadata[column_info.cohort_col])
    matched_ids = []
    graph = nx.Graph()

    case_df = metadata[metadata[column_info.cohort_col] == column_info.case_cohort_name]
    case_samples = set(case_df[column_info.sample_col])
    control_cohorts = (
        cohorts - set([column_info.case_cohort_name])
        if match_other_cohorts
        else [column_info.control_cohort_name]
    )
    for control_i, control_cohort in enumerate(control_cohorts):
        control_df = metadata[metadata[column_info.cohort_col] == control_cohort]

        if reference_data is not None:
            control_samples = set(control_df[column_info.sample_col]) & set(
                reference_data[column_info.reference_sample_col]
            )
        else:
            control_samples = set(control_df[column_info.sample_col])

        for i in range(factor):
            graph.add_nodes_from(
                map(
                    lambda s: s
                    + "{d}{i}{d}{ci}".format(d=divider_char, i=i, ci=control_i),
                    case_samples,
                ),
                bipartite=0,
                group=column_info.case_cohort_name,
            )
            if len(case_samples) * factor < len(control_samples):
                matched_ids.extend(
                    map(
                        lambda s: s
                        + "{d}{i}{d}{ci}".format(d=divider_char, i=i, ci=control_i),
                        case_samples,
                    )
                )

        graph.add_nodes_from(
            control_samples,
            bipartite=1,
            group=control_cohort,
        )
        if len(case_samples) * factor >= len(control_samples):
            matched_ids.extend(control_samples)

        edges = []
        for _, case in case_df.iterrows():
            for _, control in control_df.iterrows():
                if control[column_info.sample_col] in control_samples:
                    wght = abs(
                        case[column_info.father_age_col]
                        - control[column_info.father_age_col]
                    ) + abs(
                        case[column_info.mother_age_col]
                        - control[column_info.mother_age_col]
                    )
                    for i in range(factor):
                        edges.append(
                            [
                                case[column_info.sample_col]
                                + "{d}{i}{d}{ci}".format(
                                    d=divider_char, i=i, ci=control_i
                                ),
                                control[column_info.sample_col],
                                {"weight": wght if pandas.notna(wght) else 1000},
                            ]
                        )
        graph.add_edges_from(edges)

    matching = nx.bipartite.minimum_weight_full_matching(graph, top_nodes=matched_ids)

    extract_rid = lambda x: x.split(divider_char)[0]
    samples = set(
        [*map(extract_rid, matching.keys()), *map(extract_rid, matching.values())]
    )

    logger.info(
        "Sampled {} samples for {} control samples".format(
            len(samples), len(case_samples)
        )
    )

    sample_matching = {}
    for x, y in matching.items():
        x = extract_rid(x)
        y = extract_rid(y)
        if y.startswith("R_"):
            # Switch x and y
            z = x
            x = y
            y = z
        if x not in sample_matching:
            sample_matching[x] = set([y])
        else:
            sample_matching[x].add(y)

    downsample_df = metadata[
        metadata[column_info.sample_col].isin(samples)
        | (
            include_other_cohorts
            & ~metadata[column_info.cohort_col].isin(
                [column_info.case_cohort_name, column_info.control_cohort_name]
            )
        )
    ]
    return downsample_df, sample_matching


def impute_parental_age(
    df: pandas.DataFrame,
    column_info: DownsampleColumnInfo = default_col_info,
    min_parental_age: ParentalAge = default_min_parental_age,
) -> pandas.DataFrame:
    """Impute missing paternal or maternal ages.

    If the age of one parent is missing in the dataframe (as per the column
    info), or reported incorrectly, we impute it with the age of the partner
    here.

    Rows where both ages are reported incorrectly or missing will be removed
    from the resulting table.

    Parameters
    ----------
    df : pandas.DataFrame
         Pandas DataFrame where the age columns will be set.
    column_info : DownsampleColumnInfo
                  Column names that will be used to set or retrieve all the
                  different values.
    min_parental_age : ParentalAge
                       Minimum parental age cutoff below which values will
                       be treated as missing.

    Returns
    -------
    pandas.DataFrame
        Copy of `df` with the parental age columns adjusted according to the
        imputing guidelines.
    """
    copy_df = df.copy(deep=True)
    copy_df = copy_df[
        ~(
            (copy_df[column_info.father_age_col] < min_parental_age.father)
            & (copy_df[column_info.mother_age_col] < min_parental_age.mother)
        )
    ]
    copy_df[column_info.father_age_col] = np.where(
        copy_df[column_info.father_age_col] < min_parental_age.father,
        copy_df[column_info.mother_age_col],
        copy_df[column_info.father_age_col],
    )
    copy_df[column_info.mother_age_col] = np.where(
        copy_df[column_info.mother_age_col] < min_parental_age.mother,
        copy_df[column_info.father_age_col],
        copy_df[column_info.mother_age_col],
    )
    return copy_df


def filter_sites(
    df: pandas.DataFrame,
    valid_sites: pandas.DataFrame,
    locus_cols: typing.Iterable[str] = ["locus", "previous.locus"],
    valid_col: str = "downsample_cov",
    valid_locus_col: str = "locus",
    sample_col: str = "s",
    cohort_col: str = "cohort",
) -> pandas.DataFrame:
    is_valid = valid_sites[valid_sites[valid_col]].set_index(
        [
            cohort_col,
            sample_col,
            valid_locus_col,
        ],
        drop=False,
    )

    copy_df = df.copy(deep=True)
    for locus_col in locus_cols:
        copy_df = copy_df.set_index([cohort_col, sample_col, locus_col], drop=False)
        copy_df = copy_df[copy_df.index.isin(is_valid.index)]

    copy_df.reset_index(drop=True, inplace=True)
    return copy_df
