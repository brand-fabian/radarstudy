import functools
import logging
import os
import typing
from collections.abc import Sequence

import pandas
from scipy.stats import chi2_contingency, mannwhitneyu, wilcoxon
from statsmodels.stats import weightstats

logger = logging.getLogger(__name__)


class StatisticsRegistrationError(Exception):
    pass


class Statistics:
    def __init__(self):
        self._rows = []

    def chi2_contingency(
        self,
        name: str,
        data: pandas.DataFrame,
    ):
        """Adds a new chi2_contingency test to the object.

        This function can be used to apply the chi2_contingency test
        on the input table. No further sanitizations are applied to
        the input.

        Parameters
        ----------
        name : str
               Name of the test to report in the results.
        data : pandas.DataFrame
               DataFrame containing values to test contingency for. Should
               only contain numerical columns with count values for the
               respective events.
        """

        chi2 = chi2_contingency(data.to_numpy())
        self._rows.append(
            [
                name,
                list(data.columns),
                list(data.index),
                *map(lambda x: None, range(8)),
                "chi2",
                chi2[0],
                chi2[1],
            ]
        )

    def add(
        self,
        name: str,
        data: pandas.DataFrame,
        group_col: str = "cohort",
        value_col: str = "count",
        sample_col: str = "sample",
        groups: typing.Optional[typing.List[str]] = None,
        alternative: str = "two-sided",
        order: typing.List[str] = [],
    ):
        """Add a new row to the statistics object.

        This function computes a ttest and mannwhitneyu-test on the given
        input dataframe. Two groups are expected to be presented using the
        input parameters to compute these statistics.

        Parameters
        ----------
        name : str
               Name of the statistic in the resulting dataframe.
        data : DataFrame
               Input data.
        group_col : str
                    Column name of the column used to group the data.
        value_col : str
                    Column name of the column used to retrieve data values.
        sample_col : str
                     Column name of the column used to retrieve sample names.
        groups : List[str] | None
                 List of column group names to use for aggregating data for
                 statistical tests. If None is set, this defaults to all
                 unique values of the `group_col` column.
        alternative : str
                      Alternative setting for the mannwhitneyu test.
        order : List[str] | None
                List of sample names that will be used to order the input
                data for paired tests.
        """
        if group_col not in data.columns:
            logger.error(
                "missing group column {} in data columns {}".format(
                    group_col, data.columns
                )
            )
            return

        if groups is None:
            groups = list(set(data[group_col]))

        if len(groups) < 2:
            logger.error(
                "can not compute statistics for {} due to missing groups.".format(name)
            )
            return
        elif len(groups) > 2:
            logger.info("got more than two groups in {}.".format(groups))

        if value_col not in data.columns:
            logger.error(
                "missing value column {} in data columns {}".format(
                    value_col,
                    data.columns,
                )
            )
            return

        if order is not None and len(order) > 0 and sample_col in set(data.columns):
            sample_order = {s: i for i, s in enumerate(order)}
            data = pandas.DataFrame(
                map(
                    lambda x: x[1],
                    sorted(
                        data.iterrows(),
                        key=lambda s: sample_order[s[1][sample_col]]
                        if s[1][sample_col] in sample_order
                        else len(order) + 2,
                    ),
                )
            )

        data_groups = list(
            map(
                lambda g: (g, data[data[group_col] == g][value_col]),
                groups,
            )
        )

        descriptives = data.groupby(group_col).describe()
        medians = data.groupby(group_col).median()

        for i, d in enumerate(data_groups):
            for d2 in data_groups[i + 1 :]:
                ttest = weightstats.ttest_ind(d[1], d2[1])
                self._rows.append(
                    [
                        name,
                        d[0],
                        d2[0],
                        descriptives.loc[d[0]][value_col]["mean"],
                        medians.loc[d[0]][value_col],
                        descriptives.loc[d[0]][value_col]["std"],
                        descriptives.loc[d[0]][value_col]["count"],
                        descriptives.loc[d2[0]][value_col]["mean"],
                        medians.loc[d2[0]][value_col],
                        descriptives.loc[d2[0]][value_col]["std"],
                        descriptives.loc[d2[0]][value_col]["count"],
                        "ttest",
                        ttest[0],
                        ttest[1],
                    ]
                )

                utest = mannwhitneyu(d[1], d2[1], alternative=alternative)
                self._rows.append(
                    [
                        name,
                        d[0],
                        d2[0],
                        descriptives.loc[d[0]][value_col]["mean"],
                        medians.loc[d[0]][value_col],
                        descriptives.loc[d[0]][value_col]["std"],
                        descriptives.loc[d[0]][value_col]["count"],
                        descriptives.loc[d2[0]][value_col]["mean"],
                        medians.loc[d2[0]][value_col],
                        descriptives.loc[d2[0]][value_col]["std"],
                        descriptives.loc[d2[0]][value_col]["count"],
                        "utest",
                        utest.statistic,
                        utest.pvalue,
                    ]
                )

                if len(d[1]) == len(d2[1]):
                    # Paired mann-whitney test
                    paired = wilcoxon(d[1], d2[1], alternative=alternative)
                    self._rows.append(
                        [
                            name,
                            d[0],
                            d2[0],
                            descriptives.loc[d[0]][value_col]["mean"],
                            medians.loc[d[0]][value_col],
                            descriptives.loc[d[0]][value_col]["std"],
                            descriptives.loc[d[0]][value_col]["count"],
                            descriptives.loc[d2[0]][value_col]["mean"],
                            medians.loc[d2[0]][value_col],
                            descriptives.loc[d2[0]][value_col]["std"],
                            descriptives.loc[d2[0]][value_col]["count"],
                            "paired-wilcoxon",
                            paired.statistic,
                            paired.pvalue,
                        ]
                    )

    def to_dataframe(
        self,
        column_names: typing.List[str] = [
            "name",
            "group1",
            "group2",
            "group1.mean",
            "group1.median",
            "group1.std",
            "group1.count",
            "group2.mean",
            "group2.median",
            "group2.std",
            "group2.count",
            "test",
            "statistic",
            "p_value",
        ],
    ):
        """Get statistical data as pandas.DataFrame

        Parameters
        ----------
        column_names : List[str]
                       List of size 12 containing the names of the columns in the
                       resulting dataframe.
        """
        return pandas.DataFrame(self._rows, columns=column_names)

    def __len__(self):
        return len(self._rows)


class StatisticsFactory:
    def __init__(self):
        self._statistics: typing.Dict[str, Statistics] = {}

    def register(
        self,
        *,
        at: str = "default",
        name: typing.Optional[str] = None,
        test: str = "default",
        index: int = 0,
        **register_kwargs,
    ):
        """Wrap a function that returns a DataFrame to compute some statistics on it.

        This decorator can be used to apply the statistics class transparently
        on some function return value. \*\*kwargs are pushed as parameters to the
        statistics function used to compute the stats.

        Parameters
        ----------
        at : str
             Name of the statistics instance. Can also be used with get.
        name : str
               Name of the statistic computed. See also description of Statistics.add
        test : str
               Name of the test to conduct. Defaults to ttest/utest for the subgroups.
               Recognized values: ["default", "chi2"]
        index : int
                If the return value is an iterable, we use this index to retrieve the
                dataframe.
        **register_kwargs: Dict
                           Arguments forwarded to the underlying statistics function.


        Returns
        -------
        The original function return values.
        """
        if name is None:
            raise StatisticsRegistrationError("missing name for statistic")

        if at not in self._statistics:
            self._statistics[at] = Statistics()

        def decorator(fn):
            @functools.wraps(fn)
            def wrapper(*args, **kwargs):
                ret_val = fn(*args, **kwargs)
                df = None
                if isinstance(ret_val, Sequence) and len(ret_val) > index:
                    df = ret_val[index]
                elif isinstance(ret_val, pandas.DataFrame):
                    df = ret_val
                else:
                    raise TypeError(
                        "can not use type {} in statistics".format(type(ret_val))
                    )

                if test == "chi2":
                    self._statistics[at].chi2_contingency(name=name, data=df)
                else:
                    self._statistics[at].add(
                        name=name,
                        data=df,
                        **register_kwargs,
                    )

                return ret_val

            return wrapper

        return decorator

    def __getitem__(self, key: str):
        return self._statistics[key]

    def get(self, key: str = "default"):
        try:
            return self.__getitem__(key)
        except KeyError as err:
            logger.error('could not find statistics "{}"'.format(key))
            return Statistics()
