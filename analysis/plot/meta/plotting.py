import json
import logging
import os
import typing
from collections import namedtuple

import plotly.express as px

import meta.config as config

PlotlyStyles = namedtuple(
    "PlotlyStyles",
    [
        "template",
        "color_continuous_scale",
        "color_discrete_sequence",
        "width",
        "height",
        "symbol_sequence",
        "category_orders",
    ],
)
DEFAULT_PLOTLY_STYLES = PlotlyStyles(
    template="ggplot2",
    color_continuous_scale=px.colors.sequential.Blackbody,
    color_discrete_sequence=px.colors.qualitative.G10,
    width=1000,
    height=1000,
    symbol_sequence=[0, 4, 2, 3, 1, 5],
    category_orders={
        "cohort": config.CATEGORY_ORDER_BASE,
        "exposure_cohort": config.CATEGORY_ORDER_SUBGROUPS,
    },
)


def apply_styles(styles: PlotlyStyles = DEFAULT_PLOTLY_STYLES) -> None:
    """Apply plotly express styling."""
    px.defaults.template = styles.template
    px.defaults.color_continuous_scale = styles.color_continuous_scale
    px.defaults.color_discrete_sequence = styles.color_discrete_sequence
    px.defaults.width = styles.width
    px.defaults.height = styles.height
    px.defaults.symbol_sequence = styles.symbol_sequence
    px.defaults.category_orders = styles.category_orders


def merge_figures(fig, other):
    """Combine traces of two plotly figures.

    To merge two figures, we copy the traces from the figure `other` to the
    figure `fig`.

    Parameters
    ----------
    fig, other : plotly.Figure
                 The plotly figures that will be merged.

    Returns
    -------
    plotly.Figure
        The first figure, with traces from the second figure added.
    """
    for trace in other.data:
        fig.add_trace(trace)

    return fig


def update_legend(fig):
    """Style the figure legend."""
    fig.update_layout(
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "xanchor": "right",
            "x": 1,
            "y": 1.02,
        }
    )
    return fig


base_labels = {}
with open(
    os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "labels.json",
    )
) as json_f:
    base_labels = json.load(json_f)


def get_labels(
    lang: str = "de",
    labels: typing.Dict[str, str] = {},
    base_labels: typing.Dict[str, typing.Dict[str, str]] = base_labels,
) -> typing.Dict[str, str]:
    """Get plot axis labels.

    This function gets plot axis labels in the given language.
    Either the base label is returned (from the json file in adjacent
    to this file) or the label passed to the function. The label passed
    to this function takes precedence.

    Parameters
    ----------
    lang : str
           Language id in base_labels to use.
    labels : Dict[str, str]
             Dictionary of axis labels.
    base_labels : Dict[str, Dict[str, str]]
                  Dictionary of axis labels in different languages.

    Returns
    -------
    Dict[str, str]
        Axis labels derived from the language and labels passed in.
    """

    base = base_labels[lang] if lang in base_labels else {}
    if len(labels) == 0:
        return base
    else:
        return {**base, **labels}


def get_label(label: str, lang: str = "de") -> str:
    """Get a label for a given dimension in a given language.

    Parameters
    ----------
    label : str
            Label to retrieve.
    lang : str
           Language to retrieve label for.

    Returns
    -------
    str
        Label, returns the argument `label` if no match was found.
    """
    return (
        base_labels[lang][label]
        if lang in base_labels and label in base_labels[lang]
        else label
    )
