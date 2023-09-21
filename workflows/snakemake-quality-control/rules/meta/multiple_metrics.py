import typing
from enum import Enum


class METRICS_SOURCE(Enum):
    PICARD = 'PICARD'
    PARABRICKS = 'PARABRICKS'


class METRICS_OUTPUT_TYPE(Enum):
    DATA = 'txt'
    REPORT = 'pdf'
    IMAGE = 'png'


class METRIC_TOOL(Enum):
    COLLECT_ALIGNMENT_SUMMARY_METRICS = 'CollectAlignmentSummaryMetrics'
    COLLECT_INSERT_SIZE_METRICS = 'CollectInsertSizeMetrics'
    QUALITY_SCORE_DISTRIBUTION = 'QualityScoreDistribution'
    MEAN_QUALITY_BY_CYCLE = 'MeanQualityByCycle'
    COLLECT_BASE_DISTRIBUTION_BY_CYCLE = 'CollectBaseDistributionByCycle'
    COLLECT_GC_BIAS_METRICS = 'CollectGcBiasMetrics'
    COLLECT_SEQUENCING_ARTIFACT_METRICS = 'CollectSequencingArtifactMetrics'
    COLLECT_QUALITY_YIELD_METRICS = 'CollectQualityYieldMetrics'


METRICS_FILES: typing.Dict[
    METRIC_TOOL,
    typing.Dict[
        METRICS_SOURCE,
        typing.Dict[
            METRICS_OUTPUT_TYPE,
            typing.Union[str, typing.List[str], None]
        ]
    ]
] = {
    METRIC_TOOL.COLLECT_ALIGNMENT_SUMMARY_METRICS: {
        METRICS_SOURCE.PICARD: {
            METRICS_OUTPUT_TYPE.DATA: "alignment_summary_metrics",
            METRICS_OUTPUT_TYPE.REPORT: "read_length_histogram.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: None,
        },
        METRICS_SOURCE.PARABRICKS: {
            METRICS_OUTPUT_TYPE.DATA: "alignment.txt",
            METRICS_OUTPUT_TYPE.REPORT: "read_length_histogram.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: None,
        }
    },
    METRIC_TOOL.COLLECT_INSERT_SIZE_METRICS: {
        METRICS_SOURCE.PICARD: {
            METRICS_OUTPUT_TYPE.DATA: "insert_size_metrics",
            METRICS_OUTPUT_TYPE.REPORT: "insert_size_histogram.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: None,
        },
        METRICS_SOURCE.PARABRICKS: {
            METRICS_OUTPUT_TYPE.DATA: "insert_size.txt",
            METRICS_OUTPUT_TYPE.REPORT: "insert_size.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: "insert_size.png",
        }
    },
    METRIC_TOOL.QUALITY_SCORE_DISTRIBUTION: {
        METRICS_SOURCE.PICARD: {
            METRICS_OUTPUT_TYPE.DATA: "quality_distribution_metrics",
            METRICS_OUTPUT_TYPE.REPORT: "quality_distribution.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: None,
        },
        METRICS_SOURCE.PARABRICKS: {
            METRICS_OUTPUT_TYPE.DATA: "qualityscore.txt",
            METRICS_OUTPUT_TYPE.REPORT: "qualityscore.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: "qualityscore.png",
        }
    },
    METRIC_TOOL.MEAN_QUALITY_BY_CYCLE: {
        METRICS_SOURCE.PICARD: {
            METRICS_OUTPUT_TYPE.DATA: "quality_by_cycle_metrics",
            METRICS_OUTPUT_TYPE.REPORT: "quality_by_cycle.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: None,
        },
        METRICS_SOURCE.PARABRICKS: {
            METRICS_OUTPUT_TYPE.DATA: "mean_quality_by_cycle.txt",
            METRICS_OUTPUT_TYPE.REPORT: "mean_quality_by_cycle.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: "mean_quality_by_cycle.png",
        }
    },
    METRIC_TOOL.COLLECT_BASE_DISTRIBUTION_BY_CYCLE: {
        METRICS_SOURCE.PICARD: {
            METRICS_OUTPUT_TYPE.DATA: "base_distribution_by_cycle_metrics",
            METRICS_OUTPUT_TYPE.REPORT: "base_distribution_by_cycle.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: None,
        },
        METRICS_SOURCE.PARABRICKS: {
            METRICS_OUTPUT_TYPE.DATA: "base_distribution_by_cycle.txt",
            METRICS_OUTPUT_TYPE.REPORT: "base_distribution_by_cycle.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: "base_distribution_by_cycle.png",
        }
    },
    METRIC_TOOL.COLLECT_GC_BIAS_METRICS: {
        METRICS_SOURCE.PICARD: {
            METRICS_OUTPUT_TYPE.DATA: [ "gc_bias.summary_metrics", "gc_bias.detail_metrics" ],
            METRICS_OUTPUT_TYPE.REPORT: "gc_bias.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: None,
        },
        METRICS_SOURCE.PARABRICKS: {
            METRICS_OUTPUT_TYPE.DATA: [ "gcbias_summary.txt", "gcbias_detail.txt" ],
            METRICS_OUTPUT_TYPE.REPORT: "gcbias.pdf",
            METRICS_OUTPUT_TYPE.IMAGE: "gcbias_0.png",
        }
    },
    METRIC_TOOL.COLLECT_SEQUENCING_ARTIFACT_METRICS: {
        METRICS_SOURCE.PICARD: {
            METRICS_OUTPUT_TYPE.DATA: [
                "bait_bias_detail_metrics",
                "bait_bias_summary_metrics",
                "pre_adapter_detail_metrics",
                "pre_adapter_summary_metrics",
                "error_summary_metrics",
            ],
            METRICS_OUTPUT_TYPE.REPORT: None,
            METRICS_OUTPUT_TYPE.IMAGE: None,
        },
        METRICS_SOURCE.PARABRICKS: {
            METRICS_OUTPUT_TYPE.DATA: [
                "sequencingArtifact.bait_bias_detail_metrics.txt",
                "sequencingArtifact.bait_bias_summary_metrics.txt",
                "sequencingArtifact.pre_adapter_detail_metrics.txt",
                "sequencingArtifact.pre_adapter_summary_metrics.txt",
                "sequencingArtifact.error_summary_metrics.txt",
            ],
            METRICS_OUTPUT_TYPE.REPORT: None,
            METRICS_OUTPUT_TYPE.IMAGE: None,
        }
    },
    METRIC_TOOL.COLLECT_QUALITY_YIELD_METRICS: {
        METRICS_SOURCE.PICARD: {
            METRICS_OUTPUT_TYPE.DATA: "quality_yield_metrics",
            METRICS_OUTPUT_TYPE.REPORT: None,
            METRICS_OUTPUT_TYPE.IMAGE: None,
        },
        METRICS_SOURCE.PARABRICKS: {
            METRICS_OUTPUT_TYPE.DATA: None,
            METRICS_OUTPUT_TYPE.REPORT: None,
            METRICS_OUTPUT_TYPE.IMAGE: None,
        }
    },
}


def get_multiple_metrics_output(
    source: METRICS_SOURCE,
    types: typing.List[METRICS_OUTPUT_TYPE] = [METRICS_OUTPUT_TYPE.DATA],
    tools: typing.Set[METRIC_TOOL] = set([
        METRIC_TOOL.COLLECT_ALIGNMENT_SUMMARY_METRICS,
        METRIC_TOOL.COLLECT_INSERT_SIZE_METRICS,
        METRIC_TOOL.QUALITY_SCORE_DISTRIBUTION,
        METRIC_TOOL.MEAN_QUALITY_BY_CYCLE,
        METRIC_TOOL.COLLECT_BASE_DISTRIBUTION_BY_CYCLE,
        METRIC_TOOL.COLLECT_GC_BIAS_METRICS,
        METRIC_TOOL.COLLECT_SEQUENCING_ARTIFACT_METRICS,
    ])
) -> typing.Dict[str, str]:
    """Get multiple_metrics output files, respecting files already present from other tools.

    Parameters
    ----------
    types: List[METRICS_OUTPUT_TYPE]
           Output types to generate output for.
    source: METRIC_SOURCE
          Tool to execute to get the data, i.e. in a snakemake-sense: get the prefixes from
          (see the corresponding wrappers)
    tools: set[METRIC_TOOL]
           Tools to generate output file stanzas for (others will be ignored). Includes all
           but collect_quality_yield_metrics by default, since the latter is not supported
           by parabricks afaik.

    Returns
    -------
    dict
        Dictionary mapping valid wrapper prefix values to an output file, respecting existing
        files.
    """

    output = {}
    for tool in tools:
        if source in METRICS_FILES[tool].keys():
            for type in types:
                if type in METRICS_FILES[tool][source].keys():
                    val = METRICS_FILES[tool][source][type]
                    if val is None:
                        continue
                    elif isinstance(val, str):
                        suffix = "." + type.value if not METRICS_FILES[tool][METRICS_SOURCE.PARABRICKS][type].endswith(type.value) else ""
                        output[val] = "output/{{sample}}/{{sample}}.{metric}{suffix}".format(
                            metric=METRICS_FILES[tool][METRICS_SOURCE.PARABRICKS][type],
                            suffix=suffix,
                        )
                    elif isinstance(val, list):
                        for i, x in enumerate(val):
                            suffix = "." + type.value if not METRICS_FILES[tool][METRICS_SOURCE.PARABRICKS][type][i].endswith(type.value) else ""
                            output[x] = "output/{{sample}}/{{sample}}.{metric}{suffix}".format(
                            metric=METRICS_FILES[tool][METRICS_SOURCE.PARABRICKS][type][i],
                            suffix=suffix,
                        )

    return output

def detect_multiple_metrics_type(
    files: typing.Iterable[str],
    default: METRICS_SOURCE = METRICS_SOURCE.PICARD,
):
    """Detect the metric source for a list of file paths.
    
    Parameters
    ----------
    files : Iterable[str]
            List of file paths.
    default : METRICS_SOURCE
              The default metrics source, returned when no
              match is found.

    Returns
    -------
    METRICS_SOURCE
        The metrics source used to generate these files, or
        default if not file is a multiple metrics file.
    """
    names_by_source = {
        METRICS_SOURCE.PARABRICKS: set(),
        METRICS_SOURCE.PICARD: set(),
    }
    for tool, group in METRICS_FILES.items():
        for source, filetypes in group.items():
            for type, name in filetypes.items():
                if name is not None:
                    if isinstance(name, str):
                        names_by_source[source].add(name)
                    elif isinstance(name, list):
                        names_by_source[source].update(name)
    
    files_by_source = {
        METRICS_SOURCE.PARABRICKS: 0,
        METRICS_SOURCE.PICARD: 0,
    }
    for path in files:
        if any(x in path for x in names_by_source[METRICS_SOURCE.PARABRICKS]):
            files_by_source[METRICS_SOURCE.PARABRICKS] += 1
        elif any(x in path for x in names_by_source[METRICS_SOURCE.PICARD]):
            files_by_source[METRICS_SOURCE.PICARD] += 1
    
    if files_by_source[METRICS_SOURCE.PARABRICKS] > files_by_source[METRICS_SOURCE.PICARD]:
        return METRICS_SOURCE.PARABRICKS
    elif files_by_source[METRICS_SOURCE.PICARD] > files_by_source[METRICS_SOURCE.PARABRICKS]:
        return METRICS_SOURCE.PICARD
    return default