##
# Snakemake wrapper script for picard CollectMultipleMetrics
# 
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
from uuid import uuid4
from snakemake.logging import logger
import os

shell.executable('/bin/bash')

TOOLS = set([
    "CollectAlignmentSummaryMetrics",
    "CollectInsertSizeMetrics",
    "QualityScoreDistribution",
    "MeanQualityByCycle",
    "CollectBaseDistributionByCycle",
    "CollectGcBiasMetrics",
    "RnaSeqMetrics",
    "CollectSequencingArtifactMetrics",
    "CollectQualityYieldMetrics"
])

FILE_SUFFIXES = set([
    # CollectAlignmentSummaryMetrics
    "alignment_summary_metrics",
    "read_length_histogram.pdf"
    # CollectInsertSizeMetrics
    "insert_size_metrics",
    "insert_size_histogram.pdf",
    # QualityScoreDistribution
    "quality_distribution_metrics",
    "quality_distribution.pdf",
    # MeanQualityByCycle
    "quality_by_cycle_metrics",
    "quality_by_cycle.pdf",
    # CollectBaseDistributionByCycle
    "base_distribution_by_cycle_metrics"
    "base_distribution_by_cycle.pdf"
    # CollectGcBiasMetrics
    "gc_bias.summary_metrics",
    "gc_bias.pdf",
    "gc_bias.detail_metrics",
    # CollectSequencingArtifactMetrics
    "bait_bias_detail_metrics",
    "bait_bias_summary_metrics",
    "pre_adapter_detail_metrics",
    "pre_adapter_summary_metrics",
    "error_summary_metrics",
    # CollectQualityYieldMetrics
    "quality_yield_metrics",
])

bam = snakemake.input.get("bam", None)
if bam is None:
    raise Exception("Please provide an input bam file.")
bam = os.path.abspath(bam)

output = snakemake.output[0]

reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("Reference must be set.")
reference = os.path.abspath(reference)

exclude_metrics = snakemake.params.get(
    "exclude_metrics",
    ["RnaSeqMetrics"]
)
include_metrics = snakemake.params.get(
    "include_metrics",
    None
)

if include_metrics is not None:
    metrics = set(include_metrics) & TOOLS
else:
    metrics = TOOLS - set(exclude_metrics)
metrics_str = " --PROGRAM " + " --PROGRAM ".join(metrics)

extra = snakemake.params.get("extra", "")
jvm_args = snakemake.params.get("jvm_args", "")

memory = snakemake.resources.get("mem_mb", None)

if 'Xmx' not in jvm_args and (memory is not None):
    jvm_args += " -Xmx{}m".format(memory)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    prefix = os.path.join(tempdir, str(uuid4()))
    shell(dedent("""
        set -ex
        picard {jvm_args} CollectMultipleMetrics {extra} \
            {metrics_str} \
            --INPUT {bam} \
            --OUTPUT {prefix} \
            --REFERENCE_SEQUENCE {reference} \
            --TMP_DIR {tempdir} {log}
    """))

    for key, value in snakemake.output.items():
        if key not in FILE_SUFFIXES:
            logger.warning("Unknown file suffix " + str(key))
        shell(dedent("""
            set -ex
            mv -v {prefix}.{key} {value}
        """))
