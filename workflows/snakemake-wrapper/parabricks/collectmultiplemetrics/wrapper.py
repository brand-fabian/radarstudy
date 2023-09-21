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

FILE_SUFFIXES = set([
    # CollectAlignmentSummaryMetrics
    "alignment.txt",
    "read_length_histogram.pdf"
    # CollectInsertSizeMetrics
    "insert_size.txt",
    "insert_size.pdf",
    "insert_size.png",
    # QualityScoreDistribution
    "qualityscore.pdf",
    "qualityscore.png",
    "qualityscore.txt",
    # MeanQualityByCycle
    "mean_quality_by_cycle.pdf",
    "mean_quality_by_cycle.png",
    "mean_quality_by_cycle.txt",
    # CollectBaseDistributionByCycle
    "base_distribution_by_cycle.txt",
    "base_distribution_by_cycle.png",
    "base_distribution_by_cycle.pdf",
    # CollectGcBiasMetrics
    "gcbias_0.png",
    "gcbias_summary.txt",
    "gcbias.pdf",
    "gcbias_detail.txt",
    # CollectSequencingArtifactMetrics
    "sequencingArtifact.bait_bias_detail_metrics.txt",
    "sequencingArtifact.bait_bias_summary_metrics.txt",
    "sequencingArtifact.pre_adapter_detail_metrics.txt",
    "sequencingArtifact.pre_adapter_summary_metrics.txt",
    "sequencingArtifact.error_summary_metrics.txt",
])

bam = snakemake.input.get("bam", None)
if bam is None:
    raise Exception("Please provide an input bam file.")
bam = os.path.abspath(bam)

reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("Reference must be set.")
reference = os.path.abspath(reference)


extra = snakemake.params.get("extra", "")
module_name = snakemake.params.get("module_name", "parabricks")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    prefix = os.path.join(tempdir, str(uuid4()))
    os.makedirs(prefix)
    shell(dedent("""
        set -ex

        # Load parabricks if a module system is present
        if `type module 2>&1 | grep -q function`; then
            module load {module_name}
        fi

        pbrun collectmultiplemetrics {extra} \
            --gen-all-metrics \
            --bam {bam} \
            --out-qc-metrics-dir {prefix} \
            --ref {reference} \
            --tmp-dir {tempdir} {log}
    """))

    for key, value in snakemake.output.items():
        if key not in FILE_SUFFIXES:
            logger.warning("Unknown file suffix " + str(key))
        shell(dedent("""
            set -ex
            mv -v {prefix}/{key} {value}
        """))
