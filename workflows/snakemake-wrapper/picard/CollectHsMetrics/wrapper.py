##
# Snakemake wrapper script for picard CollectHsMetrics
# 
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")
jvm_args = snakemake.params.get("jvm_args", "")

memory = snakemake.resources.get("mem_mb", None)

if 'Xmx' not in jvm_args and (memory is not None):
    jvm_args += " -Xmx{}m".format(memory)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

per_base = snakemake.output.get("per_base", "")
per_target = snakemake.output.get("per_target", "")

bait = snakemake.input.get("bait_intervals", snakemake.input.target_intervals)

if snakemake.input.get("target_intervals", "") == "":
    raise Exception("CollectHsMetrics at least needs a target intervals.")

shell(r"""
VALUE=${{SCRATCH_DIR:-}}

if [ -z "$VALUE" ]
then
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT
    echo "Using temp dir $TMPDIR."
else
    export TMPDIR=$SCRATCH_DIR/gatk
    mkdir -p "$TMPDIR"
    echo "Using temp dir $TMPDIR."
fi

PER_BASE="{per_base}"
if [ -z "$PER_BASE" ]
then
    PER_BASE="$TMPDIR/per_base_coverage.txt"
fi
PER_TARGET="{per_target}"
if [ -z "$PER_TARGET" ]
then
    PER_TARGET="$TMPDIR/per_target_coverage.txt"
fi

picard {jvm_args} CollectHsMetrics {extra} \
    INPUT={snakemake.input.bam} \
    OUTPUT={snakemake.output.metrics} \
    PER_BASE_COVERAGE="$PER_BASE" \
    PER_TARGET_COVERAGE="$PER_TARGET" \
    REFERENCE_SEQUENCE={snakemake.input.ref} \
    BAIT_INTERVALS={bait} \
    TARGET_INTERVALS={snakemake.input.target_intervals} {log}
""")
