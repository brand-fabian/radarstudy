__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
import os

shell.executable('/bin/bash')

fqs = snakemake.input.get('fq', None)
if fqs is None:
    raise Exception("No input provided")
if len(snakemake.input.fq) % 2 != 0:
    raise Exception("Please provide paired end reads.")

ref = snakemake.params.get('reference', None)
if ref is None:
    raise Exception("No reference genome found.")

sample = snakemake.params.get('sample', None)
if sample is None:
    raise Exception("Please provide a sample name.")

bam_out = snakemake.output.get('bam', None)
if bam_out is None:
    raise Exception("Please provide an output bam path")

recal_out = snakemake.output.get('recal', '')
dups_metrics = snakemake.output.get('dups_metrics', '')

extra = snakemake.params.get("extra", "")
module_name = snakemake.params.get("module_name", "parabricks")
threads = max(1, min(snakemake.threads, 16))
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

input_str = ""
for fq1, fq2 in zip(fqs[::2], fqs[1::2]):
    input_str += f' --in-fq {fq1} {fq2}'


shell(dedent("""
    # Setup scratch dir (local ssd/generic temp)
    SSD=${{SCRATCH_DIR:-}}
    if [ -n "$SSD" ]; then
        export TMPDIR=$SCRATCH_DIR
    fi
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    RECAL_OUT="{recal_out}"
    if [ -z "$RECAL_OUT" ]; then
        RECAL_OUT=$(mktemp)
    fi
    DUPS_METRICS="{dups_metrics}"
    if [ -z "$DUPS_METRICS" ]; then
        DUPS_METRICS=$(mktemp)
    fi

    # Load parabricks if a module system is present
    if `type module 2>&1 | grep -q function`; then
        module load {module_name}
    fi

    pbrun fq2bam --ref {ref} \\
        {input_str} \\
        --out-bam {bam_out} \\
        --out-recal-file $RECAL_OUT \\
        --out-duplicate-metrics $DUPS_METRICS \\
        --tmp-dir $TMPDIR \\
        --num-cpu-threads {threads} \\
        --gpusort \\
        --gpuwrite \\
        --read-group-sm {sample} \\
        --read-group-lb {sample} \\
        --read-group-pl Illumina \\
        --read-group-id-prefix {sample} \\
        {extra} {log}
"""))