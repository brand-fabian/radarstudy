##
# Snakemake wrapper script for picard CollectMultipleMetrics
# 
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
import os

shell.executable('/bin/bash')

bam = snakemake.input.get("bam", None)
if bam is None:
    raise Exception("Please provide an input bam file.")
bam = os.path.abspath(bam)

output = os.path.abspath(snakemake.output[0])

reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("Reference must be set.")
reference = os.path.abspath(reference)

extra = snakemake.params.get("extra", "")
module_name = snakemake.params.get("module_name", "parabricks")
threads = snakemake.threads
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    shell(dedent("""
        set -ex

        # Load parabricks if a module system is present
        if `type module 2>&1 | grep -q function`; then
            module load {module_name}
        fi

        pbrun bammetrics {extra} \
            --bam {bam} \
            --ref {reference} \
            --out-metrics-file {output} \
            --num-threads {threads} \
            --tmp-dir {tempdir} {log}
    """))
