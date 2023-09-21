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

output = snakemake.output[0]

reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("Reference must be set.")
reference = os.path.abspath(reference)

extra = snakemake.params.get("extra", "")
jvm_args = snakemake.params.get("jvm_args", "")

memory = snakemake.resources.get("mem_mb", None)

if 'Xmx' not in jvm_args and (memory is not None):
    jvm_args += " -Xmx{}m".format(memory)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    shell(dedent("""
        set -ex
        picard {jvm_args} CollectWgsMetrics {extra} \
            --INPUT {bam} \
            --OUTPUT {output} \
            --REFERENCE_SEQUENCE {reference} \
            --TMP_DIR {tempdir} {log}
    """))
