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

bam = snakemake.input.get("bam", None)
if bam is None:
    raise Exception("Please provide an input bam file.")
bam = os.path.abspath(bam)

output_bam = snakemake.output.get("bam", None)
if output_bam is None:
    raise Exception("Please provide an output bam path.")
output_bam = os.path.abspath(output_bam)

rgid = snakemake.params.get("rgid", None)
rglb = snakemake.params.get("rglb", None)
rgpl = snakemake.params.get("rbpl", "Illumina")
rgpu = snakemake.params.get("rgpu", None)
rgsm = snakemake.params.get("rgsm", None)
if any(x is None for x in [rgid, rglb, rgpl, rgpu, rgsm]):
    raise Exception("Please provide read group information in params.")

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
        {{
            picard {jvm_args} AddOrReplaceReadGroups {extra} \
                --INPUT {bam} \
                --OUTPUT {output_bam} \
                --RGID {rgid} \
                --RGLB {rglb} \
                --RGPL {rgpl} \
                --RGPU {rgpu} \
                --RGSM {rgsm} \
                --TMP_DIR {tempdir}
            samtools index {output_bam}
        }} {log}
    """))
