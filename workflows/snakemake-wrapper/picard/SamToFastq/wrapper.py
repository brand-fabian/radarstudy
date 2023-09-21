__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
from uuid import uuid4
import os

shell.executable('/bin/bash')

bam = snakemake.input.get('bam', None)
if bam is None:
    raise Exception("No input provided")
bam = os.path.abspath(bam)

ref = snakemake.params.get('reference', None)
if ref is None:
    raise Exception("No reference genome found.")
ref = os.path.abspath(ref)

fq_out = snakemake.output.get('fastqs', None)
if fq_out is None:
    raise Exception("Please provide an output fq path")
if len(fq_out) != 2:
    raise Exception("Please specify exactly two fq output paths")

threads = snakemake.threads
extra = snakemake.params.get("extra", "")
jvm_args = snakemake.params.get("jvm_args", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

memory = snakemake.resources.get("mem_mb", None)

if 'Xmx' not in jvm_args and (memory is not None):
    jvm_args += " -Xmx{}m".format(memory)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    first_read = os.path.abspath(fq_out[0])
    paired_read = os.path.abspath(fq_out[1])

    shell(dedent("""
        set -ex

        export TMPDIR={tempdir}

        picard {jvm_args} SamToFastq {extra} \
            --INPUT {bam} \
            --FASTQ {first_read} \
            --SECOND_END_FASTQ {paired_read} {log}
    """))
