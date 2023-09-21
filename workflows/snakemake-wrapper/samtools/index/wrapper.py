__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
import os

shell.executable('/bin/bash')

input_bam = snakemake.input.get("bam", None)
if input_bam is None:
    raise Exception("Error: Missing .bam")

index_out = snakemake.output.get("bai", None)
if index_out is None:
    raise Exception("Error: Missing output .bai")

extra = snakemake.params.get("extra", "")
threads = snakemake.threads
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    shell(dedent("""
        set -ex
        export TMPDIR={tempdir}

        samtools index -@{threads} {extra} {input_bam} {index_out} {log}
    """))
