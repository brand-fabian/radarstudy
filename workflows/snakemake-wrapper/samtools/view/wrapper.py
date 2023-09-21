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

bam_out = snakemake.output.get("bam", None)
if bam_out is None:
    raise Exception("Error: missing output bam")
bai_out = snakemake.output.get("bai", None)

bam_fmt = ""
if bam_out.endswith("bam"):
    bam_fmt = "--bam"
elif bam_out.endswith("cram"):
    bam_fmt = "--cram"

reference = snakemake.params.get("reference", None)
reference_str = ""
if reference is not None:
    reference_str = "--reference {}".format(reference)

extra = snakemake.params.get("extra", "")
threads = snakemake.threads
log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
    shell(dedent("""
        set -ex
        export TMPDIR={tempdir}

        samtools view {extra} -@ {threads} {reference_str} {bam_fmt} -o {bam_out} {input_bam} {log}
    """))

    if bai_out is not None:
        shell(dedent("""
            set -ex
            samtools index {bam_out}
        """))

        if bam_out.endswith("bam") and bai_out != "{}.bai".format(bam_out):
            shell(dedent("""
                mv -v {bam_out}.bai {bai_out}
            """))
        elif bam_out.endswith("cram") and bai_out != "{}.csi".format(bai_out):
            shell(dedent("""
                mv -v {bam_out}.csi {bai_out}
            """))