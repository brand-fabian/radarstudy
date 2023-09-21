__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
from textwrap import dedent
from tempfile import gettempdir, TemporaryDirectory

shell.executable('/bin/bash')

input_vcf = snakemake.input.get("vcf", None)
if input_vcf is None:
    raise Exception("Missing input file")
input_vcf = os.path.abspath(input_vcf)
input_bams = snakemake.input.get("bam", None)
if input_bams is None:
    raise Exception("Missing bam files")
fam = snakemake.input.get("fam", None)
if fam is not None:
    fam_str = f"--ped {os.path.abspath(fam)}"
else:
    fam_str = ""
bam_str =  "-b " + " ".join(map(os.path.abspath, input_bams))
bam_samples = snakemake.params.get("bam_samples", "")
if isinstance(bam_samples, list):
    bam_samples = " ".join(bam_samples)
log = snakemake.log[0]
if log is None or log == "":
    log_str = ""
else:
    log_str = f">> {os.path.abspath(log)} 2>&1"
extra = snakemake.params.get("extra", "")
output = os.path.abspath(snakemake.output[0])

with TemporaryDirectory(
    dir=os.getenv("SCRATCH_DIR", gettempdir())
) as tempdir:
    shell(dedent("""
        set -ex
        cd {tempdir}
        samplot vcf {extra} \
            --vcf {input_vcf} \
            {bam_str} \
            --sample_ids {bam_samples} \
            -O png \
            -d {tempdir} \
            {fam_str} {log_str}
        tar -czvf {output} ./*
        cd -
    """))
