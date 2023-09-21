__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
import os

bam_input = snakemake.input.get("bam", None)
if bam_input is None:
    raise Exception("missing bam input file")
bam_input = os.path.abspath(bam_input)

vcf_input = snakemake.input.get("vcf", None)
if vcf_input is None:
    raise Exception("missing vcf input file")
vcf_input = os.path.abspath(vcf_input)

ped = snakemake.input.get("ped", None)
if ped is None:
    raise Exception("no pedigree file provided")
ped = os.path.abspath(ped)

reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("no reference provided")
reference = os.path.abspath(reference)

vcf_output = snakemake.output.get("vcf", None)
if vcf_output is None:
    raise Exception("no vcf output provided")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory(dir=os.getenv('SCRATCH_DIR', gettempdir())) as tempdir:
    vcf_default = os.path.join(tempdir, "phased.vcf")
    shell(dedent("""
        set -x
        export TMPDIR={tempdir}

        whatshap phase \
            --ped {ped} \
            --reference {reference} \
            -o {vcf_default} \
            {vcf_input} {bam_input} {log}
        bcftools view {vcf_default} | bgzip > {vcf_output}
        tabix {vcf_output}
    """))