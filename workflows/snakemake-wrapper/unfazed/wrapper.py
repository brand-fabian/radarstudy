__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
import os


ped = None
vcf_input = None
sites_input = None
bam_input = []
# If there is a sample named "ped", "vcf" or "sites",
# we may get a problem here...
for k, v in snakemake.input.items():
    if k == 'ped':
        ped = os.path.abspath(v)
    elif k == 'vcf':
        vcf_input = os.path.abspath(v)
    elif k == 'sites':
        sites_input = os.path.abspath(v)
    else:
        bam_input.append("{}:{}".format(
            k, os.path.abspath(v),
        ))

if any(x is None for x in [ped, vcf_input, sites_input]):
    raise Exception("missing input data for unfazed")

if len(bam_input) == 0:
    raise Exception("no bams provided for unfazed")

bam_input = "--bam-pairs " + " ".join(bam_input)

reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("no reference provided")
reference = os.path.abspath(reference)

reference_version = snakemake.params.get("reference_version", "37")
extra = snakemake.params.get("extra", "")
threads = snakemake.threads

vcf_output = snakemake.output.get("vcf", None)
if vcf_output is None:
    raise Exception("no vcf output provided")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with TemporaryDirectory(dir=os.getenv('SCRATCH_DIR', gettempdir())) as tempdir:
    vcf_default = os.path.join(tempdir, "phased.vcf")
    shell(dedent("""
        set -x
        export TMPDIR={tempdir}

        unfazed --verbose -t {threads} \
            --outfile {vcf_default} \
            --dnms {vcf_input} \
            --sites {sites_input} \
            --ped {ped} \
            --reference {reference} --build {reference_version} \
            {bam_input} {log}
        bcftools view {vcf_default} | bgzip > {vcf_output}
        tabix {vcf_output}
    """))
