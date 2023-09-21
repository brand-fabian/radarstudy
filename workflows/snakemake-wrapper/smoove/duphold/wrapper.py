__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
from textwrap import dedent
from collections.abc import Iterable

input_vcf = snakemake.input.get("vcf", None)
if input_vcf is None:
    raise Exception("No vcf file specified.")

input_bam = snakemake.input.get("bam", None)
if input_bam is None:
    raise Exception("No bam file specified.")
if isinstance(input_bam, str):
    input_bam = [ input_bam ]
if not isinstance(input_bam, Iterable):
    raise Exception("Please supply an iterable of bam files as input.")
input_bam = list(map(os.path.abspath, input_bam))
bam_str = ' '.join(input_bam)

reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("Please provide a path to a reference genome.")
reference = os.path.abspath(reference)

output = snakemake.output['vcf']

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log_append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
threads = snakemake.threads if snakemake.threads > 0 else 1
module_name = snakemake.params.get("module_name", "smoove")

if output.endswith('gz'):
    compress = True
    vcf_out = '.'.join(output.split('.')[:-1])
else:
    compress = False
    vcf_out = output

shell(dedent("""
    set -ex
    if `type module 2>&1 | grep -q function`; then
        module load {module_name}
        if [ -d "/ssd" ]; then
            export SINGULARITY_BINDPATH="/ceph01,/gpfs01,/ssd,/tmp"
        else
            export SINGULARITY_BINDPATH="/ceph01,/gpfs01,/tmp"
        fi
    fi

    smoove duphold {extra} -p {threads} \
        --fasta {reference} \
        --vcf {input_vcf} \
        --outvcf {vcf_out} \
        {bam_str} {log}
"""))

if compress:
    shell(dedent("""
        set -ex
        bgzip {vcf_out}
        tabix {output}
    """))