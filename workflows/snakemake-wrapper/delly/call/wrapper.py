__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
from textwrap import dedent

input_bam = snakemake.input.get("bam", None)
if input_bam is None:
    raise Exception("Please provide an input bam file.")
input_bam = os.path.abspath(input_bam)
if not os.path.isfile(input_bam):
    raise Exception("Could not find input bam at {}".format(input_bam))

input_sites = snakemake.input.get("sites", None)
if input_sites is not None:
    input_sites_str = '-v {}'.format(input_sites)
else:
    input_sites_str = ''

reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("Please provide the path to a human genome reference as parameter.")
reference = os.path.abspath(reference)
if not os.path.isfile(reference):
    raise Exception("Could not find reference file at {}".format(reference))

ref_type = snakemake.params.get("reference_type", "hg19")
if ref_type not in [ "hg19", "hg38" ]:
    raise Exception("Please provide a reference_type. Valid choices: [ hg19, hg38 ].")

output = snakemake.output[0]

exclude_regions = snakemake.params.get("exclude_regions", True)
exclude_file = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    '..',
    'data',
    '{}.excl.tsv'.format(ref_type)
))
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")
exclude_str = '-x {}'.format(exclude_file) if exclude_regions else ''
module_name = snakemake.params.get("module_name", "delly")
threads = snakemake.threads if snakemake.threads > 0 else 1

if not os.path.isdir(os.path.dirname(output)):
    os.makedirs(os.path.dirname(output), exist_ok=True)


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
    export OMP_NUM_THREADS={threads}
    delly call {extra} -g {reference} \
        -o {output} {exclude_str} {input_sites_str} \
        {input_bam} {log}
"""))