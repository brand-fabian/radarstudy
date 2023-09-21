__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
from textwrap import dedent
from collections.abc import Iterable
from tempfile import TemporaryDirectory, gettempdir

input_vcf = snakemake.input.get("vcf", None)
if input_vcf is None:
    raise Exception("No vcf file specified.")
if isinstance(input_vcf, str):
    input_vcf = [ input_vcf ]
if not isinstance(input_vcf, Iterable):
    raise Exception("Please supply an iterable of vcf files as input.")
input_vcf = list(map(os.path.abspath, input_vcf))
vcf_str = ' '.join(input_vcf)

sample = snakemake.params.get("sample", "merged")
reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("Please provide a path to a reference genome.")
reference = os.path.abspath(reference)

output_files = {
    "vcf": "{}.sites.vcf.gz".format(sample)
}

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log_append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
module_name = snakemake.params.get("module_name", "smoove")

with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR", gettempdir())) as tempdir:
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

        smoove merge --outdir {tempdir} \
            {extra} --name {sample} \
            --fasta {reference} {vcf_str} {log}
    """))

    for key, value in output_files.items():
        out = snakemake.output.get(key, None)
        if out is not None:
            if out.endswith('vcf.gz'):
                sort_out = f"-o {out} -O z"
            elif out.endswith('bcf'):
                sort_out = f"-o {out} -O b"
            elif out.endswith('vcf'):
                sort_out = f"-o {out} -O v"
            else:
                sort_out = f"-o {out} -O z"
            shell(dedent("""
                set -ex
                if [ ! -d $(dirname "{out}") ]; then
                    mkdir $(dirname "{out}")
                fi
                bcftools sort -T {tempdir} {sort_out} {tempdir}/{value} {log_append}
                tabix {out}
            """))