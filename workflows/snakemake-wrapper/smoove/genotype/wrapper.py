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

input_bam = snakemake.input.get("bam", None)
if input_bam is None:
    raise Exception("No bam file specified.")

sample = snakemake.params.get("sample", input_bam.split('.')[0])
reference = snakemake.params.get("reference", None)
if reference is None:
    raise Exception("Please provide a path to a reference genome.")
reference = os.path.abspath(reference)

output_files = {
    "vcf": "{}-smoove.genotyped.vcf.gz".format(sample)
}

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log_append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
threads = snakemake.threads if snakemake.threads > 0 else 1
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

        smoove genotype --outdir {tempdir} \
            {extra} -p {threads} --name {sample} \
            --vcf {input_vcf} \
            --fasta {reference} \
            {input_bam} {log}
    """))

    for key, value in output_files.items():
        out = snakemake.output.get(key, None)
        if out is not None:
            shell(dedent("""
                set -ex
                if [ ! -d $(dirname "{out}") ]; then
                    mkdir $(dirname "{out}")
                fi
                mv -v {tempdir}/{value} {out} {log_append}
            """))

            if any(out.endswith(e) for e in [ 'vcf.gz', 'bcf', 'vcf' ]):
                shell(dedent("""
                    tabix {out}
                """))