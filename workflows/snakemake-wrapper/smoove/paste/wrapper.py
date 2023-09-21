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

name = snakemake.params.get("name", "cohort")

output_files = {
    "vcf": "{}.smoove.square.vcf.gz".format(name)
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

        smoove paste --outdir {tempdir} \
            {extra} --name {name} \
            {vcf_str} {log}
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