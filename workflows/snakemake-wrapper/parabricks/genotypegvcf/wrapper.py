__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
import os

shell.executable('/bin/bash')

gvcf = snakemake.input.get('gvcf', None)
if gvcf is None:
    raise Exception("No input provided")
input_str = f"--in-gvcf {gvcf}"

ref = snakemake.params.get('reference', None)
if ref is None:
    raise Exception("No reference genome found.")
ref_str = f"--ref {ref}"

vcf_out = snakemake.output.get('vcf', None)
if vcf_out is None:
    raise Exception('Please provide a vcf output path.')
compress = vcf_out.endswith('.gz')

uncompressed_vcf = '.'.join(vcf_out.split('.')[:-1])
if compress:
    output_str = f"--out-vcf {uncompressed_vcf}"
else:
    output_str = f"--out-vcf {vcf_out}"


extra = snakemake.params.get("extra", "")
module_name = snakemake.params.get("module_name", "parabricks")
threads = snakemake.threads
threads_str = f"--num-threads {snakemake.threads}"
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(dedent("""
    # Setup scratch dir (local ssd/generic temp)
    SSD=${{SCRATCH_DIR:-}}
    if [ -n "$SSD" ]; then
        export TMPDIR=$SCRATCH_DIR
    fi
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    # Load parabricks if a module system is present
    if `type module 2>&1 | grep -q function`; then
        module load {module_name}
    fi

    pbrun genotypegvcf {extra} \\
        {ref_str} \\
        {input_str} \\
        {output_str} \\
        {threads_str} {log}
"""))

if compress:
    shell(dedent("""
        # Load parabricks if a module system is present
        if `type module 2>&1 | grep -q function`; then
            module load {module_name}
        fi

        pbgzip -n {threads} {uncompressed_vcf}
        pbrun indexgvcf --input {vcf_out}
        rm -rf {uncompressed_vcf}
    """))