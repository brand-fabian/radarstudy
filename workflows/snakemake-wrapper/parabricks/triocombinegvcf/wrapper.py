__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
import os

shell.executable('/bin/bash')

gvcfs = snakemake.input.get('gvcfs', None)
if gvcfs is None:
    raise Exception("No input provided")
input_str = "--in-gvcf " + " --in-gvcf ".join(gvcfs)

ref = snakemake.params.get('reference', None)
if ref is None:
    raise Exception("No reference genome found.")
ref_str = f"--ref {ref}"

gvcf_out = snakemake.output.get('gvcf', None)
if gvcf_out is None:
    raise Exception('Please provide a vcf output path.')
output_str = f"--out-variants {gvcf_out}"

module_name = snakemake.params.get("module_name", "parabricks")
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

    pbrun triocombinegvcf \\
        {ref_str} \\
        {input_str} \\
        {output_str} {log}
    pbrun indexgvcf --input {gvcf_out}
"""))