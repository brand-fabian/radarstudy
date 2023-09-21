__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
from textwrap import dedent
from tempfile import TemporaryDirectory, gettempdir
from uuid import uuid4
import os

shell.executable('/bin/bash')

bam = snakemake.input.get('bam', None)
if bam is None:
    raise Exception("No input provided")
input_str = f"--in-bam {bam}"

ref = snakemake.params.get('reference', None)
if ref is None:
    raise Exception("No reference genome found.")
ref_str = f"--ref {ref}"

sample = snakemake.params.get('sample', None)
if sample is None:
    raise Exception("Please provide a sample name.")

vcf_out = snakemake.output.get('vcf', None)
gvcf_out = snakemake.output.get('gvcf', None)
if vcf_out is None and gvcf_out is None:
    raise Exception('Please provide either a vcf or gvcf output path.')
if vcf_out is not None:
    output_str = ""
elif gvcf_out is not None:
    output_str = "--gvcf"
output_gz = gvcf_out.endswith("gz") if gvcf_out is not None else vcf_out.endswith("gz")


recal_in = snakemake.params.get('recal_in', None)
if recal_in is not None:
    recal_str = f"--in-recal-file {recal_in}"
else:
    recal_str = ""


intervals = snakemake.params.get('intervals', [])
if len(intervals) > 0:
    intervals_str = "--interval-file " + " --interval-file ".join(intervals)
else:
    intervals_str = ""


extra = snakemake.params.get("extra", "")
module_name = snakemake.params.get("module_name", "parabricks")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


with TemporaryDirectory(dir=os.getenv("SCRATCH_DIR",  gettempdir())) as tempdir:
    ending = "g.vcf" if gvcf_out is not None else "vcf"
    ending += ".gz" if output_gz else ""
    prefix = os.path.join(tempdir, f"{str(uuid4())}.{ending}")
    index = "tbi" if output_gz else "idx"
    output_file = os.path.abspath(gvcf_out) if gvcf_out is not None else os.path.abspath(vcf_out)
    shell(dedent("""
        set -ex
        export TMPDIR="{tempdir}"

        # Load parabricks if a module system is present
        if `type module 2>&1 | grep -q function`; then
            module load {module_name}
        fi

        pbrun haplotypecaller {extra} \\
            {ref_str} \\
            {input_str} \\
            {output_str} --out-variants {prefix} \\
            {recal_str} \\
            --tmp-dir {tempdir} \\
            {intervals_str} {log}

        mv -v {prefix} {output_file}
        mv -v {prefix}.{index} {output_file}.{index}
    """))