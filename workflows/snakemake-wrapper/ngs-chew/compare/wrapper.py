##
# Snakemake wrapper script for ngs-chew compare
#
__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os
import sys
from tempfile import TemporaryDirectory, gettempdir

shell.executable('/bin/bash')

extra = snakemake.params.get("extra", "")

sample_name = snakemake.input[0].split("/")[-1].split(".")[0]
title_aab = snakemake.params.get("title_aab", "Alternate allele balance")
title_var_het = snakemake.params.get("title_var_het", "Variance of het alternate allele balance")

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

input_npz = snakemake.input.get("npz", None)
input_vcf = snakemake.input.get("vcf", None)
if input_npz is None:
    raise Exception("Error: No npz fingerprints found.")
if input_vcf is None:
    raise Exception("Error: No vcf fingerprints found.")

output_relatedness = snakemake.output.get("relatedness", None)
if output_relatedness is None:
    raise Exception("Error: No output file for relatedness plot specified.")

with TemporaryDirectory(dir=os.getenv('SCRATCH_DIR', gettempdir())) as tempdir:
    print("Using temp dir {}".format(tempdir), file=sys.stderr)
    output_aab = snakemake.output.get("aab", None)
    if output_aab is not None:
        shell(r"""
            ngs-chew plot_aab \
                --title "{title_aab}" \
                --aab-cache "{tempdir}/aab_cache.json" \
                "{output_aab}" \
                {snakemake.input.vcf} {log}
        """)
    output_var_het = snakemake.output.get("var_het", None)
    if output_var_het is not None:
        shell(r"""
            ngs-chew plot_var_het \
                --title "{title_var_het}" \
                --var_het_cache "{tempdir}/var_het_cache.json" \
                "{output_var_het}" \
                {snakemake.input.vcf} {log}
        """)
    shell(r"""
        ngs-chew compare --output-prefix "{tempdir}/comparison" {snakemake.input.npz}
        ngs-chew plot_compare "{tempdir}/comparison.relatedness.txt" {snakemake.output.relatedness}
    """)