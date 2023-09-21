__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

import os
from textwrap import dedent
from snakemake.shell import shell

shell.executable('/bin/bash')

in_vcf = snakemake.input[0]
if not isinstance(in_vcf, str):
    raise Exception("Error: please provide exactly one input vcf.")
in_vcf_str = os.path.abspath(in_vcf)

output_vcf = snakemake.output.get("vcf", None)
if output_vcf is None:
    raise Exception("Error: Please provide an output vcf file")

extra = snakemake.params.get("extra", "")

vcf_type_out = ""
if output_vcf.split(".")[-1] == "gz":
    vcf_type_out += " -O z"
elif output_vcf.split(".")[-1] == "bcf":
    vcf_type_out += " -O b"

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(dedent("""
    set -x
    {{
        bcftools view {vcf_type_out} -o {output_vcf} {extra} {in_vcf_str}
        tabix {output_vcf}
    }} {log}
"""))
