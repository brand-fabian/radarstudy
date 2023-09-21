__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

input_vcf = snakemake.input.get('vcf', None)
if input_vcf is None:
    raise Exception("Please provide an input vcf.")
filter_str = snakemake.params.get('filter_str', None)
if filter_str is None:
    raise Exception("No filter string provided to snpsift filter.")

output_vcf = snakemake.output.get('vcf', None)
if output_vcf is None:
    raise Exception("Please provide an output vcf path")

if input_vcf.endswith('gz'):
    input_str = f"zcat {input_vcf}"
else:
    input_str = f"cat {input_vcf}"

if output_vcf.endswith('gz'):
    output_str = f"| bgzip -c > {output_vcf} && tabix {output_vcf}"
else:
    output_str = f"> {output_vcf} && tabix {output_vcf}"

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell("""
    ({input_str} | SnpSift filter {extra} "{filter_str}" {output_str}) {log}
""")