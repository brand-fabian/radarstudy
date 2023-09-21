__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

input_vcf = snakemake.input.get("vcf", None)
if input_vcf is None:
    raise Exception("Error: Please provide an input vcf path")
if not os.path.isfile(input_vcf):
    raise Exception(f"{input_vcf} could not be found.")

if input_vcf.endswith('gz'):
    input_str = f"--gzvcf {input_vcf}"
else:
    input_str = f"--vcf {input_vcf}"

output_vcf = snakemake.output.get("vcf", None)
if output_vcf is None:
    raise Exception("Error: Please provide an output vcf path")

text_output = 'vcf' not in output_vcf

if text_output:
    bgzip = '-c'
    tabix = ''
else:
    if output_vcf.endswith('gz'):
        bgzip = '--recode -c | bgzip -c'
        tabix = f"&& tabix {output_vcf}"
    elif output_vcf.endswith('bcf'):
        bgzip = '--recode-bcf -c'
        tabix = ''
    else:
        bgzip = '--recode -c'
        tabix = ''

extra = snakemake.params.get("extra", "")

shell("""
    exec 2> >(tee -a "{snakemake.log}")
    vcftools {input_str} {extra} {bgzip} > {output_vcf} {tabix}
""")